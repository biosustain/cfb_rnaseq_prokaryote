rule hisat2_build:
    input:
        'references/{reference}.fa'
    output:
      # extention may change and suffixes may be reduced based on input and paramters
      expand('hisat2/build_index/{{reference}}.{suffix}.ht2', suffix=[i for i in range(1,9)])
    log:
      'logs/hisat2/build_index/{reference}.log'
    params:
      prefix='hisat2/build_index/{reference}',
      extra=''
    threads: 4
    resources:
        walltime = 1200
    wrapper:
        f'{wrapper_source}/bio/hisat2/build_index'

rule hisat2_align: # not satisfied with the input structure of this wrapper.
    input:
        r1 = lambda wildcards: f'reads/{wildcards.read_type}/{wildcards.read_type}_{wildcards.group}_R1_paired.fastq.gz',
        r2 = lambda wildcards: f'reads/{wildcards.read_type}/{wildcards.read_type}_{wildcards.group}_R2_paired.fastq.gz',
        rU = lambda wildcards: f'reads/{wildcards.read_type}/{wildcards.read_type}_{wildcards.group}_R12_unpaired.fastq.gz' if wildcards.read_type == 'trimmed' else [],
        indexes = expand('hisat2/build_index/{{reference}}.{suffix}.ht2', suffix=[i for i in range(1,9)])
    output:
      'hisat2/align/{read_type}/{reference}/{group}/accepted_hits.sam',
      'hisat2/align/{read_type}/{reference}/{group}/align_summary.txt',
#       'hisat2/align/{read_type}/{reference}/{group}/{read_type}_{reference}_{group}.align_summary.txt'
    log:                                # optional
      'logs/hisat2/align/{read_type}/{reference}/{group}.log'
    params:                             # idx is required, extra is optional
      idx='hisat2/build_index/{reference}',
      extra=config['params']['hisat2']
    threads: 8                          # optional, defaults to 1
    resources:
        walltime = 7200 
    wrapper:
        f'{wrapper_source}/bio/hisat2/align'

rule sam_to_bam:
    input:
        'hisat2/align/{read_type}/{reference}/{group}/accepted_hits.sam'
    output:
        'hisat2/align/{read_type}/{reference}/{group}/accepted_hits.bam'
    log:
        'logs/samtools/{read_type}/{reference}/{group}_alignments_to_bam.log'
    threads: 2
    resources:
        walltime = 3600 
    params:
        '-hbu'
    wrapper:
        f'{wrapper_source}/bio/samtools/view'

rule sort_bam:
    input:
        'hisat2/align/{read_type}/{reference}/{group}/accepted_hits.bam'
    output:
        'hisat2/align/{read_type}/{reference}/{group}/accepted_hits_sorted.bam'
    log:
        'logs/samtools/{read_type}/{reference}/{group}_sort_bam.log'
    threads: 2
    resources:
        walltime = 3600 
    params:
        ''
    wrapper:
        f'{wrapper_source}/bio/samtools/sort'

rule index_bam:
    input:
        'hisat2/align/{read_type}/{reference}/{group}/accepted_hits_sorted.bam'
    output:
        'hisat2/align/{read_type}/{reference}/{group}/accepted_hits_sorted.bam.bai'
    log:
        'logs/samtools/{read_type}/{reference}/{group}_index_bam.log'
    threads: 1
    resources:
        walltime = 600 
    params:
        ''
    wrapper:
        f'{wrapper_source}/bio/samtools/index'

rule extract_mask_alignment:
    input:
        'hisat2/align/{read_type}/{reference}/{group}/accepted_hits_sorted.bam',
        'references/combined_mask.bed'
    output:
        'hisat2/align/{read_type}/{reference}/{group}/accepted_hits_sorted_mask.bam'
    log:
        'logs/samtools/{read_type}/{reference}/{group}_alignments_to_mask.log'
    threads: 1
    resources:
        walltime = 1800
    params:
        config['params']['samtools']['view']['mask']
    wrapper:
        f'{wrapper_source}/bio/samtools/view'

rule flagstat:
    input:
        'hisat2/align/{read_type}/{reference}/{group}/{alignment}.bam',
    output:
        'hisat2/align/{read_type}/{reference}/{group}/{alignment}.flagstat'
    log:
        'logs/samtools/{read_type}/{reference}/{group}_{alignment}.flagstat.log'
    threads: 1
    params:
        ''
    wrapper:
        f'{wrapper_source}/bio/samtools/flagstat'

rule featureCounts:
    input:
        reference = 'references/{reference}.gff3',
        alignment = 'hisat2/align/{read_type}/{reference}/{group}/accepted_hits_sorted.bam'
    output:
        output = 'hisat2/align/{read_type}/{reference}/{group}/featureCounts.out',
        summary = 'hisat2/align/{read_type}/{reference}/{group}/featureCounts.out.summary'
    threads: 1
    params:
        config['params']['featureCounts']
    log:
        'logs/feature_counts/{read_type}/{reference}/{group}/featureCounts.log'
    wrapper:
        f'{wrapper_source}/bio/subread/featureCounts'
        
rule featureCountsCollector:
    input:
        expand('hisat2/align/{{read_type}}/{{reference}}/{group}/featureCounts.{type}',
                               group=group_to_files.keys(), 
                               type=['out', 'out.summary'])
    output:
        counts = 'hisat2/align/{read_type}/{reference}/collected_counts.csv',
        summary = 'hisat2/align/{read_type}/{reference}/collected_counts.summary.csv'
    log:
        'logs/feature_counts/{read_type}/{reference}/featureCountsCollector.log'
    run:
        pd.concat([pd.read_table(f'hisat2/align/{wildcards.read_type}/{wildcards.reference}/{group}/featureCounts.out',
                                 header=1,
                                 skiprows=0,
                                 index_col=0,
                                 usecols=[0,6],
                                 names=['Feature', group]) \
                                                   for group in group_to_files.keys()], 
                                 axis=1).to_csv(output.counts)
        pd.concat([pd.read_table(f'hisat2/align/{wildcards.read_type}/{wildcards.reference}/{group}/featureCounts.out.summary',
                                 header=None,
                                 skiprows=1,
                                 index_col=0,
                                 names=['Status', group]) \
                                                  for group in group_to_files.keys()], 
                                 axis=1).to_csv(output.summary)
