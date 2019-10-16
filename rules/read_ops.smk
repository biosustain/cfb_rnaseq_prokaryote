rule concat_paired_reads:
    input:
        lambda wildcards: sample_to_files[wildcards.sample_id]
    params:
        space_sepped_R1_fastqs=lambda wildcards: ' '.join([read for read in sample_to_files[wildcards.sample_id] if read.find('_R1_')>-1]),
        space_sepped_R2_fastqs=lambda wildcards: ' '.join([read for read in sample_to_files[wildcards.sample_id] if read.find('_R2_')>-1]),
    output:
        'reads/raw/raw_{sample_id}_R1_paired.fastq.gz',
        'reads/raw/raw_{sample_id}_R2_paired.fastq.gz',
    threads: 1
    log:
        'logs/reads/raw/{sample_id}_concat_paired.log'
    shell:
        'cat {params.space_sepped_R1_fastqs} > reads/raw/raw_{wildcards.sample_id}_R1_paired.fastq.gz; '
        'cat {params.space_sepped_R2_fastqs} > reads/raw/raw_{wildcards.sample_id}_R2_paired.fastq.gz; '


rule trimmomatic_pe:
    input:
        r1='reads/raw/raw_{sample_id}_R1_paired.fastq.gz',
        r2='reads/raw/raw_{sample_id}_R2_paired.fastq.gz',
    output:
        r1='reads/trimmed/trimmed_{sample_id}_R1_paired.fastq.gz',
        r2='reads/trimmed/trimmed_{sample_id}_R2_paired.fastq.gz',
        # reads where trimming entirely removed the mate
        r1_unpaired='reads/trimmed/trimmed_{sample_id}_R1_unpaired.fastq.gz',
        r2_unpaired='reads/trimmed/trimmed_{sample_id}_R2_unpaired.fastq.gz',
    threads: 4
    resources:
        walltime = 600
    log:
        'logs/trimmomatic/{sample_id}.log'
    params:
        # list of trimmers (see manual)
        trimmer=[trimmer.replace('illumina_adapters',
                                 config['illumina_adapters'])
                 for trimmer in config['params']['trimmomatic']['trimmer']],
        # optional parameters
        extra=config['params']['trimmomatic']['extra'],
    wrapper:
        f'{wrapper_version}/bio/trimmomatic/pe'

rule fastqc:
    input:
        '{folder}/{read}'
    output:
        html='read_QC/{folder}/{read}.html',
        zip='read_QC/{folder}/{read}.zip'
    params: ''
    threads: 1
    resources:
        walltime = 600
    log:
        'logs/read_qc/{folder}/{read}.log'
    wrapper:
        f'{wrapper_version}/bio/fastqc'
