rule hisat2_index:
    input:
        reference_fasta
    output:
        directory('hisat2/index')
    conda:
        '../envs/hisat2.yaml'
    log:
        'logs/hisat2/index.log'
    params:
        prefix = reference_id,
        extra = ''
    threads: 4
    shell:
        """
(mkdir -p {output} && \
hisat2-build {params.extra} \
-p {threads} \
-f {input} {output}/{params.prefix}) \
> {log} 2>&1
        """

rule hisat2_align: # not satisfied with the input structure of this wrapper.
    input:
        idx = 'hisat2/index',
        r1 = 'reads/trimmed/trimmed_{sample_id}_R1_paired.fastq.gz',
        r1u = 'reads/trimmed/trimmed_{sample_id}_R1_unpaired.fastq.gz',
        r2 = 'reads/trimmed/trimmed_{sample_id}_R2_paired.fastq.gz',
        r2u = 'reads/trimmed/trimmed_{sample_id}_R2_unpaired.fastq.gz'
    output:
        bam = 'hisat2/align/{sample_id}.accepted_hits.bam',
        summary = 'hisat2/align/{sample_id}.accepted_hits.align_summary.txt',
    conda:
        '../envs/hisat2.yaml'
    log:
        'logs/hisat2/align_{sample_id}.log'
    params:
         reference_id = reference_id,
         extra = config['params']['hisat2']['align']['extra']
    threads: 4
    shell:
        """
(hisat2 {params.extra} \
--threads {threads} \
-x {input.idx}/{params.reference_id} \
-1 {input.r1} \
-2 {input.r2} \
-U {input.r1u},{input.r2u} \
--summary-file {output.summary} \
| samtools sort \
| samtools view -hb \
1> {output.bam} ) \
> {log} 2>&1
        """


rule index_bam:
    input:
        'hisat2/align/{sample_id}.accepted_hits.bam'
    output:
        'hisat2/align/{sample_id}.accepted_hits.bam.bai'
    log:
        'logs/samtools/{sample_id}.accepted_hits.index.log'
    threads: 1
    params:
        ''
    wrapper:
        f'{wrapper_version}/bio/samtools/index'


rule extract_mask_alignment:
    input:
        'hisat2/align/{sample_id}.accepted_hits.bam',
        f'references/{reference_id}_mask.bed'
    output:
        'hisat2/align/{sample_id}.mask.bam'
    log:
        'logs/samtools/{sample_id}_alignments_to_mask.log'
    threads: 1
    params:
        f'-hbL references/{reference_id}_mask.bed'
    wrapper:
        f'{wrapper_version}/bio/samtools/view'


rule flagstat:
    input:
        'hisat2/align/{alignment}.bam',
    output:
        'hisat2/align/{alignment}.flagstat'
    log:
        'logs/samtools/{alignment}.flagstat.log'
    threads: 1
    params:
        ''
    wrapper:
        f'{wrapper_version}/bio/samtools/flagstat'
