rule featureCounts:
    input:
        reference = f'references/{reference_id}_masked.gff',
        alignment = 'hisat2/align/{sample_id}.accepted_hits.bam'
    output:
        output = 'hisat2/align/{sample_id}.featureCounts.out',
        summary = 'hisat2/align/{sample_id}.featureCounts.out.summary'
    conda:
        '../envs/featureCounts.yaml'
    threads: 4
    params:
        config['params']['featureCounts']
    log:
        'logs/feature_counts/{sample_id}.featureCounts.log'
    shell:
        """
(featureCounts \
{params} \
-T {threads} \
-a {input.reference} \
-o {output.output} \
{input.alignment}) \
> {log} 2>&1
        """


rule counts_for_db:
    input:
        'hisat2/align/{sample_id}.featureCounts.out'
    output:
        'hisat2/align/{sample_id}.counts_for_db.csv'
    threads: 1
    log:
        'logs/feature_counts/{sample_id}.counts_for_db.csv'
    run:
        import pandas as pd
        df = pd.read_csv(input[0],
                         sep='\t',
                         skiprows=1,
                         index_col=0,
                         )
        df.index.rename('Feature', inplace=True)
        df.rename(columns=lambda col: 'Count' if col.startswith('hisat2')
                                      else col,
                  inplace=True)
        df['Sample'] = sample_id
        df[['Sample', 'Count']].to_csv(output[0])
