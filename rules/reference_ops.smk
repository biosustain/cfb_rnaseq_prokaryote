rule generate_gff3_reference:
    input:
        reference_gff
    output:
        'references/{reference_id}_processed.gff'
    conda: '../envs/gffread.yaml'
    threads: 1
    log:
        'logs/gffread/{reference_id}_processed.log'
    params:
        '-F --keep-exon-attrs'
    shell:
        'gffread {input} -o {output} {params} > {log} 2>&1'


rule generate_masked_reference:
    input:
        'references/{reference_id}_processed.gff'
    output:
        'references/{reference_id}_filtered.gff'
    conda: '../envs/gffread.yaml'
    threads: 1
    log:
        'logs/gffread/{reference_id}_filtered.log'
    params:
        '-CF --no-pseudo --keep-exon-attrs'
    shell:
        'gffread {input} -o {output} {params} > {log} 2>&1'


# Some pseudogenes are not marked. They don't seem to have a protein_id, which
# is the attribute we want to use. Here, I eliminate those lines.
rule generate_filtered_reference:
    input:
         'references/{reference_id}_filtered.gff'
    output:
        'references/{reference_id}_masked.gff'
    threads: 1
    log:
        'logs/gffread/{reference_id}_masked.log'
    run:
        with open(f'{input}', 'r') as infile:
            with open(f'{output}', 'w') as outfile:
                for line in infile.readlines():
                    if ('\tCDS' in line) and (';protein_id=' not in line):
                        continue
                    outfile.writelines(line)


rule generate_mask_reference:
    input:
        in1='references/{reference_id}_processed.gff',
        in2='references/{reference_id}_masked.gff'
    output:
        'references/{reference_id}_mask.gff'
    conda:
        '../envs/bedtools.yaml'
    threads: 1
    log:
        'logs/bedtools/subtract/{reference_id}.log'
    params:
        ''
    shell:
         """
bedtools subtract {params} \
-a {input.in1} -b {input.in2} \
1> {output} \
2> {log}
        """


rule generate_bed:
    input:
        'references/{reference}.gff'
    output:
        'references/{reference}.bed'
    conda:
        '../envs/bedops.yaml'
    threads: 1
    log:
        'logs/bedops/gff2bed/{reference}.log'
    params:
        ''
    shell:
        """
(gff2bed {params} \
< {input} > {output}) \
> {log}
        """
