rule pl_convert_gb_to_gff:
    input:
        lambda wildcards: [ref for ref in config['references']['gb_refs'] 
                                 if splitext(basename(ref))[0] == wildcards.reference][0]
    output:
        'references/gff3_from_gb/{reference}.gff3'
    conda:
        '../envs/pl_gb2gff3.yaml'
    threads: 1
    log:
        'logs/pl_convert_gb_to_gff/{reference}.log'
    shell:
    # have to include the BIO and LWP modules from 5.22.0 
    # since conda perl-bioperl kind of broken.
        'perl -I ${{CONDA_PREFIX}}/lib/perl5/site_perl/5.22.0 '
        '${{CONDA_PREFIX}}/bin/bp_genbank2gff3.pl {input} '
        '-o stdout 1> {output} '
        '2> {log} '

rule fasta_from_genbank:
    input:
        lambda wildcards: [ref for ref in config['references']['gb_refs'] 
                                 if splitext(basename(ref))[0] == wildcards.reference][0]
    output:
        'references/fasta_from_gb/{reference}.fa'
    threads: 1
    log:
        'logs/readseq/{reference}.log'
    params:
        '-C -f fa'
    wrapper:
        f'{wrapper_source}/bio/readseq'
    
rule concat_reference_gff3s:
    input:
        ([ref for ref in config['references']['gff_refs']] 
              if config['references']['gff_refs'] != ['None'] else []) + 
        (['references/gff3_from_gb/'+splitext(basename(ref))[0]+'.gff3' 
              for ref in  config['references']['gb_refs']] 
              if config['references']['gb_refs'] != ['None'] else []) 
    output:
        'references/combined.gff3'
    shell:
        'cat {input} > references/combined.gff3'

rule concat_reference_fastas:
    input:
        ([ref for ref in config['references']['fa_refs']] 
              if config['references']['fa_refs'] != ['None'] else []) + 
        (['references/fasta_from_gb/'+splitext(basename(ref))[0]+'.fa' 
              for ref in  config['references']['gb_refs']] 
              if config['references']['gb_refs'] != ['None'] else []) 
    output: # generates a full sequence, but the name has to match the gff reference used.
        'references/combined_only_CDS_no_pseudo.fa'
    threads: 1
    log:
        'logs/concat_reference_fastas.log'
    shell:
        'cat {input} > references/combined_only_CDS_no_pseudo.fa'

rule generate_gff3_reference:
    input:
        'references/combined.gff3'
    output:
        'references/combined_gffread.gff3'
    threads: 1
    log:
        'logs/gffread/combined_gffread.log'
    params:
        '-F'
    wrapper:
        f'{wrapper_source}/bio/gffread'

rule generate_masked_reference:
    input:
        'references/combined_gffread.gff3'
    output:
        'references/combined_only_CDS_no_pseudo.gff3'
    threads: 1
    log:
        'logs/gffread/combined_only_CDS_no_pseudo.log'
    params:
        '-CF --no-pseudo'
    wrapper:
        f'{wrapper_source}/bio/gffread'    

rule generate_mask_reference:
    input:
        in1='references/combined_gffread.gff3',
        in2='references/combined_only_CDS_no_pseudo.gff3'
    output:
        'references/combined_mask.gff3'
    threads: 1
    log:
        'logs/bedtools/subtract/combined_mask.log'
    params:
        ""
    wrapper:
        f'{wrapper_source}/bio/bedtools/subtract'

rule generate_mask_bed:
    input:
        'references/{reference}.gff3'
    output:
        'references/{reference}.bed'
    threads: 1
    log:
        'logs/bedops/gff2bed/{reference}.log'
    params:
        ""
    wrapper:
        f'{wrapper_source}/bio/bedops/gff2bed'

rule generate_mask_fasta:
    input:
        in1='references/combined_only_CDS_no_pseudo.fa',
        in2='references/combined_mask.gff3'
    output:
        'references/combined_mask.fa'
    threads: 1
    log:
        'logs/bedtools/getfasta/combined_mask.log'
    params:
        ""
    wrapper:
        f'{wrapper_source}/bio/bedtools/getfasta'
