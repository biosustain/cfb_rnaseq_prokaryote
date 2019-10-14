# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import yaml
# from glob import glob
# from os.path import basename, splitext
# from os import symlink, makedirs, chdir
# import pandas as pd

configfile: "config.yaml"
# report: "report/workflow.rst"

with open('sample.yaml', 'r') as file:
    sample_details = yaml.load(file, yaml.FullLoader)
print (sample_details)

# sample = dict()
# # re pattern for filenames in samplename_S#_L00#_R#_00#.fastq.gz format (recent illumina filename format)
# p1 = re.compile('^(.*)(_S[0-9]+)(_L00[1-4])(_R[12])(_[0-9]{3}[\S]*\.fastq\.gz)$')
# # older file name patterns - not supported any longer:
# # re pattern for filenames in samplename_S#_R#_00#.fastq.gz format. (probably pre_nextseq)
# # p2 = re.compile('^(.*)(_S[0-9]+)(_R[12])(_[0-9]{3}[\S]*\.fastq\.gz)$')
# # re pattern for filenames in samplename.f.slx.gz format. (probably very early solexa reads)
# # p3 = re.compile('^(.*)(\.f\.slx[\S]*\.gz)$')
# # [\S]* is to accommodate for trimmed reads
# for x in glob('{}/**/*.gz'.format(config['raw_read_path']), recursive=True):
#     filename = x.split('/')[-1]
#     pm = p1.match(filename)
#     if pm:
#         group = pm.group(1)
#     else:
#         raise ValueError('Filename does not match recognized format for grouping:\n{}'.format(x))
#         # pm = p2.match(filename)
#         # if pm:
#         #     group = pm.group(1)
#         # else:
#         #     pm = p3.match(filename)
#         #     if pm:
#         #         group = pm.group(1)
#         #     else:
#         #         raise ValueError('Filename does not match recognized format for grouping:\n{}'.format(x))
#     if group in group_to_files:
#         group_to_files[group].append(x)
#     else:
#         group_to_files[group] = [x]
# # keep things sorted
# for group, reads in group_to_files.items():
#     group_to_files[group] = sorted(reads)

rule all:
    input:
        'test.out'
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

rule test:
    input:
        'sample_template.yaml'
    output:
        'test.out'
    shell:
        'touch test.out'


# include: 'rules/read_ops.smk'
# include: 'rules/reference_ops.smk'
# include: 'rules/alignment_ops.smk'

