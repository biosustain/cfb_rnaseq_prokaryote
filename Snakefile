# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import yaml
import re
from glob import glob


# it's better to assign configfile from the command line with "--configfile"
# configfile: "config.yaml"

# report: "report/workflow.rst"  # not working of this yet

# get snakemake wrapper version to be used in this pipeline
wrapper_version = config['wrapper_version']

# get sample and reference details through config input
sample_id = config['sample_id']
raw_read_path = config['raw_read_path']
reference_id = config['reference_id']
reference_table = config['reference_table']

# find reference fasta and gff
with open(reference_table, 'r') as file:
    reference_table = yaml.load(file, yaml.FullLoader)
    reference_fasta = reference_table[reference_id]['fasta']
    reference_gff = reference_table[reference_id]['gff']

# find sample read files under raw_read_path
sample_to_files = {sample_id: []}
# re pattern for filenames in samplename_S#_L00#_R#_00#.fastq.gz format
# (recent illumina filename format)
p1 = re.compile('^(.*)(_S[0-9]+)(_L00[1-4])(_R[12])(_[0-9]{3}[\S]*\.fastq\.gz)$')
# older file name patterns below are not supported any longer:
# re pattern for filenames in samplename_S#_R#_00#.fastq.gz format.
# (probably pre_nextseq)
# p2 = re.compile('^(.*)(_S[0-9]+)(_R[12])(_[0-9]{3}[\S]*\.fastq\.gz)$')
# re pattern for filenames in sampleName.f.slx.gz format.
# (probably very early solexa reads) ([\S]* is to accommodate for trimmed reads)
# p3 = re.compile('^(.*)(\.f\.slx[\S]*\.gz)$')
sampleNumbers = set()
for x in glob('{}/**/*.gz'.format(raw_read_path), recursive=True):
    filename = x.split('/')[-1]
    pm = p1.match(filename)

    if pm:
        sampleName = pm.group(1)
        if sampleName == sample_id:
            sampleNumbers.add(pm.group(2))
            sample_to_files[sampleName].append(x)

# make sure there are read files for the sample:
assert (len(sample_to_files[sample_id]) > 0), \
    f'No read files matching {sample_id} found.'

# make sure read files have the same sample illumina sample number (..._S#_...)
assert (len(sampleNumbers) == 1), \
    f'Illumina sample numbers are not unique for {sample_id}: {sampleNumbers}'

# keep things sorted
sample_to_files[sample_id] = sorted(sample_to_files[sample_id])


rule all:
    input:
        # hisat2 alignments
        f'hisat2/align/{sample_id}.accepted_hits.bam',
        # counts ready for db
        f'hisat2/align/{sample_id}.counts_for_db.csv',


# The first rule should define the default target files
# Subsequent target rules can be specified below. They should start with all_*.


rule all_read_ops:
    input:
        # Concat reads from different lanes (somewhat encumbered names - helps later)
        expand('reads/raw/raw_{sample_id}_{end}_paired.fastq.gz',
               sample_id=[sample_id],
               end=['R1', 'R2']),
        # Get trimmed reads
        expand('reads/trimmed/trimmed_{sample_id}_{end}_{pairing}.fastq.gz',
               sample_id=[sample_id],
               end=['R1', 'R2'],
               pairing=['paired', 'unpaired']),
        # QC for raw concat'ed reads
        expand('read_QC/{folder}/{read}.{type}',
               folder=['reads/raw'],
               read=['_'.join(['raw', sample_id, end, 'paired'])+'.fastq.gz'
                     for end in ['R1', 'R2']],
               type=['html', 'zip']),
        # QC for trimmed reads
        expand('read_QC/{folder}/{read}.{type}',
               folder=['reads/trimmed'],
               read=['_'.join(['trimmed', sample_id, end, pairing])+'.fastq.gz'
                     for end in ['R1', 'R2']
                     for pairing in ['paired', 'unpaired']],
               type=['html', 'zip']),


rule all_reference_ops:
    input:
        # gffread processed reference
        f'references/{reference_id}_processed.gff',
        # filtered reference
        f'references/{reference_id}_filtered.gff',
        # masked reference
        f'references/{reference_id}_masked.gff',
        # mask gff reference
        f'references/{reference_id}_mask.gff',
        # bed files
        expand('references/{reference_id}_{type}.bed',
               reference_id=[reference_id],
               type=['processed', 'masked', 'mask']),


rule all_alignment_ops:
    input:
        # index files will be hidden under this folder (eg. ".1.ht2")
        f'hisat2/index',
        # hisat2 alignments
        f'hisat2/align/{sample_id}.accepted_hits.bam',
        # hisat2 alignments summary
        f'hisat2/align/{sample_id}.accepted_hits.align_summary.txt',
        # index alignments
        f'hisat2/align/{sample_id}.accepted_hits.bam.bai',
        # alignments to mask
        f'hisat2/align/{sample_id}.mask.bam',
        # generate flagstat for mask and reference
        f'hisat2/align/{sample_id}.accepted_hits.flagstat',
        f'hisat2/align/{sample_id}.mask.flagstat',


rule all_count_ops:
    input:
        # featureCounts output
        f'hisat2/align/{sample_id}.featureCounts.out',
        f'hisat2/align/{sample_id}.featureCounts.out.summary',
        # counts ready for db
        f'hisat2/align/{sample_id}.counts_for_db.csv',


include: 'rules/read_ops.smk'
include: 'rules/reference_ops.smk'
include: 'rules/alignment_ops.smk'
include: 'rules/count_ops.smk'


