FROM snakemake/snakemake:v5.7.0
MAINTAINER Emre Ozdemir <emoz@biosustain.dtu.dk>
WORKDIR /
RUN git clone https://github.com/meono/cfb_rnaseq_prokaryote.git
WORKDIR /cfb_rnaseq_prokaryote
RUN snakemake all \
              --snakefile /cfb_rnaseq_prokaryote/Snakefile \
              --configfile /cfb_rnaseq_prokaryote/config.yaml \
              --directory /cfb_rnaseq_prokaryote/test \
              --config sample_id=ecoli_01 \
                       raw_read_path=/cfb_rnaseq_prokaryote/test/data \
                       reference_id=GCF_000005845.2 \
                       reference_table=reference_table.yaml \
              --use-conda \
              --conda-prefix /conda_envs \
              --create-envs-only