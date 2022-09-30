# conda install mamba
# mamba install -c bioconda snakemake

# snakemake -s /work/sduvarcall/G56-2016-Glioblastom/WES/scripts/sequenza.smk -j 30 --use-conda --conda-frontend mamba --cluster "sbatch -c 1 --mem 6g"
snakemake -s /work/sduvarcall/G56-2016-Glioblastom/WES/scripts/sequenza.smk -j 30 --cluster "sbatch -c 1 --mem 6g"