snakemake -s scripts/varscan2-joint-calling-wes-somatic-hg38.smk -j 999 --cluster "sbatch -c 1 --mem 6g"