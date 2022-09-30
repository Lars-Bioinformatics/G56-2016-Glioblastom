WORK="/work/sduvarcall/"
INPUT = WORK+"G56-2016-Glioblastom/fastq/"
OUTPUT = WORK+"G56-2016-Glioblastom/bam_noAlt/"
YAML = WORK+"G56-2016-Glioblastom/ReadGroupsYAML/"

SAMPLES, = glob_wildcards(INPUT+"{sample}_1.fastq.gz")
# SAMPLES = ["G56-sampleC3_truseq-nano-genome_HGL2LDSXX_S12", "G56-blod-pt3_truseq-nano-genome_HGL2LDSXX_S14"]
SAMPLES = "G56-blod-pt3_truseq-nano-genome_HGL2LDSXX_S14"
# SAMPLES = "G56-sampleA8_truseq-nano-genome_HGL2LDSXX_S20"

# Resources - paths inside docker
# ref = WORK+"G37-2016-BRCAX-Genomes/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
ref = WORK+"Resources/hg38/GRCh38_no_alt_plus_hd38d1_reference"
# ref = WORK+"G37-2016-BRCAX-Genomes/reference_files/core_ref_GRCh37d5"

rule all:
    input:
        expand(OUTPUT+"{sample}/{sample}.bam", sample=SAMPLES)

## IMPORTANT NOTE: Make sure fastq files ends on _1.fastq.gz and _2.fastq.gz
## - otherwise faulty bam files are produced
rule cgpmap:
    input:
        f1=INPUT+"{sample}_1.fastq.gz",
        f2=INPUT+"{sample}_2.fastq.gz",
        rg=YAML+"ReadGroups_{sample}.yaml"
    output:
        OUTPUT+"{sample}/{sample}.bam"
    threads: 24
    shell:
        """
        python2 ~/udocker-1.1.4/udocker run \
        --user=laran \
        --volume={WORK}:{WORK} \
        cgpmap \
        ds-cgpmap.pl \
        -reference {ref} \
        -bwa_idx {ref} \
        -sample {wildcards.sample} \
        -groupinfo {input.rg} \
        -threads {threads} \
        -outdir {OUTPUT}{wildcards.sample} \
        {input.f1} {input.f2}
        """