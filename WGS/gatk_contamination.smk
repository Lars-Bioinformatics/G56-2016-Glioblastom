# Rules in Snakefile based on https://github.com/cancerit/dockstore-cgpwgs/blob/develop/scripts/analysisWGS.sh

WORK="/work/sduvarcall/"
INPUT_BAM = WORK+"G56-2016-Glioblastom/bam"
OUTPUT = WORK+"G56-2016-Glioblastom/gatk_contamination"

# SAMPLES, = glob_wildcards(INPUT+"{sample}_R1_001.fastq.gz")

# NORMAL = "G45-ECV2-31-blod_truseq-nano-genome_HGL2LDSXX_S27"
# TUMOR = "G45-ECV2-31-biopsi-F1_truseq-nano-genome_HGL2LDSXX_S22"

# NORMAL = "G56-blod-pt1_truseq-nano-genome_HGL2LDSXX_S5"
# TUMOR = "G56-sampleA1_truseq-nano-genome_HGL2LDSXX_S1"

configfile: WORK+"G56-2016-Glioblastom/samples.yaml"

ref_build = "GRCh38"

# Resources - paths inside docker
if ref_build == "GRCh37":
    # GRCh37 paths
    ref = WORK+"G37-2016-BRCAX-Genomes/reference_files/core_ref_GRCh37d5"
    res = WORK+"G37-2016-BRCAX-Genomes/reference_files"
    common_variants = "/work/sduvarcall/knownSNPs/gnomad/somatic-b37_small_exac_common_3.vcf"
else:
    # GRCh38 paths
    ref = WORK+"G37-2016-BRCAX-Genomes/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
    res = WORK+"G37-2016-BRCAX-Genomes/reference_files_GRCh38"
    common_variants = "/work/sduvarcall/Resources/hg38/knownSNPs/small_exac_common_3.hg38.vcf.gz"

mem = "-Xmx12g"

onstart:
    shell("mkdir -p " + OUTPUT)


rule all:
    input:
        [expand(OUTPUT+"/{tumor}_vs_{normal}_contamination.table", tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in config]


###############################################################################
#### Calculate GATK contamination
###############################################################################
rule GetPileupSummaries:
    input:
        bam=INPUT_BAM+"/{sample}/{sample}.bam"
    output:
        pileup=OUTPUT+"/{sample}_pileup.table"
    shell:
        """
        gatk --java-options {mem} GetPileupSummaries \
        -I {input.bam} \
        -V {common_variants} \
        -L {common_variants} \
        -O {output}
        """

rule CalculateContamination:
    input:
        normal=OUTPUT+"/{normal}_pileup.table",
        tumor=OUTPUT+"/{tumor}_pileup.table"
    output:
        contamination=OUTPUT+"/{tumor}_vs_{normal}_contamination.table"
    shell:
        """
        gatk --java-options {mem} CalculateContamination \
        -I {input.tumor} \
        -matched {input.normal} \
        -O {output}
        """
