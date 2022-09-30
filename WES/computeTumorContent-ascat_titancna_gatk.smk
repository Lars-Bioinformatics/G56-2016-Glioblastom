# Rules in Snakefile based on https://github.com/cancerit/dockstore-cgpwgs/blob/develop/scripts/analysisWGS.sh

WORK="/work/sduvarcall/"
INPUT_BAM = WORK+"G56-2016-Glioblastom/WES/bam"
OUTPUT_ASCAT = WORK+"G56-2016-Glioblastom/WES/ascat"
OUTPUT_GATK = WORK+"G56-2016-Glioblastom/WES/gatk_contamination"
OUTPUT_TitanCNA = WORK+"G56-2016-Glioblastom/WES/TitanCNA"

configfile: "/work/sduvarcall/G56-2016-Glioblastom/WES/samples_b37.yaml"
titan_samples = "samples_G56-Glioblastom-Exomes.yaml"

REF_FILES = "/work/sduvarcall/G37-2016-BRCAX-Genomes/"
TITAN_DIR="/work/sduvarcall/G37-2016-BRCAX-Genomes/TitanCNA/scripts/snakemake/config"
ref_build = "GRCh37"

# Resources - paths inside docker
if ref_build == "GRCh37":
    # GRCh37 paths
    ref = REF_FILES+"reference_files/core_ref_GRCh37d5"
    res = REF_FILES+"reference_files"
else:
    # GRCh38 paths
    ref = REF_FILES+"/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
    res = REF_FILES+"/reference_files_GRCh38"

# Global variables
if ref_build == "GRCh37":
    SPECIES = "Human"   # Specify according to bam file e.g. Human or human
    ASSEMBLY = "NCBI37" # Specify according to bam file e.g. GRCh38 or NCBI37
    common_variants = "/work/sduvarcall/knownSNPs/gnomad/somatic-b37_small_exac_common_3.vcf"
    titan_config = "config_b37.yaml"
else:
    SPECIES = "human"   # Specify according to bam file e.g. Human or human
    ASSEMBLY = "GRCh38" # Specify according to bam file e.g. GRCh38 or NCBI37
    common_variants = "/work/sduvarcall/Resources/hg38/knownSNPs/small_exac_common_3.hg38.vcf.gz"
    titan_config = "config_hg38.yaml"

PROTOCOL = "WXS"

ASCAT_ADD_ARGS='' # Let ASCAT compute purity and ploidy
# ASCAT_ADD_ARGS='-pu ASCAT_PURITY -pi ASCAT_PLOIDY'
# ASCAT_ADD_ARGS='-pu 0.8 -pi 3'

# Get sample names from config file
SAMPLES = [sample for sample in config]
# SAMPLES = ["G90-VHL521-10-CNS1"]
print(SAMPLES)

mem = "-Xmx12G"

onstart:
    shell("mkdir -p " + OUTPUT_ASCAT)
    shell("mkdir -p " + OUTPUT_GATK)
    shell("mkdir -p " + OUTPUT_TitanCNA)


rule all:
    input:
        [expand(OUTPUT_ASCAT+"/{tumor}/timings/"+PROTOCOL+"_{tumor}_vs_{normal}.time.ascat", 
            tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in SAMPLES],
        OUTPUT_GATK+"/GATK4_merged_contamination.table",
        OUTPUT_TitanCNA+"/results/titan/hmm/optimalClusterSolution.txt"



###############################################################################
#### Ascat
###############################################################################
rule ascat:
    input:
        normal=INPUT_BAM+"/{normal}.bam",
        tumor=INPUT_BAM+"/{tumor}.bam"
    output:
        res3_1=OUTPUT_ASCAT+"/{tumor}/timings/"+PROTOCOL+"_{tumor}_vs_{normal}.time.ascat"
    threads: 24
    shell:
        """
        echo -e "Running ASCAT"; \
        python2 ~/udocker-1.1.4/udocker run \
        --user=pi \
        --volume={WORK}:{WORK} \
        cgpwgs \
        bash -c '/usr/bin/time -v ascat.pl \
        -o {OUTPUT_ASCAT}/{wildcards.tumor}/{PROTOCOL}_{wildcards.tumor}/ascat \
        -t {input.tumor} \
        -n {input.normal} \
        -sg {res}/ascat/SnpGcCorrections.tsv \
        -r {ref}/genome.fa \
        -q 20 \
        -g L \
        -l {res}/qcGenotype/gender.tsv \
        -rs {SPECIES} \
        -ra {ASSEMBLY} \
        -pr {PROTOCOL} \
        -pl ILLUMINA \
        -c {threads} \
        -force \
        {ASCAT_ADD_ARGS}' \
        2>&1 | tee {output.res3_1}
        """


###############################################################################
#### Calculate GATK contamination
###############################################################################
rule GetPileupSummaries:
    input:
        bam=INPUT_BAM+"/{sample}.bam"
    output:
        pileup=OUTPUT_GATK+"/{sample}_pileup.table"
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
        normal=OUTPUT_GATK+"/{normal}_pileup.table",
        tumor=OUTPUT_GATK+"/{tumor}_pileup.table"
    output:
        cont_table=OUTPUT_GATK+"/{tumor}_vs_{normal}_contamination.table"
    shell:
        """
        gatk --java-options {mem} CalculateContamination \
        -I {input.tumor} \
        -matched {input.normal} \
        -O {output}
        """

rule JoinContaminationTables:
    input:
        cont_tables=[expand(OUTPUT_GATK+"/{tumor}_vs_{normal}_contamination.table", 
                    normal=config[s]["normal"], tumor=config[s]["tumor"]) for s in config]
    output:
        joined_table=OUTPUT_GATK+"/GATK4_merged_contamination.table"
    shell:
        """
        echo -e 'sample\tcontamination\terror' > {output};
        awk -F'\t' 'FNR == 2' {input} >> {output}
        """


###############################################################################
#### Calculate GATK contamination
###############################################################################
rule TitanCNA:
    input:
        TITAN_DIR+"/"+titan_samples,
        TITAN_DIR+"/"+titan_config
    output:
        OUTPUT_TitanCNA+"/results/titan/hmm/optimalClusterSolution.txt"
    shell:
        """
        \cp {TITAN_DIR}/{titan_samples} {TITAN_DIR}/samples.yaml
        \cp {TITAN_DIR}/{titan_config} {TITAN_DIR}/config.yaml
        cd {OUTPUT_TitanCNA}
        bash /work/sduvarcall/G37-2016-BRCAX-Genomes/TitanCNA/scripts/snakemake/runTitanCNA.sh
        """