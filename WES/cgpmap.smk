WORK="/work/sduvarcall/"
INPUT = WORK+"G56-2016-Glioblastom/WES/fastq/"
OUTPUT = WORK+"G56-2016-Glioblastom/WES/bam/"
YAML = WORK + "G56-2016-Glioblastom/WES/ReadGroupsYAML/"

SAMPLES, = glob_wildcards(INPUT+"{sample}_1.fastq.gz")
# SAMPLES = "G56-sampleC1_truseq-nano-genome_HGL2LDSXX_S3"
# SAMPLES = "G56-sampleE4_truseq-nano-genome_HGL2LDSXX_S18"
# SAMPLES = "G56-sampleA8_truseq-nano-genome_HGL2LDSXX_S20"

# Resources
ref = WORK+"G37-2016-BRCAX-Genomes/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv"
ref_fa = WORK+"G37-2016-BRCAX-Genomes/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa"

# ref = WORK+"G37-2016-BRCAX-Genomes/reference_files/core_ref_GRCh37d5"
resource_path = WORK+"Resources/hg38/"
dbsnp = resource_path + "Homo_sapiens_assembly38.dbsnp138.vcf"
mills_1000g = resource_path + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
phase1_1000g = resource_path + "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
# bed = resource_path + "MedExome_hg38_capture_targets.bed"
# interval_list = resource_path + "MedExome_hg38_capture_targets.interval_list"

mem = "12"

rule all:
    input:
        # expand(OUTPUT+"{sample}.bam", sample=SAMPLES)
        expand(OUTPUT+"quality_control/{sample}_post_recalibration.grp", sample=SAMPLES)

## IMPORTANT NOTE: Make sure fastq files ends on _1.fastq.gz and _2.fastq.gz
## - otherwise faulty bam files are produced
rule cgpmap:
    input:
        f1=INPUT+"{sample}_1.fastq.gz",
        f2=INPUT+"{sample}_2.fastq.gz",
        rg=YAML+"ReadGroups_{sample}.yaml"
    output:
        bam=OUTPUT+"{sample}.bam",
        bai=OUTPUT+"{sample}.bam.bai"
    threads: 30
    shell:
        """
        python2 /work/sduvarcall/udocker/udocker run \
        --user=laran \
        --volume={WORK}:{WORK} \
        cgpmap \
        ds-cgpmap.pl \
        -reference {ref} \
        -bwa_idx {ref} \
        -sample {wildcards.sample} \
        -groupinfo {input.rg} \
        -threads {threads} \
        -outdir {OUTPUT} \
        {input.f1} {input.f2}
        """


## BaseRecalibrator

rule BaseRecalibrator_markdup:
    input:
        bam = OUTPUT+"{sample}.bam",
        bai = OUTPUT+"{sample}.bam.bai",
    output:
        OUTPUT+"quality_control/{sample}_pre_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx{mem}G BaseRecalibrator \
        -R {ref_fa} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output}
        """

rule ApplyBQSR_markdup:
    input:
        bam = OUTPUT+"{sample}.bam",
        bai = OUTPUT+"{sample}.bam.bai",
        recal = OUTPUT+"quality_control/{sample}_pre_recalibration.grp"
    output:
        bam = OUTPUT+"{sample}.recalibrated.bam",
        bai = OUTPUT+"{sample}.recalibrated.bam.bai"
    params:
        bai = OUTPUT+"{sample}.recalibrated.bai"
    shell:
        """
        gatk --java-options -Xmx{mem}G ApplyBQSR \
        -R {ref_fa} \
        -I {input.bam} \
        --bqsr {input.recal} \
        -O {output.bam}

        cp {params.bai} {output.bai}
        """

rule post_recalibration_table_markdup:
    input:
        bam = OUTPUT+"{sample}.recalibrated.bam",
        bai = OUTPUT+"{sample}.recalibrated.bam.bai"
    output:
        post = OUTPUT+"quality_control/{sample}_post_recalibration.grp"
    shell:
        """
        gatk --java-options -Xmx{mem}G BaseRecalibrator \
        -R {ref_fa} \
        -I {input.bam} \
        --known-sites {dbsnp} \
        --known-sites {mills_1000g} \
        --known-sites {phase1_1000g} \
        -O {output}
        """
