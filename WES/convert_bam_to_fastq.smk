INPUT="bam_GRCh37/"
OUTPUT="fastq/"

SAMPLES, = glob_wildcards(INPUT+"{sample}.bam")

onstart:
    shell("mkdir -p fastq")

rule all:
    input:
        expand(OUTPUT+"{sample}_1.fastq.gz", sample=SAMPLES)


rule samtools_bam_to_fastq:
    input:
        bam=INPUT+"{sample}.bam"
    output:
        fq1=OUTPUT+"{sample}_1.fastq.gz",
        fq2=OUTPUT+"{sample}_2.fastq.gz"
    params:
        bam=INPUT+"{sample}_sorted.bam"
    threads: 2
    shell:
        """
        samtools sort -n {input.bam} -o {params.bam}
        samtools fastq -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null -@ {threads} {params.bam}
        """

# rule bamToFastq:
#     input:
#         bam=INPUT+"{sample}.bam"
#     output:
#         fq1=OUTPUT+"{sample}_1.fastq.gz",
#         fq2=OUTPUT+"{sample}_2.fastq.gz"
#     params:
#         fq1=OUTPUT+"{sample}_1.fastq",
#         fq2=OUTPUT+"{sample}_2.fastq"
#     threads: 1
#     shell:
#         """
#         bamToFastq -i {input.bam} -fq {params.fq1} -fq2 {params.fq2}

#         bgzip {params.fq1}
#         bgzip {params.fq2}
#         """
