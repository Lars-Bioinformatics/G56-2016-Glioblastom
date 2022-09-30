
WGS_BAMS, = glob_wildcards("bam/{sample}/{sample}.bam")
WES_BAMS, = glob_wildcards("WES/bam/{sample}.bam")

print(WGS_BAMS)
print(WES_BAMS)

rule all:
    input:
        expand("{sample}_merged_wes_wgs.bam", 

rule create_rg_tag_file:
    input:
        wgs_bam="bam/{sampleid}_{protocol_wgs}_{flowcell_wgs}_{sid}/{sampleid}_{protocol_wgs}_{flowcell_wgs}_{sid}.bam",
        wes_bam="WES/bam/{sampleid}_{protocol_wes}_{flowcell_wes}.bam"
    output:

    shell:
        """
        \"@RG\\tID:{params.rgid}\\tLB:{params.rglb}\\tSM:{params.rgsm}\\tPL:{params.rgpl}\\tPU:{params.rgpu}\"
        perl -e 'print "@RG\\tID:ga\\tSM:hs\\tLB:ga\\tPL:Illumina\\n@RG\\tID:454\\tSM:hs\\tLB:454\\tPL:454\\n"' > {wildcards.sample}_rg.txt
        """

rule merge_bam:
    input:
        wgs_bam="bam/{sampleid}_{protocol_wgs}_{flowcell_wgs}_{sid}/{sampleid}_{protocol_wgs}_{flowcell_wgs}_{sid}.bam",
        wes_bam="WES/bam/{sampleid}_{protocol_wes}_{flowcell_wes}.bam"
    output:

    shell:
        """
        samtools merge -@2 -rh rg.txt -o {output} {input}
        """