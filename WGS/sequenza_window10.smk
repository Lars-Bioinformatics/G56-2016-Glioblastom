from snakemake.utils import R

base_window = "10" # default is 50 - and used here: https://github.com/mbourgey/EBI_cancer_workshop_CNV
bin_size = "10" # Can be varied

WORK = "/work/sduvarcall/G56-2016-Glioblastom/"
INPUT_BAM = WORK+"bam/"
OUTPUT_SEQZ = WORK+"sequenza_window"+base_window+"/sequenza_seqz/"
OUTPUT = WORK+"sequenza_window"+base_window+"/sequenza_bin"+bin_size+"/"

configfile: WORK+"samples.yaml"

ref_build = "hg38"

# Resources - paths inside docker
if ref_build == "GRCh37":
    # GRCh37 paths
    ref = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.fasta"
    wig = "/work/sduvarcall/bwa-0.7.13/reference/human_g1k_v37_decoy.gc"+base_window+"Base.wig.gz"
    chr_prefix = ""
else:
    # GRCh38 paths
    ref = "/work/sduvarcall/G37-2016-BRCAX-Genomes/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.fa"
    wig = "/work/sduvarcall/G37-2016-BRCAX-Genomes/reference_files_GRCh38/core_ref_GRCh38_hla_decoy_ebv/genome.gc"+base_window+"Base.wig.gz"
    chr_prefix = "chr"


# Get sample names from config file
SAMPLES = [sample for sample in config]
# SAMPLES = ["G56-sampleA1_nimblegen-medexome_HGL2LDSXX"]
print(SAMPLES)

onstart:
    shell("mkdir -p "+OUTPUT)

rule all:
    input:
        [expand(OUTPUT+"Sequenza_{tumor}_vs_{normal}_bin"+bin_size+"/Sequenza_{tumor}_vs_{normal}_chromosome_view.pdf", 
            tumor=config[s]["tumor"], normal=config[s]["normal"]) for s in SAMPLES],



###############################################################################
#### Sequenza
###############################################################################
rule gc_wiggle:
    input:
        ref = ref
    output:
        wig = wig
    conda: "envs/sequenza.yaml"
    shell:
        """
        sequenza-utils gc_wiggle --fasta {input.ref} -w {base_window} -o {output.wig}
        """

rule bam2seqz:
    input:
        normal=INPUT_BAM+"{normal}/{normal}.bam",
        tumor=INPUT_BAM+"{tumor}/{tumor}.bam",
        wig = wig
    output:
        seqz=OUTPUT_SEQZ+"{tumor}_vs_{normal}.seqz.gz"
    conda: "envs/sequenza.yaml"
    shell:
        """
        sequenza-utils bam2seqz -n {input.normal} -t {input.tumor} \
        -q 20 -N 20 -f "illumina" --fasta {ref} -gc {wig} -o {output}
        """

rule seqz_binning:
    input:
        seqz=OUTPUT_SEQZ+"{tumor}_vs_{normal}.seqz.gz"
    output:
        seqz_bin=OUTPUT+"{tumor}_vs_{normal}.seqz.bin{bin_size}.gz"
    conda: "envs/sequenza.yaml"
    shell:
        """
        sequenza-utils seqz_binning --seqz {input.seqz} -w {bin_size} -o {output.seqz_bin}
        """

rule sequenza_analysis:
    input:
        seqz_bin=OUTPUT+"{tumor}_vs_{normal}.seqz.bin{bin_size}.gz"
    output:
        res=OUTPUT+"Sequenza_{tumor}_vs_{normal}_bin{bin_size}/Sequenza_{tumor}_vs_{normal}_chromosome_view.pdf"
    conda: "envs/sequenza.yaml"
    shell:
        """
        Rscript -e '
        library(sequenza)
        seqzdata = sequenza.extract("{input.seqz_bin}", chromosome.list = paste0("{chr_prefix}",c(1:22,"X","Y")))

        CP = sequenza.fit(seqzdata, XY = c(X="chrX",Y="chrY"))
        
        id = "Sequenza_{wildcards.tumor}_vs_{wildcards.normal}"
        out.dir = paste0("{OUTPUT}", id, "_bin{wildcards.bin_size}")
        sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = id, out.dir=out.dir, XY = c(X="chrX",Y="chrY"))
        '
        """
