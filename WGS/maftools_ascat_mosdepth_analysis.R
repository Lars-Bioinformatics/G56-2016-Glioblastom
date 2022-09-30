# Setup R project using renv
setwd("/work/sduvarcall/G56-2016-Glioblastom/ascat_maftools")
renv::init()
Sys.setenv(RENV_PATHS_CACHE = '/work/sduvarcall/G56-2016-Glioblastom/ascat_maftools/renv/cache')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("remotes")
BiocManager::install("PoisonAlien/maftools")
BiocManager::install("GenomicRanges")
remotes::install_github(repo = 'VanLoo-lab/ascat/ASCAT')
install.packages("R.utils")
install.packages("tidyverse")

# Activate project
renv::activate(project = "/work/sduvarcall/G56-2016-Glioblastom/ascat_maftools")

library(ASCAT)
library(maftools)

project = "/work/sduvarcall/G56-2016-Glioblastom/"

tumors = paste0("G56-sample",
                # c("A1","B1","C1","D1","A2","B2","A3","B3","C3","D3","A4","B4","D4","E4","A8","B8","C8","A9","B9","C9","D9"),
                c("A9","B9","C9","D9"),
                "_truseq-nano-genome_HGL2LDSXX_",
                # c("S1","S2","S3","S4","S6","S7","S10","S11","S12","S13","S15","S16","S17","S18","S20","S31","S32","S35","S36","S37","S38"))
                c("S35","S36","S37","S38"))
# normals = c(rep("G56-blod-pt1_truseq-nano-genome_HGL2LDSXX_S5",4),rep("G56-blod-pt2_truseq-nano-genome_HGL2LDSXX_S9",2),
#             rep("G56-blod-pt3_truseq-nano-genome_HGL2LDSXX_S14",4),rep("G56-blod-pt4_truseq-nano-genome_HGL2LDSXX_S19",4),
#             rep("G56-blod-pt8_truseq-nano-genome_HGL2LDSXX_S34",3),rep("G56-blod-pt10_truseq-nano-genome_HGL2LDSXX_S34",4))
normals = rep("G56-blod-pt10_truseq-nano-genome_HGL2LDSXX_S39",4)
# Could not find solution for ploidy and purity
# D8,S33
# D9,S38

for (i in 1:length(tumors)){
  setwd("/work/sduvarcall/G56-2016-Glioblastom/ascat_maftools")
  
  tumor = tumors[i]
  normal = normals[i]
  
  dir.create(tumor, showWarnings = F)
  setwd(tumor)
  
  # counts = maftools::gtMarkers(t_bam = "../bam/G56-sampleA1_truseq-nano-genome_HGL2LDSXX_S1/G56-sampleA1_truseq-nano-genome_HGL2LDSXX_S1.bam",
  #                              n_bam = "../bam/G56-blod-pt1_truseq-nano-genome_HGL2LDSXX_S5/G56-blod-pt1_truseq-nano-genome_HGL2LDSXX_S5.bam",
  counts = maftools::gtMarkers(t_bam = paste0(project,"bam/",tumor,"/",tumor,".bam"),
                               n_bam = paste0(project,"bam/",normal,"/",normal,".bam"),
                               fa = "/work/sduvarcall/resources/hg38/Homo_sapiens_assembly38.fasta",
                               prefix = "chr",
                               build = "hg38",
                               loci = paste0(project,"ascat_maftools/GRCh38_SNP6.tsv.gz"),
                               nthreads = 32)
  
  # tum_counts = read.table("G56-sampleC3_truseq-nano-genome_HGL2LDSXX_S12_nucleotide_counts.tsv", header = T)
  # tum_counts = read.table("G56-sampleA1_truseq-nano-genome_HGL2LDSXX_S1_nucleotide_counts.tsv", header = T)
  
  
  # Generate BAF and LogR files
  ascat.bc = maftools::prepAscat(t_counts = paste0(tumor,"_nucleotide_counts.tsv"),
                                 n_counts = paste0(normal,"_nucleotide_counts.tsv"))
  # ascat.bc = maftools::prepAscat(t_counts = "G56-sampleA1_truseq-nano-genome_HGL2LDSXX_S1_nucleotide_counts.tsv",
  #                                n_counts = "G56-blod-pt1_truseq-nano-genome_HGL2LDSXX_S5_nucleotide_counts.tsv")
  # ascat.bc = maftools::prepAscat(t_counts = "G56-sampleC3_truseq-nano-genome_HGL2LDSXX_S12_nucleotide_counts.tsv",
  #                                n_counts = "G56-blod-pt3_truseq-nano-genome_HGL2LDSXX_S14_nucleotide_counts.tsv")
  
  # Run ASCAT
  ascat.bc = ASCAT::ascat.loadData(
    Tumor_LogR_file = paste0(tumor, "_nucleotide_counts.tumour.logR.txt"),
    Tumor_BAF_file = paste0(tumor, "_nucleotide_counts.tumour.BAF.txt"),
    Germline_LogR_file = paste0(tumor, "_nucleotide_counts.normal.logR.txt"),
    Germline_BAF_file = paste0(tumor, "_nucleotide_counts.normal.BAF.txt"),
    chrs = c(1:22, "X", "Y"),
    sexchromosomes = c("X", "Y")
  )
  
  ASCAT::ascat.plotRawData(ASCATobj = ascat.bc, img.prefix = tumor)
  ascat.bc = ASCAT::ascat.aspcf(ascat.bc)
  ASCAT::ascat.plotSegmentedData(ascat.bc)
  ascat.output = ASCAT::ascat.runAscat(ascat.bc) 
  
  # CBS segmentation
  maftools::segmentLogR(tumor_logR = paste0(tumor, "_nucleotide_counts.tumour.logR.txt"), sample_name = tumor)
}





# Mosdepth analysis
# library(tidyverse)
# p = plotMosdepth(
#   t_bed = "mosdepth/sampleA1.regions.bed.gz",
#   n_bed = "mosdepth/pt1_blod.regions.bed.gz",
#   segment = TRUE,
#   sample_name = "sampleA1"
# )
# # pdf(file = "sampleA1_cbs.pdf", height = 8, width = 16)
# png(file = "sampleA1_cbs.png", height = 2, width = 3)
# plotMosdepth(
#   t_bed = "mosdepth/sampleA1.regions.bed.gz",
#   n_bed = "mosdepth/pt1_blod.regions.bed.gz",
#   segment = TRUE,
#   sample_name = "sampleA1"
# )
# dev.off()
# ggsave(filename = "sampleA1_cbs.png", width = 16, height = 8)
# ggsave(filename = "sampleA1_cbs.jpg", width = 16, height = 8)
# ggsave(filename = "sampleA1_cbs.pdf", width = 16, height = 8)



