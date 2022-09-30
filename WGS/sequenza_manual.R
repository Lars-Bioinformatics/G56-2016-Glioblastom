# setRepositories(graphics = FALSE, ind = 1:6)
# install.packages("sequenza")

setwd("/work/sduvarcall/G56-2016-Glioblastom/sequenza_default/")

library(sequenza)

seqzdata <- sequenza.extract("G56-sampleA1_truseq-nano-genome_HGL2LDSXX_S1_vs_G56-blod-pt1_truseq-nano-genome_HGL2LDSXX_S5.seqz.gz",
                                   chromosome.list = paste0("chr",c(1:22,"X","Y")))

CP <- sequenza.fit(seqzdata)

dir.create("../sequenza_default/", showWarnings = F)
sequenza.results(sequenza.extract = seqzdata,
                 cp.table = CP, 
                 cellularity = 0.8,
                 ploidy = 2, female = T, XY = c("chrX","chrY"),
                 sample.id = "A1", out.dir = "../sequenza_manual/")

# Grid search
cp.plot(CP)
cp.plot.contours(CP, add = TRUE,
                 likThresh = c(0.999, 0.95),
                 col = c("lightsalmon", "red"), pch = 20)

# Chromosome view
chromosome.view(mut.tab = seqzdata$mutations[[1]], baf.windows = seqzdata$BAF[[1]],
                ratio.windows = seqzdata$ratio[[1]],  min.N.ratio = 2,
                segments = seqzdata$segments[[1]], 
                main = seqzdata$chromosomes[1], CNn = 2,
                cellularity = 0.8, ploidy = 6,
                avg.depth.ratio = 0.5)
