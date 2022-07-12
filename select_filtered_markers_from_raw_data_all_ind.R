library(qtl)
library(tidyr)
library(dplyr)
library(ggplot2)

# read new cross object after adding and organizing 
# newly-added markers
cross.allchrs.round3 <- 
  read.cross(file = "../round3/cross20.allchrs.round3-0.2_gen.csv",
             format = "csv",
             dir = "./",
             estimate.map = FALSE,
             genotypes = c("AA", "AB", "BB"))

cross.nocparam.allind <- read.cross("csv", dir = "./",
                    file = "batch_1.genotypes_1.merge_sire3492.sorted.transposed.genotyped_flipped.rqtl.csv",
                    estimate.map = FALSE, genotypes = c("AA","AB","BB"))

# read cross data with -c param after building the map
cross_withc <- read.csv("../round3/cross20.allchrs.round3-0.2_gen.csv",
         header = FALSE)
cross_withc <- t(cross_withc[1:3,-1])
colnames(cross_withc) <- c("marker_name", "scaffold", "pos")
cross_withc <- as.data.frame(cross_withc)
cross_withc[1:10,]

# read cross data without -c parameter, raw data
cross_noc <- read.csv("batch_1.genotypes_1.merge_sire3492.sorted.transposed.genotyped_flipped.rqtl.csv",
                      header = FALSE)
cross_noc_ids <- read.csv("batch_1.genotypes_1.merge_sire3492.sorted.transposed.genotyped_flipped.rqtl.csv",
                          header = FALSE)[,1]
cross_noc_ids <- as.character(cross_noc_ids[-1:-2])
head(cross_noc_ids)
cross_noc <- t(cross_noc[,-1])
colnames(cross_noc) <- c("marker_name", 
                         "scaffold",
                         cross.nocparam.allind$pheno$ID)
cross_noc <- as.data.frame(cross_noc)
cross_noc[1:10,1:10]

# filter data of cross without -c and filter it
cross_noc.filt <- left_join(cross_withc, 
                            cross_noc,
                            by = "marker_name")
cross_noc.filt[1:10,1:10]

# reformat the data of the cross without -c back to rqtl format
cross_noc.filt <- t(cross_noc.filt[, -4])
cross_noc.filt[1:10,1:10]

# export filtered data of the cross without -c parameter
write.csv(cross_noc.filt,
          "cross20.allchrs.round3-0.2_gen.nocparameterstacks.csv",
          row.names = TRUE,
          quote = FALSE)
