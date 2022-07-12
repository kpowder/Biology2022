library(qtl)
library(tidyr)
library(dplyr)
library(ggplot2)

# FUNCTIONS ---------------------------------------------------------------
chr_info <- function(cross, chr = "NC_036799.1"){
  # overall
  tmp <- pull.geno(cross, chr = chr)
  tmp.map <- pull.map(cross)
  overall <- apply(tmp, MARGIN = 2, 
                   function(x) length(which(is.na(x))))
  AA <- apply(tmp, MARGIN = 2, 
              function(x) round(length(which(x == 1))/length(which(!is.na(x))),4))
  AB <- apply(tmp, MARGIN = 2,
              function(x) round(length(which(x == 2))/length(which(!is.na(x))),4))
  BB <- apply(tmp, MARGIN = 2, 
              function(x) round(length(which(x == 3))/length(which(!is.na(x))),4))
  df <- rbind(names(tmp.map[[chr]]),
              as.numeric(tmp.map[[chr]]), 
              overall, AA, AB, BB)
  df <- t(df) %>% 
    as.data.frame() %>% 
    mutate(LG = gsub("_[0-9]*$","",row.names(.)),
           bp = as.numeric(gsub(".*_([0-9]*$)","\\1",row.names(.))))
  colnames(df) <- c("ID", "pos", "missing", "AA","AB","BB","LG", "bp")
  df <- df %>% select(ID, pos, LG, bp, missing, AA, AB, BB)
  df$pos <- as.character(df$pos) %>% as.numeric()
  df$bp <- as.character(df$bp) %>% as.numeric()
  df$missing <- as.character(df$missing) %>% as.numeric()
  df$AA <- as.character(df$AA) %>% as.numeric()
  df$AB <- as.character(df$AB) %>% as.numeric()
  df$BB <- as.character(df$BB) %>% as.numeric()
  return(df)
}
detect_outlier <- function(x, start.win){
  iqr.vec <- IQR(x)
  third_quantile <- quantile(x)[4]
  return(which(x > third_quantile+(iqr.vec*1.5)) + (start.win-1))
}
SlidingWindow.outliers <- function(data, win = 20, steps = 10){
  total <- length(data)
  spots <- seq(from = 1, to = (total - win), by = steps)
  result <- list()
  for (i in 1:length(spots)) {
    result[[i]] <- detect_outlier(
      data[spots[i]:(spots[i] + win - 1)], start.win = spots[i])
  }
  return(unique(unlist(result)))
}


# chr1 --------------------------------------------------------------------

chr <- 1

# read cross object with chr info after cleaning-round2
load(paste("../round2/cross20.chr", 
           chr,
           ".rmOutliers.rmRM.round2.Rdata",
           sep = ""))
cross20.chr.rmOutliers.rmRM

# the offset of this chr is not 0. These lines fixes the problem.
newmap <- est.map(cross20.chr.rmOutliers.rmRM, 
                  error.prob = 0.005,
                  offset = 0)
cross20.chr.rmOutliers.rmRM <- 
  replace.map(cross20.chr.rmOutliers.rmRM, newmap)
summaryMap(newmap)

par(mfrow = c(1,2))
plotRF(cross20.chr.rmOutliers.rmRM, what = "rf")

# get max genetic position in the chr
#max_position <- max(cross20.chr.rmOutliers.rmRM$geno$NC_036794.1$map)
#max_position

# read markers to add
#toadd.df <- read.table("../round2/missassembled_markers/from17tochr15.txt",
#                       header = FALSE,
#                       row.names = 1)
#toadd.df[,1:10]

# add markers one by one
cross20.chr.rmOutliers.rmRM.addmarkers <- cross20.chr.rmOutliers.rmRM
#for(i in 1:nrow(toadd.df)){
#  cross20.chr.rmOutliers.rmRM.addmarkers <-
#    addmarker(cross20.chr.rmOutliers.rmRM.addmarkers, 
#              as.numeric(toadd.df[i,]),
#              rownames(toadd.df)[i],
#              chr = chrnames(cross20.chr.rmOutliers.rmRM.addmarkers)[1],
#              pos = max_position+(i*0.1))
#}

#cross20.chr.rmOutliers.rmRM.addmarkers <- 
#  est.rf(cross20.chr.rmOutliers.rmRM.addmarkers)
#plotRF(cross20.chr.rmOutliers.rmRM.addmarkers, what = "rf")

# export cross
write.cross(cross = cross20.chr.rmOutliers.rmRM.addmarkers,
            format = "csvs",
            filestem = paste("cross20.chr", 
                             chr,
                             ".rmOutliers.rmRM.addmarkers", 
                             sep = ""))


# read new cross object after adding and organizing 
# newly-added markers
cross20.chr.rmOutliers.rmRM.addmarkers.organized <- 
  read.cross(file = paste0("cross20.chr",
                           chr,
                           ".round3_gen.organized-v0.2.csv"),
             format = "csv",
             dir = "./",
             estimate.map = FALSE,
             genotypes = c("AA", "AB", "BB"))

cross20.chr.rmOutliers.rmRM.addmarkers.organized <- 
  est.rf(cross20.chr.rmOutliers.rmRM.addmarkers.organized)

plotRF(cross20.chr.rmOutliers.rmRM.addmarkers.organized, what = "rf")

# estimate map in chr with newly added and organized markers
summaryMap(cross20.chr.rmOutliers.rmRM)
newmap <- est.map(cross20.chr.rmOutliers.rmRM.addmarkers.organized, 
                  error.prob = 0.005,
                  offset = 0)

# read new cross object after adding and organizing 
# newly-added markers
cross20.chr.rmOutliers.rmRM.addmarkers.organized <- 
  read.cross(file = paste0("cross20.chr",
                           chr,
                           ".rmOutliers.rmRM.addmarkers_gen.organized.csv"),
             format = "csv",
             dir = "./",
             estimate.map = FALSE,
             genotypes = c("AA", "AB", "BB"))

cross20.chr.rmOutliers.rmRM.addmarkers.organized <- 
  est.rf(cross20.chr.rmOutliers.rmRM.addmarkers.organized)

plotRF(cross20.chr.rmOutliers.rmRM.addmarkers.organized, what = "rf")

# estimate map in chr with newly added and organized markers
summaryMap(cross20.chr.rmOutliers.rmRM)
newmap <- est.map(cross20.chr.rmOutliers.rmRM.addmarkers.organized, 
                  error.prob = 0.005,
                  offset = 0)
cross20.chr.rmOutliers.rmRM.addmarkers.organized <- 
  replace.map(cross20.chr.rmOutliers.rmRM.addmarkers.organized, newmap)
summaryMap(newmap)

plotMap(cross20.chr.rmOutliers.rmRM)
plotMap(cross20.chr.rmOutliers.rmRM.addmarkers.organized)

# export cross object with newly-added and organized markers, and 
# re-estimated map
write.cross(cross = cross20.chr.rmOutliers.rmRM.addmarkers.organized,
            format = "csvs",
            filestem = paste("cross20.chr", 
                             chr,
                             ".round3.3", 
                             sep = ""))

cross20.chr.rmOutliers.rmRM.addmarkers.organized <- 
  replace.map(cross20.chr.rmOutliers.rmRM.addmarkers.organized, newmap)
summaryMap(newmap)

plotMap(cross20.chr.rmOutliers.rmRM)
plotMap(cross20.chr.rmOutliers.rmRM.addmarkers.organized)

tmp.chr_info.rmOut.rmRM <- chr_info(cross20.chr.rmOutliers.rmRM, 
                                    chr = chrnames(cross)[chr])
ggplot(tmp.chr_info.rmOut.rmRM, aes(x = bp, y = pos)) +
  geom_point()

tmp.chr_info.rmOut.rmRM.organized <- 
  chr_info(cross20.chr.rmOutliers.rmRM.addmarkers.organized, 
                                    chr = chrnames(cross)[chr])
ggplot(tmp.chr_info.rmOut.rmRM.organized, aes(x = bp, y = pos)) +
  geom_point()

# export cross object with newly-added and organized markers, and 
# re-estimated map
write.cross(cross = cross20.chr.rmOutliers.rmRM.addmarkers.organized,
            format = "csvs",
            filestem = paste("cross20.chr", 
                             chr,
                             ".round3-0.2", 
                             sep = ""))
