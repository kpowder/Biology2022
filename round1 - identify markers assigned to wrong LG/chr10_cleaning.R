library(qtl)
library(tidyr)
library(dplyr)
library(ggplot2)
#library(rcompanion)
#library(evobiR)

# This code uses 4 new functions: chr_info, detect_outlier,
# SlidingWindow.outliers, and dropRedundantMarkers. The code
# for the last one is in a separate file.
# 
# Additionally, this code requires to have a cross object
# (already filtered for missing data) loaded in the environment.
# The second section of this file, RF symmetrical matrix,
# creates a symmetrical matrix with the recombination frequency
# for any pair of markers across all chromosomes. Since
# this calculation takes a few minutes, I will only include
# this piece of code on the file for chr1.
# 
# The purpose of this code is to perform a first pass of
# marker cleaning in order to build an accurate genetic map.
# This first pass consists of the following steps:
# 1) Identify and remove markers whose recombination frequency
# profile does not agree to that of the rest of the markers
# in the chromosome.
# 2) Classify the discordant markers in two different categories:
# a) misassembled markers (markers that belong to a different chr),
# and b) outlier markers (markers that do not belong to the
# current or any other chr).
# 3) Remove redundant markers (consecutive markers that have 
# the same genotypic information and only differ in missing data).
# 
# All the markers that are removed in this first round of
# cleaning are stored in different files, depending on the
# step at which they were removed.


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
SlidingWindow.outliers <- function(data, win = 20, steps =10){
  total <- length(data)
  spots <- seq(from = 1, to = (total - win), by = steps)
  result <- list()
  for (i in 1:length(spots)) {
    result[[i]] <- detect_outlier(
      data[spots[i]:(spots[i] + win - 1)], start.win = spots[i])
  }
  return(unique(unlist(result)))
}

# RF symmetrical matrix ---------------------------------------------------
# RF are stored in rqtl in a matrix where the top-left triangle of the matrix
# contains the RF values and the bottom-right triangle of the matrix contains
# the LOD scores for a statistical test that evaluates whether two markers segragate
# independently or not.
# 
# The following commands extract the RF values from this matrix and makes a new
# symmetrical matrix with these values.

#rf_allchrs <- cross$rf
#diag(rf_allchrs) <- 1
#rf_allchrs[upper.tri(rf_allchrs, diag = FALSE)] <- NA

#for(i in 1:(ncol(rf_allchrs)-1)){
#  for(j in (i+1):nrow(rf_allchrs)){
#    rf_allchrs[i,j] <- rf_allchrs[j,i]
#  }
#}

#rf_allchrs.df <- data.frame(index=1:ncol(rf_allchrs), 
#                         median = apply(rf_allchrs, 1, median), 
#                         sd = apply(rf_allchrs, 1, sd),
#                         names = colnames(rf_allchrs))



# chr10 --------------------------------------------------------------------

chr <- 10

load("Robjects_round1_genetic_map_building_LcLt3492.RData")

cross20.chr <- subset(cross, chr = chrnames(cross)[chr])
#cross20.chr$geno$NC_036780.1$map %>% names %>% View

## Create pairwise RF plots of chr 1 with chr1:22
# Plots showing RF values
#pdf(paste("RFplots_362ind-maxmissingdata-0.2.chr", 
#          chr, "vsAll.rf.pdf", sep = ""))
#for(i in 1:22){
#  plotRF(cross, what = "rf", chr = c(chrnames(cross)[chr], chrnames(cross)[i]))
#}
#dev.off()

# Plots showing LOD scores for Ho: rf = 1/2
#pdf(paste("RFplots_362ind-maxmissingdata-0.2.chr",
#          chr, "vsAll.lod.pdf", sep = ""))
#for(i in 1:22){
#  plotRF(cross, what = "lod", chr = c(chrnames(cross)[chr], chrnames(cross)[i]))
#}
#dev.off()

# RF are stored in rqtl in a matrix where the top-left triangle of the matrix
# contains the RF values and the bottom-right triangle of the matrix contains
# the LOD scores for a statistical test that evaluates whether two markers segragate
# independently or not.
# 
# The following commands extract the RF values from this matrix and makes a new
# symmetrical matrix with these values.
rf_chr <- cross20.chr$rf
diag(rf_chr) <- 1
rf_chr[upper.tri(rf_chr, diag = FALSE)] <- NA

for(i in 1:(ncol(rf_chr)-1)){
  for(j in (i+1):nrow(rf_chr)){
    rf_chr[i,j] <- rf_chr[j,i]
  }
}

# calculating the median and sd of each row of the rf_chr matrix
# and saving it in a data-frame
rf_chr.df <- data.frame(index=1:ncol(rf_chr), 
                        median = apply(rf_chr, 1, median), 
                        sd = apply(rf_chr, 1, sd),
                        names = colnames(rf_chr))

# adding outlier column, ussing my own function to detect them
rf_chr.df$outliers <- 0
rf_chr.df$outliers[SlidingWindow.outliers(rf_chr.df$median)] <- 1

#ggplot(rf_chr.df, aes(x = index, y = median)) +
#  geom_point() +
#  geom_text(aes(label = ifelse(outliers==1, as.character(index),'')), )

#ggplot(rf_chr.df, aes(x = index, y = median)) +
#  geom_point() +
#  geom_text(aes(label = as.character(index)), size = 3, hjust = 2)

# filter outliers manually (the sliding window can
# be somewhat stringent, it is important to compare the RF plots
# and the sliding window plot and see in what markers they agree
# or don't agree)
SlidingWindow.outliers(rf_chr.df$median)
rf_chr.df$outliers <- 0
rf_chr.df$outliers[c(26, 95, 97, 134, 135,
                     167, 174, 190)] <- 1

# asking whether the outliers belong to a different chromosome or not
#par(mfrow = c(3,4))
#for(marker in rf_chr.df[rf_chr.df$outliers == 1, "names"]){
#  mean_rf <- numeric()
#  for(i in 1:22){
#    tmp <- rf_allchrs[marker,grep(pattern = chrnames(cross)[i],
#                                  x = colnames(rf_allchrs))]
#    mean_rf[i] <- (mean(tmp))
#    print(marker)
#    print(mean_rf)
#  }
#  plot(x = c(1:20,22,23), y = mean_rf, 
#       main = paste("INDEX:", rf_chr.df[rf_chr.df$names == marker, "index"]))
#}

mean_rf.df <- data.frame(marker = "x", index = 0, chr1 = 0, chr2 = 0, chr3 = 0, chr4 = 0,
                         chr5 = 0, chr6 = 0, chr7 = 0, chr8 = 0, chr9 = 0, chr10 = 0,
                         chr11 = 0, chr12 = 0, chr13 = 0, chr14 = 0, chr15 = 0, chr16 = 0,
                         chr17 = 0, chr18 = 0, chr19 = 0, chr20 = 0, chr22 = 0, chr23 = 0, 
                         stringsAsFactors = FALSE)

par(mfrow = c(1,1))
for(marker in rf_chr.df[rf_chr.df$outliers == 1, "names"]){
  mean_rf <- numeric()
  for(i in 1:22){
    tmp <- rf_allchrs[marker,grep(pattern = chrnames(cross)[i],
                                  x = colnames(rf_allchrs))]
    mean_rf[i] <- (mean(tmp))
    print(marker)
    print(mean_rf)
  }
  mean_rf.df <- rbind(mean_rf.df, c(marker, rf_chr.df[rf_chr.df$names == marker, "index"],
                                    mean_rf))
}

# reshape data frame with mean RF's of outlier markers
mean_rf.df <- gather(data = mean_rf.df, 
                     key = "chr_number", 
                     value = "mean_rf", -marker, -index) %>%
  filter(marker != "x") %>% 
  filter(chr_number != paste("chr", chr, sep = ""))

# make a boxplot of mean RF's on each chromosome for outlier markers
#mean_rf.df %>% 
#  ggplot(aes(x = index, y = as.numeric(mean_rf)), data = .) +
#  geom_boxplot()

# markers that have a mean RF less than 0.4 with the markers of a different 
# chromosome are likely to be missassembled.
# These commands save the info of those markers in the missassembled_markers dir.
missassembled_markers <- mean_rf.df[mean_rf.df$mean_rf < 0.4,]
print("These are the detected misassembled markers:")
print(missassembled_markers)

#for(i in 1:nrow(missassembled_markers)){
#  tmp.out <- cross20.chr$geno[[1]]$data[,missassembled_markers$marker[i]]
#  tmp.out <- c(missassembled_markers[i, "marker"], tmp.out)
#  write.table(t(tmp.out),
#              file = paste("missassembled_markers\\from", chr, "to", missassembled_markers$chr[i], ".txt", sep = ""),
#              append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
#}

# The info of the remaining outlier markers are then saved in the outlier_markers dir
tmp <- rf_chr.df[rf_chr.df$outliers == 1, "names"]
outlier_markers <- as.character(tmp[!(tmp %in% missassembled_markers$marker)])

print(paste("There are ", length(outlier_markers), " remaining outlier markers", sep = ""))
#for(i in 1:length(outlier_markers)){
#  tmp.out <- cross20.chr$geno[[1]]$data[,outlier_markers[i]]
#  tmp.out <- c(outlier_markers[i], tmp.out)
#  write.table(t(tmp.out),
#              file = paste("outlier_markers\\", chr, ".txt", sep = ""),
#              append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
#}

# save outliers stats
tmp.out <- chr_info(cross20.chr, chr = chrnames(cross)[chr])[rf_chr.df[rf_chr.df$outliers == 1, "index"],]
#write.table(tmp.out,
#            file = paste("all_outlier_markers_chr", chr, "_round1.txt", sep = ""),
#            quote = FALSE, sep = "\t",
#            row.names = TRUE, col.names = TRUE)

### drop all outlier markers from the chromosome
cross20.chr.rmOutliers <- cross20.chr
cross20.chr.rmOutliers
todrop <- as.character(rf_chr.df[rf_chr.df$outliers == 1, "names"])
cross20.chr.rmOutliers <- drop.markers(cross20.chr.rmOutliers, todrop)
cross20.chr.rmOutliers

# estimate map in chr without outliers
newmap <- est.map(cross20.chr.rmOutliers, error.prob=0.005)
cross20.chr.rmOutliers <- replace.map(cross20.chr.rmOutliers, newmap)
summaryMap(newmap)

# plot Map
par(mfrow=c(1,2))
plotMap(newmap, main = paste("CHR: ", chr, " after removing outliers", sep = ""))

### drop all redundant markers from the chromosome
output.dRM <- dropRedundantMarkers(cross20.chr.rmOutliers)
cross20.chr.rmOutliers.rmRM <- drop.markers(cross20.chr.rmOutliers, output.dRM[[1]])

# estimate map in chr without outliers
newmap <- est.map(cross20.chr.rmOutliers.rmRM, error.prob = 0.005)
cross20.chr.rmOutliers.rmRM <- replace.map(cross20.chr.rmOutliers.rmRM, newmap)
summaryMap(newmap)

# plot Map
plotMap(newmap, main = paste("CHR: ", chr, " after removing outliers and RM", sep = ""))

# save redundant markers (RM) stats
tmp.out <- data.frame(todrop = names(cross20.chr$geno[[1]]$map)[as.numeric(output.dRM[[2]]$todrop)], 
                      delegate = names(cross20.chr$geno[[1]]$map)[as.numeric(output.dRM[[2]]$delegate)])
write.table(tmp.out,
            file = paste("redundant_markers_chr", chr, "_round1.txt", sep = ""),
            quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# generate summary plots
par(mfrow=c(1,1))
rf <- pull.rf(cross20.chr.rmOutliers.rmRM)
lod <- pull.rf(cross20.chr.rmOutliers.rmRM, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
plotRF(cross20.chr.rmOutliers.rmRM)

# saving chr info after removing outliers and RM in a cross object
save(cross20.chr.rmOutliers.rmRM, file = paste("cross20.chr", chr, ".rmOutliers.rmRM.Rdata", sep = ""))
