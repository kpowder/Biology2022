library(qtl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rcompanion)
library(evobiR)

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

chr <- 1
rf_allchrs <- cross$rf
diag(rf_allchrs) <- 1
rf_allchrs[upper.tri(rf_allchrs, diag = FALSE)] <- NA

for(i in 1:(ncol(rf_allchrs)-1)){
  for(j in (i+1):nrow(rf_allchrs)){
    rf_allchrs[i,j] <- rf_allchrs[j,i]
  }
}

#rf_allchrs.df <- data.frame(index=1:ncol(rf_allchrs), 
#                         median = apply(rf_allchrs, 1, median), 
#                         sd = apply(rf_allchrs, 1, sd),
#                         names = colnames(rf_allchrs))


# chr1 --------------------------------------------------------------------

cross20.chr1 <- subset(cross, chr = chrnames(cross)[1])
#cross20.chr1$geno$NC_036780.1$map %>% names %>% View

## Create pairwise RF plots of chr 1 with chr1:22
# Plots showing RF values
pdf("RFplots_362ind-maxmissingdata-0.2.chr1vsAll.rf.pdf")
for(i in 1:22){
  plotRF(cross, what = "rf", chr = c(chrnames(cross)[1], chrnames(cross)[i]))
}
dev.off()

# Plots showing LOD scores for Ho: rf = 1/2
pdf("RFplots_362ind-maxmissingdata-0.2.chr1vsAll.lod.pdf")
for(i in 1:22){
  plotRF(cross, what = "lod", chr = c(chrnames(cross)[1], chrnames(cross)[i]))
}
dev.off()

# RF are stored in rqtl in a matrix where the top-left triangle of the matrix
# contains the RF values and the bottom-right triangle of the matrix contains
# the LOD scores for a statistical test that evaluates whether two markers segragate
# independently or not.
# 
# The following commands extract the RF values from this matrix and makes a new
# symmetrical matrix with these values.
rf_chr1 <- cross20.chr1$rf
diag(rf_chr1) <- 1
rf_chr1[upper.tri(rf_chr1, diag = FALSE)] <- NA

for(i in 1:(ncol(rf_chr1)-1)){
  for(j in (i+1):nrow(rf_chr1)){
    rf_chr1[i,j] <- rf_chr1[j,i]
  }
}

rf_chr1.df <- data.frame(index=1:ncol(rf_chr1), 
                         median = apply(rf_chr1, 1, median), 
                         sd = apply(rf_chr1, 1, sd),
                         names = colnames(rf_chr1))

# adding outlier column, ussing my own function to detect them
rf_chr1.df$outliers <- 0
rf_chr1.df$outliers[SlidingWindow.outliers(rf_chr1.df$median)] <- 1

ggplot(rf_chr1.df, aes(x = index, y = median)) +
  geom_point() +
  geom_text(aes(label = ifelse(outliers==1, as.character(index),'')), )

# asking whether the outliers belong to a different chromosome or not

par(mfrow = c(3,4))
for(marker in rf_chr1.df[rf_chr1.df$outliers == 1, "names"]){
  mean_rf <- numeric()
  for(i in 1:22){
    tmp <- rf_allchrs[marker,grep(pattern = chrnames(cross)[i],
                                x = colnames(rf_allchrs))]
    mean_rf[i] <- (mean(tmp))
    print(marker)
    print(mean_rf)
  }
  plot(x = c(1:20,22,23), y = mean_rf, 
       main = paste("INDEX:", rf_chr1.df[rf_chr1.df$names == marker, "index"]))
}

mean_rf.df <- data.frame(marker = "x", index = 0, chr1 = 0, chr2 = 0, chr3 = 0, chr4 = 0,
           chr5 = 0, chr6 = 0, chr7 = 0, chr8 = 0, chr9 = 0, chr10 = 0,
           chr11 = 0, chr12 = 0, chr13 = 0, chr14 = 0, chr15 = 0, chr16 = 0,
           chr17 = 0, chr18 = 0, chr19 = 0, chr20 = 0, chr22 = 0, chr23 = 0, 
           stringsAsFactors = FALSE)
par(mfrow = c(1,1))
for(marker in rf_chr1.df[rf_chr1.df$outliers == 1, "names"]){
  mean_rf <- numeric()
  for(i in 1:22){
    tmp <- rf_allchrs[marker,grep(pattern = chrnames(cross)[i],
                                  x = colnames(rf_allchrs))]
    mean_rf[i] <- (mean(tmp))
    print(marker)
    print(mean_rf)
  }
  mean_rf.df <- rbind(mean_rf.df, c(marker, rf_chr1.df[rf_chr1.df$names == marker, "index"],
        mean_rf))
}

# reshape data frame with mean RF's of outlier markers
mean_rf.df <- gather(data = mean_rf.df, 
       key = "chr", 
       value = "mean_rf", -marker, -index) %>%
  filter(marker != "x") %>% 
  filter(chr != paste("chr", chr, sep = ""))

# make a boxplot of mean RF's on each chromosome for outlier markers
mean_rf.df %>% 
  ggplot(aes(x = index, y = as.numeric(mean_rf)), data = .) +
  geom_boxplot()

# markers that have a mean RF less than 0.4 with the markers of a different 
# chromosome are likely to be missassembled.
# These commands save the info of those markers in the missassembled_markers dir.
missassembled_markers <- mean_rf.df[mean_rf.df$mean_rf < 0.4,]

for(i in 1:nrow(missassembled_markers)){
  tmp.out <- cross20.chr1$geno$NC_036780.1$data[,missassembled_markers$marker[i]]
  tmp.out <- c(missassembled_markers[i, "marker"], tmp.out)
  write.table(t(tmp.out),
            file = paste("missassembled_markers\\from", chr, "to", missassembled_markers$chr[i], ".txt", sep = ""),
            append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}

# The info of the remaining outlier markers are then saved in the outlier_markers dir
tmp <- rf_chr1.df[rf_chr1.df$outliers == 1, "names"]
outlier_markers <- as.character(tmp[!(tmp %in% missassembled_markers$marker)])

for(i in 1:length(outlier_markers)){
  tmp.out <- cross20.chr1$geno$NC_036780.1$data[,outlier_markers[i]]
  tmp.out <- c(outlier_markers[i], tmp.out)
  write.table(t(tmp.out),
              file = paste("outlier_markers\\", chr, ".txt", sep = ""),
              append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}

# save outliers stats
tmp.out <- chr_info(cross20.chr1, chr = chrnames(cross)[1])[rf_chr1.df[rf_chr1.df$outliers == 1, "index"],]
write.table(tmp.out,
            file = paste("filtered_markers_chr", chr, "_round1.txt", sep = ""),
            quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = TRUE)

# Saving a bunch of cleaned chromosomes as small cross objects to be merged later
cross20.chr1.cleaned <- cross20.chr1
cross20.chr1.cleaned

# drop all outlier markers from the chromosome
todrop <- as.character(rf_chr1.df[rf_chr1.df$outliers == 1, "names"])
cross20.chr1.cleaned <- drop.markers(cross20.chr1.cleaned, todrop)
cross20.chr1.cleaned

# plot LOD scores vs estimated RF for all marker pairs
rf <- pull.rf(cross20.chr1.cleaned)
lod <- pull.rf(cross20.chr1.cleaned, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

# estimate map in clean chr
newmap <- est.map(cross20.chr1.cleaned, error.prob=0.005)
#cross20.chr1.cleaned <- replace.map(cross20.chr1.cleaned, newmap)

# print map summary and generate summary plots
summaryMap(newmap)
plotMap(newmap)
plotRF(cross20.chr1.cleaned)

# prune redundant markers
output.prunning <- dropSimilarMarkers(cross = cross20.chr1.cleaned,
                                      rf.threshold = 1e-10,
                                      sd.weight = 0,
                                      blockSize =400, 
                                      keepEnds = TRUE,
                                      verbose = TRUE)
a <- est.map(output.prunning, error.prob = 0.005)
output.prunning <- replace.map(output.prunning, a)
summaryMap(newmap)
summaryMap(a)
par(mfrow=c(1,2))
plotMap(newmap)
plotMap(a)
names_newmap <- cross20.chr1.cleaned$geno$NC_036780.1$map %>% names
names_a <- output.prunning$geno$NC_036780.1$map %>% names
setdiff(names_newmap, names_a)

chr_i



test <- dropRedundantMarkers(cross20.chr1.cleaned)
output.redundant <- drop.markers(cross20.chr1.cleaned, test[[2]])
a <- est.map(output.redundant, error.prob = 0.005)
output.redundant <- replace.map(output.redundant, a)
summaryMap(newmap)
summaryMap(a)nfo(output.prunning, chrnames(output.prunning)[1])

par(mfrow=c(1,2))
plotMap(newmap)
plotMap(a)





output.prunning <- dropSimilarMarkers_modified(cross = cross20.chr1.cleaned, 
                            sd.weight = 0,
                            keepEnds = TRUE,
                            verbose = TRUE)

# saving the cleaned small cross object
save(cross20.chr1.cleaned, file = "cross20.chr1.cleaned.Rdata")

a <- dropSimilarMarkers(cross20.chr1.cleaned,
                   chr = chrnames(cross20.chr1.cleaned)[1],
                   rf.threshold = 1e-10,
                   blockSize = 400)

## remove duplicates and missing data > 20?
## make a diagnostic plot to make sure the filtering
## is not too stringengt.
## 
## Maybe make a scatterplot: on the x-axis the bp 
## of the marker and in the y-axis the distance 
## with respect to the next marker
## 
## We could do this before and after prunning
## and color code the filtering steps: before/after

## I think at this point, prunning by hand may be
## the best option. At most I could use
## removeDuplicates.






colnames(rf_chr1) %>% grep("NC_036780.1_3064708", .)
colnames(rf_chr1)[5]    #NC_036780.1_664890
colnames(rf_chr1)[18]   #NC_036780.1_3064708
colnames(rf_chr1)[21]   #NC_036780.1_3580165
colnames(rf_chr1)[26]   #NC_036780.1_4499331
colnames(rf_chr1)[152]  #NC_036780.1_25936185
colnames(rf_chr1)[183]  #NC_036780.1_30033678

chr_info(cross20.chr1, chr = chrnames(cross20.chr1)[1])


mean_rf <- numeric()
for(i in 1:22){
  tmp <- cross$rf[colnames(rf_chr1)[5],grep(pattern = chrnames(cross)[i],
                                     x = colnames(cross$rf))]
  mean_rf[i] <- (mean(tmp))
}

plot(mean_rf)

colnames(rf_chr1)[5]

max(rf_chr1[,1], na.rm = TRUE)
rf_chr1[1,1]
View(rf_chr1)
ncol(rf_chr1)-1)










females_indexes <- read.table("females_LcLt.txt")[,1]
cross$pheno$sex <- 1
cross$pheno$sex[cross$pheno$ID %in% females_indexes] <- 0

tmp <- data.frame(ID = as.numeric(names(ntyped(cross))),
                  ntyped = ntyped(cross),
                  sex = cross$pheno$sex,
                    sequencer = "gray30", stringsAsFactors = FALSE)
tmp$sequencer[tmp$ID < 1272] <- "gray30"
tmp$sequencer[tmp$ID >= 1272] <- "blue"
#tmp %>% 
ggplot(tmp, aes(x = ID, y = ntyped, color = as.factor(sex))) +
  geom_point(size = 3) +
  geom_text(aes(label = ifelse(ntyped<4500 & as.factor(sex)==0, 
                               as.character(ID),'')), size = 20) +
  theme_light()
