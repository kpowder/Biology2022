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


# MAIN --------------------------------------------------------------------
# 
# read new cross object after adding and organizing 
# newly-added markers
cross.allchrs.allind <- 
  read.cross(file = "../round4.2/cross20.allchrs.round3-0.2_gen.withcparameterstacks.csv",
             format = "csv",
             dir = "./",
             estimate.map = FALSE,
             genotypes = c("AA", "AB", "BB"))

# look for problematic individuals
plot(countXO(cross.allchrs.allind), ylab="Number of crossovers")

#cross.allchrs.allind <- subset(cross.allchrs.allind,
#                               ind=(countXO(cross.allchrs.allind) < 350))

# estimate map again
#newmap <- est.map(cross.allchrs.allind,
#                  error.prob=0.005, 
#                  offset = 0)
#cross.allchrs.allind <- replace.map(cross.allchrs.allind, newmap)
summaryMap(cross.allchrs.allind)

# look for genotyping errors
cross.allchrs.allind <- calc.errorlod(cross.allchrs.allind, 
                                      error.prob=0.010)

print(toperr <- top.errorlod(cross.allchrs.allind, cutoff=3))
write.csv(toperr, "genotyping_errors.csv", quote = FALSE)

i <- 9
tmp.ind <- toperr$id[toperr$chr==chrnames(cross.allchrs.allind)[i]]
tmp.ind <- as.character(tmp.ind)
plotGeno(cross.allchrs.allind,
         chr = chrnames(cross.allchrs.allind)[i],
         ind = tmp.ind,
         cutoff = 8, 
         include.xo = FALSE)

dim(toperr)

cross.allchrs.allind.clean <- cross.allchrs.allind
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- as.character(toperr$id[i])
  id_index <- which(cross.allchrs.allind$pheno$ID == id)
  mar <- toperr$marker[i]
  cross.allchrs.allind.clean$geno[[chr]]$data[id_index, mar] <- NA
}

plot(countXO(cross.allchrs.allind.clean), ylab="Number of crossovers")

i <- 9
tmp.ind <- toperr$id[toperr$chr==chrnames(cross.allchrs.allind.clean)[i]]
tmp.ind <- as.character(tmp.ind)
plotGeno(cross.allchrs.allind,
         chr = chrnames(cross.allchrs.allind.clean)[i],
         ind = tmp.ind,
         cutoff = 3, 
         include.xo = FALSE)

# remove the genotypes with error lod >= 3
dim(toperr)
round6.filt.clean <- round6.filt
for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  print(i)
  round6.filt.clean$geno[[chr]]$data[id, mar] <- NA
}

# estimate genotyping error rate
loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075,
                   0.01, 0.0125, 0.015, 0.0175, 0.02)
for(i in seq(along=err)) {
  cat(i, "of", length(err), "\n")
  tempmap <- est.map(cross.allchrs.allind,
                     error.prob=err[i])
  loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)

plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02),
     ylab=expression(paste(log[10], " likelihood")))

pdf("genotyping error rate.pdf")
plot(err[1:6], lod[1:6], xlab="Genotyping error rate", xlim=c(0,0.02),
     ylab=expression(paste(log[10], " likelihood")))
dev.off()



par(mfrow = c(1,2), las = 1)
plot(ntyped(cross.allchrs.allind.clean), 
     ylab = "No. typed markers", 
     main = "No. genotypes by individual", 
     ylim = c(0,3200))
plot(ntyped(cross.allchrs.allind.clean, "mar"), 
     ylab = "No. typed individuals", 
     main = "No. genotypes by marker",
     ylim = c(0,360))

# estimate final map
pdf("recombination_plots_map_cross20-LcLt3492.pdf")
par(mfrow=c(1,1))
for(i in c(1:22))
{
  plotRF(cross.allchrs.allind.clean,
         alternate.chrid = TRUE, 
         chr = chrnames(cross.allchrs.allind.clean)[i],
         what = "rf")
}
dev.off()

tmp.newmap <- est.map(cross.allchrs.allind.clean,
                      offset = 0,
                      error.prob = 0.010,
                      maxit = 10000,
                      verbose = TRUE)
cross.clean <- replace.map(cross.allchrs.allind.clean, tmp.newmap)

par(mfrow=c(1,1))
summary(cross.clean)
summaryMap(cross.clean)
plotMap(cross.clean)
plotMissing(cross.clean)

# export cross object with newly-added and organized markers, and 
# re-estimated map
write.cross(cross = cross.clean,
            format = "csvs",
            filestem = "LcLtF2.cross20.finalmap-v0.2.041320")

pdf("dotplot_cross20.allchrs.round3-0.2_gen.withcparameterstacks.pdf")
for(chr in 1:22){
print("hello!")
tmp.chr_info.cross_clean <- chr_info(cross.clean, 
         chr = chrnames(cross.clean)[chr])
ggplot(data = tmp.chr_info.cross_clean,
       aes(x=bp, y=pos)) + geom_point()
}
dev.off()
