# script to process popoolation stats
setwd("~/Desktop/FST")
library(tidyr)
library(vegan)
library(fossil)
# read the main output file and set population labels
FST_main <- read.csv("population_first.fst", 
                     sep = '\t', col.names = 
                       c("Scaffold", "Window", "nSNPs", "fraction_covered", "coverage_ave", "N25_N75", "N25_N125", "N25_N175", 
                         "N25_N225", "N25_S25", "N25_S75", "N25_S125", "N25_S175", "N25_S225", 
                         "N25_UZ", "N75_N125", "N75_N175", "N75_N225", "N75_S25", "N75_S75", "N75_S125", 
                         "N75_S175", "N75_S225", "N75_UZ", "N125_N175", "N125_N225", "N125_S25", "N125_S75", 
                         "N125_S125", "N125_S175", "N125_S225", "N125_UZ", "N175_N225", "N175_S25", 
                         "N175_S75", "N175_S125", "N175_S175", "N175_S225", "N175_UZ", "N225_S25", 
                          "N225_S75", "N225_S125", "N225_S175", "N225_S225", "N225_UZ", "S25_S75", 
                         "S25_S125", "S25_S175", "S25_S225", "S25_UZ", "S75_S125", 
                         "S75_S175", "S75_S225", "S75_UZ", "S125_S175", "S125_S225", "S125_UZ",
                         "S175_S225", "S175_UZ", "S225_UZ"))
z_chr_scaffolds <- c('scaffold4435|size4078', 'scaffold751|size275148', 'scaffold35|size722381', 'scaffold572|size163719', 'scaffold740|size235087', 'scaffold1546|size33109', 'scaffold247|size627634', 'scaffold196|size868841', 'scaffold779|size92251', 'scaffold2358|size14287', 'scaffold600|size417131', 'scaffold216|size461870', 'scaffold390|size298204', 'scaffold364|size585662', 'scaffold1895|size21730', 'scaffold409|size222884', 'scaffold169|size242750', 'scaffold444|size384115', 'scaffold325|size886211', 'scaffold348|size845032', 'scaffold865|size221853', 'scaffold15|size673834', 'scaffold546|size598492', 'scaffold218|size453784', 'scaffold645|size139123', 'scaffold1044|size193178', 'scaffold1660|size63095', 'scaffold581|size122005', 'scaffold1237|size50895', 'scaffold1246|size194460', 'scaffold1756|size66030', 'scaffold1161|size202793', 'scaffold155|size350779', 'scaffold282|size321509', 'scaffold518|size233254', 'scaffold789|size91266', 'scaffold502|size248588', 'scaffold3599|size6641', 'scaffold804|size458571', 'scaffold1040|size105542', 'scaffold275|size1004024', 'scaffold532|size371787', 'scaffold221|size590889', 'scaffold140|size987430', 'scaffold1038|size208510', 'scaffold598|size233855', 'scaffold48|size965936', 'scaffold1486|size80768', 'scaffold137|size260670', 'scaffold193|size439059', 'scaffold50|size611115v02', 'scaffold87|size517236', 'scaffold229|size374464', 'bab_gene_region_CCv02', 'scaffold728|size179971v02', 'scaffold330|size554422', 'scaffold673|size283851', 'scaffold812|size88987', 'scaffold1046|size116232', 'scaffold1204|size193440', 'scaffold1618|size30317', 'scaffold412|size197301', 'scaffold820|size139877')
#cleanup table 
A <- function(x) gsub(".*=","",x)
FST_m <- cbind(FST_main[1:5], apply(FST_main[6:60],2, A))
B <- function(x) as.numeric(as.character(x))
FST_m <- cbind(FST_main[1:5], apply(FST_main[6:60],2, A))
# write output file for later use
write.table(FST_m, "FST_by_pop", sep = '\t', quote = FALSE)

# if restarting analysis, begin here
FST_m <- read.table("FST_by_pop", sep = '\t', header = TRUE)
#filter out Z and autosomes
FST_m_z <- FST_m[FST_m$Scaffold %in% z_chr_scaffolds, ]
FST_m_auto <- FST_m[!FST_m$Scaffold %in% z_chr_scaffolds, ]
fst_out_auto <- matrix(nrow = 9, ncol = 9)
fst_out_z <- matrix(nrow = 9, ncol = 9)

# filtering steps, can be skipped if not filtering
FST_m_z <- FST_m_z[FST_m_z$nSNPs > 20, ]
FST_m_z <- FST_m_z[FST_m_z$fraction_covered > 0.20, ]
FST_m_auto <- FST_m_auto[FST_m_auto$nSNPs > 20, ]
FST_m_auto <- FST_m_auto[FST_m_auto$fraction_covered > 0.20, ]
#rename columns
colnames(fst_out_auto) <- c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ")
row.names(fst_out_auto) <- c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ")
#write renamed and possibly filtered columns
write.table(fst_out_auto, "~/Desktop/fst_out_autosomes.csv", quote = FALSE, sep = '\t')
write.table(fst_out_z, "~/Desktop/fst_out_z_linked.csv", quote = FALSE, sep = '\t')
#clunky script to reorganize fst matrix
diag(fst_out_auto) <- 0
fst_out_auto[1, 2] <- mean(FST_m_auto[, 33])
fst_out_auto[1, 3] <- mean(FST_m_auto[, 26])
fst_out_auto[1, 4] <- mean(FST_m_auto[, 18])
fst_out_auto[1, 5] <- mean(FST_m_auto[, 40])
fst_out_auto[1, 6] <- mean(FST_m_auto[, 42])
fst_out_auto[1, 7] <- mean(FST_m_auto[, 43])
fst_out_auto[1, 8] <- mean(FST_m_auto[, 44])
fst_out_auto[1, 9] <- mean(FST_m_auto[, 45])
fst_out_auto[2, 3] <- mean(FST_m_auto[, 25])
fst_out_auto[2, 4] <- mean(FST_m_auto[, 17])
fst_out_auto[2, 5] <- mean(FST_m_auto[, 34])
fst_out_auto[2, 6] <- mean(FST_m_auto[, 36])
fst_out_auto[2, 7] <- mean(FST_m_auto[, 37])
fst_out_auto[2, 8] <- mean(FST_m_auto[, 38])
fst_out_auto[2, 9] <- mean(FST_m_auto[, 39])
fst_out_auto[3, 4] <- mean(FST_m_auto[, 16])
fst_out_auto[3, 5] <- mean(FST_m_auto[, 27])
fst_out_auto[3, 6] <- mean(FST_m_auto[, 29])
fst_out_auto[3, 7] <- mean(FST_m_auto[, 30])
fst_out_auto[3, 8] <- mean(FST_m_auto[, 31])
fst_out_auto[3, 9] <- mean(FST_m_auto[, 32])
fst_out_auto[4, 5] <- mean(FST_m_auto[, 19])
fst_out_auto[4, 6] <- mean(FST_m_auto[, 21])
fst_out_auto[4, 7] <- mean(FST_m_auto[, 22])
fst_out_auto[4, 8] <- mean(FST_m_auto[, 23])
fst_out_auto[4, 9] <- mean(FST_m_auto[, 24])
fst_out_auto[5, 6] <- mean(FST_m_auto[, 47])
fst_out_auto[5, 7] <- mean(FST_m_auto[, 48])
fst_out_auto[5, 8] <- mean(FST_m_auto[, 49])
fst_out_auto[5, 9] <- mean(FST_m_auto[, 50])
fst_out_auto[6, 7] <- mean(FST_m_auto[, 55])
fst_out_auto[6, 8] <- mean(FST_m_auto[, 56])
fst_out_auto[6, 9] <- mean(FST_m_auto[, 57])
fst_out_auto[7, 8] <- mean(FST_m_auto[, 58])
fst_out_auto[7, 9] <- mean(FST_m_auto[, 59])
fst_out_auto[8, 9] <- mean(FST_m_auto[, 60])
# more clunky script to rebind columns
fst_s_auto <- cbind((FST_m_auto[, 33]), (FST_m_auto[, 26]), (FST_m_auto[, 18]), (FST_m_auto[, 40]), (FST_m_auto[, 42]), (FST_m_auto[, 43]), (FST_m_auto[, 44]), (FST_m_auto[, 45]), (FST_m_auto[, 25]), (FST_m_auto[, 17]), (FST_m_auto[, 34]), (FST_m_auto[, 36]), (FST_m_auto[, 37]), (FST_m_auto[, 38]), (FST_m_auto[, 39]), (FST_m_auto[, 16]), (FST_m_auto[, 27]), (FST_m_auto[, 29]), (FST_m_auto[, 30]), (FST_m_auto[, 31]), (FST_m_auto[, 32]), (FST_m_auto[, 19]), (FST_m_auto[, 21]), (FST_m_auto[, 22]), (FST_m_auto[, 23]), (FST_m_auto[, 24]), (FST_m_auto[, 47]), (FST_m_auto[, 48]), (FST_m_auto[, 49]), (FST_m_auto[, 50]), (FST_m_auto[, 55]), (FST_m_auto[, 56]), (FST_m_auto[, 57]), (FST_m_auto[, 58]), (FST_m_auto[, 59]), (FST_m_auto[, 60]))
fst_noUZ_auto <- cbind(FST_m_auto[, 1:5], (FST_m_auto[, 33]), (FST_m_auto[, 26]), (FST_m_auto[, 18]), (FST_m_auto[, 40]), (FST_m_auto[, 42]), (FST_m_auto[, 43]), (FST_m_auto[, 44]), (FST_m_auto[, 25]), (FST_m_auto[, 17]), (FST_m_auto[, 34]), (FST_m_auto[, 36]), (FST_m_auto[, 37]), (FST_m_auto[, 38]), (FST_m_auto[, 16]), (FST_m_auto[, 27]), (FST_m_auto[, 29]), (FST_m_auto[, 30]), (FST_m_auto[, 31]), (FST_m_auto[, 19]), (FST_m_auto[, 21]), (FST_m_auto[, 22]), (FST_m_auto[, 23]), (FST_m_auto[, 47]), (FST_m_auto[, 48]), (FST_m_auto[, 49]), (FST_m_auto[, 55]), (FST_m_auto[, 56]), (FST_m_auto[, 58]))
fst_UZ_auto <- cbind(FST_m_auto[, 1:5], (FST_m_auto[, 45]), (FST_m_auto[, 39]), (FST_m_auto[, 32]), (FST_m_auto[, 24]), (FST_m_auto[, 50]), (FST_m_auto[, 57]), (FST_m_auto[, 59]), (FST_m_auto[, 60]))
fst_UZ_z <- cbind(FST_m_z[, 1:5], (FST_m_z[, 45]), (FST_m_z[, 39]), (FST_m_z[, 32]), (FST_m_z[, 24]), (FST_m_z[, 50]), (FST_m_z[, 57]), (FST_m_z[, 59]), (FST_m_z[, 60]))

# with UZ filtered out, fix scaffold names
fst_noUZ_auto <- separate(data = fst_noUZ_auto, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
fst_noUZ_auto$Scaffold<-gsub("v02","",as.character(fst_noUZ_auto$Scaffold))
fst_noUZ_auto$length<-gsub("v02","",as.character(fst_noUZ_auto$length))
fst_noUZ_auto$length<-replace_na(fst_noUZ_auto$length,"414000")
colnames(fst_noUZ_auto) <- c("Scaffold", "Length", "Window", "nSNPs", "fraction_covered", "coverage_ave",
                         "N75_N125", "N75_N175", "N75_N225", "N75_S25", "N75_S125", 
                         "N75_S175", "N75_S225", "N125_N175", "N125_N225", "N125_S25", 
                         "N125_S125", "N125_S175", "N125_S225", "N175_N225", "N175_S25", 
                          "N175_S125", "N175_S175", "N175_S225", "N225_S25", 
                          "N225_S125", "N225_S175", "N225_S225", 
                         "S25_S125", "S25_S175", "S25_S225", 
                        "S125_S175", "S125_S225",
                         "S175_S225")
#optional filtering
fst_noUZ_auto_filtered <- fst_noUZ_auto[fst_noUZ_auto$nSNPs >  20, ]
fst_noUZ_auto_filtered <- fst_noUZ_auto_filtered[fst_noUZ_auto_filtered$fraction_covered >  0.2, ]
#some statistics and plotting
hist(as.numeric(as.matrix(fst_noUZ_auto[, 7:34])))
mean(as.numeric(as.matrix(fst_noUZ_auto[, 7:34])))
mean(as.numeric(as.matrix(fst_noUZ_auto_filtered[, 7:34])))
pdf("auto_fst.pdf", height = 4, width = 4)
hist(as.numeric(as.matrix(fst_noUZ_auto_filtered[, 7:34])), xlim = c(0,1), col = "black", main = "Distribution of autosomal FST values", xlab = "FST")
dev.off()

pdf("auto_fst_zoomed.pdf", height = 4, width = 6)
hist(as.numeric(as.matrix(fst_noUZ_auto_filtered[, 7:34])), xlim = c(0.4, 0.6), ylim = c(0,40), col = 'black', main = "", xlab = "", ylab = "", breaks = 50)
dev.off()
#means and t tests
mean(as.matrix(fst_noUZ_auto_filtered[, 7:34]))
mean(as.matrix(fst_noUZ_z_filtered[, 7:34]))

t.test(as.matrix(fst_noUZ_auto_filtered[, 7:34]))
t.test(as.matrix(fst_noUZ_z_filtered[, 7:34]))


##z chromosome uses repeated process of above but filtered to Z scaffolds



fst_out_auto<- Matrix::forceSymmetric(fst_out_auto)
fst_out_auto <- as.matrix(fst_out_auto)
mean(fst_out_auto[1:8, 1:8])
mean(fst_out_auto[9, 1:8])
# Z
# sorting and indexing
FST_sorted <- FST_m[order(-FST_m$N225_S225),]
FST_m$observation <- 1:nrow(FST_m) 
plot(FST_m$observation, FST_m$N25_N75, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_m$observation, FST_m$N25_N125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_m$observation, FST_m$N25_N175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_m$observation, FST_m$N25_N225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_m$observation, FST_m$N25_S25, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_m$observation, FST_m$N25_S75, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_m$observation, FST_m$N25_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_m$observation, FST_m$N25_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_m$observation, FST_m$N25_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_m$observation, FST_m$N25_UZ, cex = 0.2, pch = 19, ylim = c(0,1))


#moving average function
ma <- function(x, n = 10){stats::filter(x, rep(1 / n, n), sides = 2)}
#filtering out period scaffold
FST_period <- FST_m[which(FST_m$Scaffold == "scaffold532"), ]
#apply moving average function
FST_period$N225_UZ_ma <- ma(FST_period$N225_UZ)
FST_period$N175_UZ_ma <- ma(FST_period$N175_UZ)
FST_period$N125_UZ_ma <- ma(FST_period$N125_UZ)
FST_period$N75_UZ_ma <- ma(FST_period$N75_UZ)
FST_period$S25_UZ_ma <- ma(FST_period$S25_UZ)
FST_period$S125_UZ_ma <- ma(FST_period$S125_UZ)
FST_period$S175_UZ_ma <- ma(FST_period$S175_UZ)
FST_period$S225_UZ_ma <- ma(FST_period$S225_UZ)
# plotting individual populations and their moving averages
pdf("period_fst_N225_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$N225_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_period$Window, FST_period$N225_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("period_fst_N175_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$N175_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_period$Window, FST_period$N175_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("period_fst_N125_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$N125_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_period$Window, FST_period$N125_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("period_fst_N75_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$N75_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_period$Window, FST_period$N75_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("period_fst_N25_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$N25_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
par(op)
dev.off()
pdf("period_fst_S25_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$S25_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_period$Window, FST_period$S25_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("period_fst_S75_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$S75_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
par(op)
dev.off()
pdf("period_fst_S125_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$S125_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_period$Window, FST_period$S125_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("period_fst_S175_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$S175_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_period$Window, FST_period$S175_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("period_fst_S225_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_period$Window, FST_period$S225_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_period$Window, FST_period$S225_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()

# add codon annotations
rect(75057, 0.9, 75082, 1 , col = 'green')
rect(77479, 16, 77748, 17 , col = 'green')
rect(101450, 0.9, 101463, 1 , col = 'green')
rect(101464, 0.9, 101551, 1 , col = 'black')
rect(103643, 0.9, 103728, 1 , col = 'black')
rect(104601, 0.9, 104754, 1 , col = 'black')
rect(110704, 0.9, 110801, 1 , col = 'black')
rect(110887, 0.9, 110970, 1 , col = 'black')
rect(112874, 0.9, 113006, 1 , col = 'black')
rect(113155, 0.9, 113214, 1 , col = 'black')
rect(113215, 0.9, 113231, 1 , col = 'black')
rect(114358, 0.9, 114563, 1 , col = 'black')
rect(114924, 0.9, 115155, 1 , col = 'black')
rect(115522, 0.9, 115636, 1 , col = 'black')
rect(115734, 0.9, 115858, 1 , col = 'black')
rect(116369, 0.9, 116521, 1 , col = 'black')
rect(125722, 0.9, 125837, 1 , col = 'black')
rect(125922, 0.9, 126020, 1 , col = 'black')
rect(128758, 0.9, 128908, 1 , col = 'black')
rect(130219, 0.9, 130411, 1 , col = 'black')
rect(132239, 0.9, 132332, 1 , col = 'black')
rect(133776, 0.9, 133915, 1 , col = 'black')
rect(134917, 0.9, 135031, 1 , col = 'black')
rect(135932, 0.9, 136111, 1 , col = 'black')
rect(138505, 0.9, 138647, 1 , col = 'black')
rect(139391, 0.9, 139498, 1 , col = 'black')
rect(140310, 0.9, 140515, 1 , col = 'black')

dev.off()
#calculate means within genic region
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N225_UZ)
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N175_UZ)
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N125_UZ)
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N75_UZ)
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N25_UZ)
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S25_UZ)
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S75_UZ)
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S125_UZ)
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S175_UZ)
mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S225_UZ)

period_fst <- c(mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N225_UZ), mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N175_UZ), mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N125_UZ), mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N75_UZ), mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$N25_UZ), mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S25_UZ), mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S75_UZ), mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S125_UZ), mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S175_UZ), mean(FST_period[FST_period$Window > 75000 & FST_period$Window < 145000, ]$S225_UZ))
#calculate means outside genic region

mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N225_UZ)
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N175_UZ)
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N125_UZ)
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N75_UZ)
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N25_UZ)
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S25_UZ)
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S75_UZ)
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S125_UZ)
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S175_UZ)
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S225_UZ)

period_out_fst <- c(mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N225_UZ),
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N175_UZ),
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N125_UZ),
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N75_UZ),
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$N25_UZ),
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S25_UZ),
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S75_UZ),
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S125_UZ),
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S175_UZ),
mean(FST_period[FST_period$Window < 75000 | FST_period$Window > 145000, ]$S225_UZ)
)

# show genic region
rect(114000, 0.9, 115000, 1, col = 'blue')

FST_pdfr <- FST_m[which(FST_m$Scaffold == "scaffold87"), ]

FST_pdfr$N225_UZ_ma <- ma(FST_pdfr$N225_UZ)
FST_pdfr$N175_UZ_ma <- ma(FST_pdfr$N175_UZ)
FST_pdfr$N125_UZ_ma <- ma(FST_pdfr$N125_UZ)
FST_pdfr$N75_UZ_ma <- ma(FST_pdfr$N75_UZ)
FST_pdfr$S25_UZ_ma <- ma(FST_pdfr$S25_UZ)
FST_pdfr$S125_UZ_ma <- ma(FST_pdfr$S125_UZ)
FST_pdfr$S175_UZ_ma <- ma(FST_pdfr$S175_UZ)
FST_pdfr$S225_UZ_ma <- ma(FST_pdfr$S225_UZ)
# repeated steps as above but for pdfr scaffold
plot(FST_pdfr$Window, FST_pdfr$N225_S25, cex = 0.2, pch = 19, ylim = c(0,1))
pdf("pdfr_fst.pdf", height = 6, width = 8)
plot(FST_pdfr$Window, FST_pdfr$N225_UZ, cex = 0.5, pch = 19, ylim = c(0,1), ylab = "Fst", xlab = "Location on scaffold", main = "pdfr, N225 vs UZ")
plot(FST_pdfr$Window, FST_pdfr$N175_UZ, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N125_UZ, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N75_UZ, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N25_UZ, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S25_UZ, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S75_UZ, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S125_UZ, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S175_UZ, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S225_UZ, cex = 0.2, pch = 19, ylim = c(0,1))

rect(44493, 0.9, 46265, 1, col = 'black')
rect(46266, 0.9, 46349, 0.97, col = 'black')
rect(47065, 0.9, 47184, 0.97, col = 'black')
rect(52728, 0.9, 52847, 0.97, col = 'black')
rect(53062, 0.9, 53180, 0.97, col = 'black')
rect(54524, 0.9, 54707, 0.97, col = 'black')
rect(58394, 0.9, 58479, 0.97, col = 'black')
rect(60223, 0.9, 60292, 0.97, col = 'black')
rect(60605, 0.9, 60834, 0.97, col = 'black')
rect(61975, 0.9, 62125, 0.97, col = 'black')
rect(62228, 0.9, 62308, 0.97, col = 'black')
rect(66266, 0.9, 66411, 0.97, col = 'black')
rect(79859, 0.9, 80216, 0.97, col = 'green')
rect(82804, 0.9, 83222, 0.97, col = 'green')

dev.off()

mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N225_UZ)
mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N175_UZ)
mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N125_UZ)
mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N75_UZ)
mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N25_UZ)
mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S25_UZ)
mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S75_UZ)
mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S125_UZ)
mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S175_UZ)
mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S225_UZ)
#comparative statistics
pdfr_fst<- c(mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N225_UZ), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N175_UZ), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N125_UZ), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N75_UZ), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N25_UZ), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S25_UZ), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S75_UZ), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S125_UZ), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S175_UZ), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S225_UZ))
pdfr_fst_range<- c(mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N175_N225), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N125_N225), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N75_N225), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N225_S25), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N225_S125), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N225_S175), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N225_S225), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N125_N175), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N75_N175), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N175_S25), mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N175_S125),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N175_S175),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N175_S225),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N75_N125),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N125_S25),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N125_S125),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N125_S175),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N125_S225),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N75_S25),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N75_S125),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N75_S175),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$N75_S225),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S25_S125),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S25_S175),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S25_S225),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S125_S175),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S125_S225),  mean(FST_pdfr[FST_pdfr$Window > 44000 & FST_pdfr$Window < 83000, ]$S175_S225))
period_fst_range<- c(mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N175_N225), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N125_N225), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N75_N225), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N225_S25), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N225_S125), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N225_S175), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N225_S225), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N125_N175), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N75_N175), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N175_S25), mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N175_S125),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N175_S175),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N175_S225),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N75_N125),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N125_S25),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N125_S125),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N125_S175),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N125_S225),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N75_S25),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N75_S125),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N75_S175),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$N75_S225),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$S25_S125),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$S25_S175),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$S25_S225),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$S125_S175),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$S125_S225),  mean(FST_period[FST_period$Window > 44000 & FST_period$Window < 83000, ]$S175_S225))

mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N225_UZ)
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N175_UZ)
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N125_UZ)
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N75_UZ)
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N25_UZ)
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S25_UZ)
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S75_UZ)
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S125_UZ)
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S175_UZ)
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S225_UZ)
#pairwise comparative plots for pdfr
plot(FST_pdfr$Window, FST_pdfr$N175_N225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N125_N225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N75_N225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N225_S25, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N225_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N225_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N225_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N125_N175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N75_N175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N175_S25, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N175_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N175_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N175_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N75_N125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N125_S25, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N125_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N125_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N125_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N75_S25, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N75_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N75_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$N75_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S25_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S25_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S25_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S125_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S125_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_pdfr$Window, FST_pdfr$S175_S225, cex = 0.2, pch = 19, ylim = c(0,1))

#comparative pairwise plots for period
plot(FST_period$Window, FST_period$N175_N225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N125_N225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N75_N225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N225_S25, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N225_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N225_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N225_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N125_N175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N75_N175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N175_S25, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N175_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N175_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N175_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N75_N125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N125_S25, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N125_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N125_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N125_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N75_S25, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N75_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N75_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$N75_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$S25_S125, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$S25_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$S25_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$S125_S175, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$S125_S225, cex = 0.2, pch = 19, ylim = c(0,1))
plot(FST_period$Window, FST_period$S175_S225, cex = 0.2, pch = 19, ylim = c(0,1))

#covert to numeric
test_vector <- abs(as.numeric(ord_pdfr_z_all[, 36 ]) - as.numeric(ord_pdfr_z_all[, 37]))

# set up pdfr matrix
pdfr_frame <- as.matrix(cbind(ord_pdfr_z_all[, 2], test_vector))
pdfr_frame <- as.data.frame(pdfr_frame)
pdfr_frame$V1 <- as.numeric(as.character(pdfr_frame$V1))
pdfr_frame$test_vector <- as.numeric(as.character(pdfr_frame$test_vector))

plot(pdfr_frame$V1, pdfr_frame$test_vector, pch = 19, cex = 0.5)
#UZ comparisons
pdfr_out_fst <- c(mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N225_UZ),
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N175_UZ),
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N125_UZ),
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N75_UZ),
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$N25_UZ),
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S25_UZ),
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S75_UZ),
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S125_UZ),
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S175_UZ),
mean(FST_pdfr[FST_pdfr$Window < 44000 | FST_pdfr$Window > 83000, ]$S225_UZ))



#means and t tests
mean(pdfr_fst[c(1:4, 6, 8:10)])
mean(period_fst[c(1:4, 6, 8:10)])
mean(pdfr_out_fst[c(1:4, 6, 8:10)])
mean(period_out_fst[c(1:4, 6, 8:10)])

t.test(pdfr_fst[c(1:4, 6, 8:10)], pdfr_out_fst[c(1:4, 6, 8:10)])
t.test(period_fst[c(1:4, 6, 8:10)], period_out_fst[c(1:4, 6, 8:10)])

t.test(pdfr_fst[c(1:4, 6, 8:10)], period_fst[c(1:4, 6, 8:10)])
t.test(pdfr_out_fst[c(1:4, 6, 8:10)], period_out_fst[c(1:4, 6, 8:10)])


#### pi
#filter individual populations for pi
N225_pi <- read.table("~/Desktop/pi/N225_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
N225_pi <- separate(data = N225_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
N225_pi$Scaffold<-gsub("v02","",as.character(N225_pi$Scaffold))
N225_pi$length<-gsub("v02","",as.character(N225_pi$length))
N225_pi$length<-replace_na(N225_pi$length,"414000")
head(N225_pi)

N225_pi_period <- N225_pi[which(N225_pi$Scaffold == "scaffold532"), ]
N225_pi_pdfr <- N225_pi[which(N225_pi$Scaffold == "scaffold87"), ]

plot(N225_pi_period$Position, N225_pi_period$Pi, pch = 19, cex = .3)
plot(N225_pi_pdfr$Position, N225_pi_pdfr$Pi, pch = 19, cex = .3)

# repeat for other populatons
#comparisons of genic vs exgenic regions for pdfr
UZ1_exgenic_pi <- mean(UZ1_pi_pdfr[UZ1_pi_pdfr$Position < 44000 | UZ1_pi_pdfr$Position > 83000, ]$Pi, na.rm = TRUE  )
UZ1_genic_pi <- mean(UZ1_pi_pdfr[UZ1_pi_pdfr$Position > 44000 & UZ1_pi_pdfr$Position < 83000, ]$Pi, na.rm = TRUE  )
N225_exgenic_pi <- mean(N225_pi_pdfr[N225_pi_pdfr$Position < 44000 | N225_pi_pdfr$Position > 83000, ]$Pi, na.rm = TRUE  )
N225_genic_pi <- mean(N225_pi_pdfr[N225_pi_pdfr$Position > 44000 & N225_pi_pdfr$Position < 83000, ]$Pi, na.rm = TRUE  )
N175_exgenic_pi <- mean(N175_pi_pdfr[N175_pi_pdfr$Position < 44000 | N175_pi_pdfr$Position > 83000, ]$Pi, na.rm = TRUE  )
N175_genic_pi <- mean(N175_pi_pdfr[N175_pi_pdfr$Position > 44000 & N175_pi_pdfr$Position < 83000, ]$Pi, na.rm = TRUE  )
N125_exgenic_pi <- mean(N125_pi_pdfr[N125_pi_pdfr$Position < 44000 | N125_pi_pdfr$Position > 83000, ]$Pi, na.rm = TRUE  )
N125_genic_pi <- mean(N125_pi_pdfr[N125_pi_pdfr$Position > 44000 & N125_pi_pdfr$Position < 83000, ]$Pi, na.rm = TRUE  )
N75_exgenic_pi <- mean(N75_pi_pdfr[N75_pi_pdfr$Position < 44000 | N75_pi_pdfr$Position > 83000, ]$Pi, na.rm = TRUE  )
N75_genic_pi <- mean(N75_pi_pdfr[N75_pi_pdfr$Position > 44000 & N75_pi_pdfr$Position < 83000, ]$Pi, na.rm = TRUE  )
S25_exgenic_pi <- mean(S25_pi_pdfr[S25_pi_pdfr$Position < 44000 | S25_pi_pdfr$Position > 83000, ]$Pi, na.rm = TRUE  )
S25_genic_pi <- mean(S25_pi_pdfr[S25_pi_pdfr$Position > 44000 & S25_pi_pdfr$Position < 83000, ]$Pi, na.rm = TRUE  )
S125_exgenic_pi <- mean(S125_pi_pdfr[S125_pi_pdfr$Position < 44000 | S125_pi_pdfr$Position > 83000, ]$Pi, na.rm = TRUE  )
S125_genic_pi <- mean(S125_pi_pdfr[S125_pi_pdfr$Position > 44000 & S125_pi_pdfr$Position < 83000, ]$Pi, na.rm = TRUE  )
S175_exgenic_pi <- mean(S175_pi_pdfr[S175_pi_pdfr$Position < 44000 | S175_pi_pdfr$Position > 83000, ]$Pi, na.rm = TRUE  )
S175_genic_pi <- mean(S175_pi_pdfr[S175_pi_pdfr$Position > 44000 & S175_pi_pdfr$Position < 83000, ]$Pi, na.rm = TRUE  )
S225_exgenic_pi <- mean(S225_pi_pdfr[S225_pi_pdfr$Position < 44000 | S225_pi_pdfr$Position > 83000, ]$Pi, na.rm = TRUE  )
S225_genic_pi <- mean(S225_pi_pdfr[S225_pi_pdfr$Position > 44000 & S225_pi_pdfr$Position < 83000, ]$Pi, na.rm = TRUE  )

pdfr_genic_pi <- c(N225_genic_pi, N175_genic_pi, N125_genic_pi, N75_genic_pi, S25_genic_pi, S125_genic_pi, S175_genic_pi, S225_genic_pi)
pdfr_exgenic_pi <- c(N225_exgenic_pi, N175_exgenic_pi, N125_exgenic_pi, N75_exgenic_pi, S25_exgenic_pi, S125_exgenic_pi, S175_exgenic_pi, S225_exgenic_pi)
mean(pdfr_genic_pi)
t.test(pdfr_genic_pi)
mean(pdfr_exgenic_pi)
t.test(pdfr_exgenic_pi)
t.test(pdfr_genic_pi, pdfr_exgenic_pi)

#repeat for period


## tajima's D

UZ_D <- read.table("~/Desktop/D/UZ1_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))

UZ_D_period <- UZ_D[UZ_D$Scaffold == 'scaffold532|size371787', ]
plot(UZ_D_period$Position, UZ_D_period$D)
UZ_D_pdfr <- UZ_D[UZ_D$Scaffold == 'scaffold87|size517236', ]
plot(UZ_D_pdfr$Position, UZ_D_pdfr$D, pch = 19, cex = 0.7)

N225_D <- read.table("~/Desktop/D/N225_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))
N175_D <- read.table("~/Desktop/D/N175_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))
N125_D <- read.table("~/Desktop/D/N125_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))
N75_D <- read.table("~/Desktop/D/N75_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))
N25_D <- read.table("~/Desktop/D/N25_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))
S25_D <- read.table("~/Desktop/D/S25_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))
S75_D <- read.table("~/Desktop/D/S175_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))
S125_D <- read.table("~/Desktop/D/S125_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))
S175_D <- read.table("~/Desktop/D/S75_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))
S225_D <- read.table("~/Desktop/D/S225_D1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "D"))



N225_D_period <- N225_D[N225_D$Scaffold == 'scaffold532|size371787', ]
plot(N225_D_period$Position, N225_D_period$D, pch = 19, cex = 0.7)
N225_D_pdfr <- N225_D[N225_D$Scaffold == 'scaffold87|size517236', ]
pdf("N225_D.pdf", height = 4, width = 6)
plot(N225_D_pdfr$Position, N225_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3), ylab = "Tajima's D", xlab = "Position on scaffold", main = "N225, Tajima's D")
rect(44493, 2.7, 46265, 2.9, col = 'black')
rect(46266, 2.7, 46349, 2.9, col = 'black')
rect(47065, 2.7, 47184, 2.9, col = 'black')
rect(52728, 2.7, 52847, 2.9, col = 'black')
rect(53062, 2.7, 53180, 2.9, col = 'black')
rect(54524, 2.7, 54707, 2.9, col = 'black')
rect(58394, 2.7, 58479, 2.9, col = 'black')
rect(60223, 2.7, 60292, 2.9, col = 'black')
rect(60605, 2.7, 60834, 2.9, col = 'black')
rect(61975, 2.7, 62125, 2.9, col = 'black')
rect(62228, 2.7, 62308, 2.9, col = 'black')
rect(66266, 2.7, 66411, 2.9, col = 'black')
rect(79859, 2.7, 80216, 2.9, col = 'green')
rect(82804, 2.7, 83222, 2.9, col = 'green')
dev.off()
# repeat for other populations and for period below

N225_D_period <- N225_D[N225_D$Scaffold == 'scaffold532|size371787', ]

pdf("N225_D_period.pdf", height = 4, width = 6)
plot(N225_D_period$Position, N225_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3), ylab = "Tajima's D", xlab = "Position on scaffold", main = "N225, Tajima's D")

rect(75057, 2.7, 75082, 2.9, col = 'green')
rect(77479, 2.7, 77748, 2.9, col = 'green')
rect(101450, 2.7, 101463, 2.9, col = 'green')
rect(101464, 2.7, 101551, 2.9, col = 'black')
rect(103643, 2.7, 103728, 2.9, col = 'black')
rect(104601, 2.7, 104754, 2.9, col = 'black')
rect(110704, 2.7, 110801, 2.9, col = 'black')
rect(110887, 2.7, 110970, 2.9, col = 'black')
rect(112874, 2.7, 113006, 2.9, col = 'black')
rect(113155, 2.7, 113214, 2.9, col = 'black')
rect(113215, 2.7, 113231, 2.9, col = 'black')
rect(114358, 2.7, 114563, 2.9, col = 'black')
rect(114924, 2.7, 115155, 2.9, col = 'black')
rect(115522, 2.7, 115636, 2.9, col = 'black')
rect(115734, 2.7, 115858, 2.9, col = 'black')
rect(116369, 2.7, 116521, 2.9, col = 'black')
rect(125722, 2.7, 125837, 2.9, col = 'black')
rect(125922, 2.7, 126020, 2.9, col = 'black')
rect(128758, 2.7, 128908, 2.9, col = 'black')
rect(130219, 2.7, 130411, 2.9, col = 'black')
rect(132239, 2.7, 132332, 2.9, col = 'black')
rect(133776, 2.7, 133915, 2.9, col = 'black')
rect(134917, 2.7, 135031, 2.9, col = 'black')
rect(135932, 2.7, 136111, 2.9, col = 'black')
rect(138505, 2.7, 138647, 2.9, col = 'black')
rect(139391, 2.7, 139498, 2.9, col = 'black')
rect(140310, 2.7, 140515, 2.9, col = 'black')
dev.off()
# repeat for other populations




#  taj D comparisons

UZ1_exgenic_D <- mean(UZ1_D_period[UZ1_D_period$Position < 75000 | UZ1_D_period$Position > 145000, ]$D, na.rm = TRUE  )
UZ1_genic_D <- mean(UZ1_D_period[UZ1_D_period$Position > 75000 & UZ1_D_period$Position < 145000, ]$D, na.rm = TRUE  )
N225_exgenic_D <- mean(N225_D_period[N225_D_period$Position < 75000 | N225_D_period$Position > 145000, ]$D, na.rm = TRUE  )
N225_genic_D <- mean(N225_D_period[N225_D_period$Position > 75000 & N225_D_period$Position < 145000, ]$D, na.rm = TRUE  )
N175_exgenic_D <- mean(N175_D_period[N175_D_period$Position < 75000 | N175_D_period$Position > 145000, ]$D, na.rm = TRUE  )
N175_genic_D <- mean(N175_D_period[N175_D_period$Position > 75000 & N175_D_period$Position < 145000, ]$D, na.rm = TRUE  )
N125_exgenic_D <- mean(N125_D_period[N125_D_period$Position < 75000 | N125_D_period$Position > 145000, ]$D, na.rm = TRUE  )
N125_genic_D <- mean(N125_D_period[N125_D_period$Position > 75000 & N125_D_period$Position < 145000, ]$D, na.rm = TRUE  )
N75_exgenic_D <- mean(N75_D_period[N75_D_period$Position < 75000 | N75_D_period$Position > 145000, ]$D, na.rm = TRUE  )
N75_genic_D <- mean(N75_D_period[N75_D_period$Position > 75000 & N75_D_period$Position < 145000, ]$D, na.rm = TRUE  )
S25_exgenic_D <- mean(S25_D_period[S25_D_period$Position < 75000 | S25_D_period$Position > 145000, ]$D, na.rm = TRUE  )
S25_genic_D <- mean(S25_D_period[S25_D_period$Position > 75000 & S25_D_period$Position < 145000, ]$D, na.rm = TRUE  )
S125_exgenic_D <- mean(S125_D_period[S125_D_period$Position < 75000 | S125_D_period$Position > 145000, ]$D, na.rm = TRUE  )
S125_genic_D <- mean(S125_D_period[S125_D_period$Position > 75000 & S125_D_period$Position < 145000, ]$D, na.rm = TRUE  )
S175_exgenic_D <- mean(S175_D_period[S175_D_period$Position < 75000 | S175_D_period$Position > 145000, ]$D, na.rm = TRUE  )
S175_genic_D <- mean(S175_D_period[S175_D_period$Position > 75000 & S175_D_period$Position < 145000, ]$D, na.rm = TRUE  )
S225_exgenic_D <- mean(S225_D_period[S225_D_period$Position < 75000 | S225_D_period$Position > 145000, ]$D, na.rm = TRUE  )
S225_genic_D <- mean(S225_D_period[S225_D_period$Position > 75000 & S225_D_period$Position < 145000, ]$D, na.rm = TRUE  )
#means and t tests
period_genic_D <- c(N225_genic_D, N175_genic_D, N125_genic_D, N75_genic_D, S25_genic_D, S125_genic_D, S175_genic_D, S225_genic_D)
period_exgenic_D <- c(N225_exgenic_D, N175_exgenic_D, N125_exgenic_D, N75_exgenic_D, S25_exgenic_D, S125_exgenic_D, S175_exgenic_D, S225_exgenic_D)
mean(period_genic_D)
mean(period_exgenic_D)
t.test(period_genic_D, period_exgenic_D)

t.test(period_genic_D)
t.test(period_exgenic_D)


N225_exgenic_D <- mean(N225_D_pdfr[N225_D_pdfr$Position < 44000 | N225_D_pdfr$Position > 83000, ]$D, na.rm = TRUE  )
N225_genic_D <- mean(N225_D_pdfr[N225_D_pdfr$Position > 44000 & N225_D_pdfr$Position < 83000, ]$D, na.rm = TRUE  )
N175_exgenic_D <- mean(N175_D_pdfr[N175_D_pdfr$Position < 44000 | N175_D_pdfr$Position > 83000, ]$D, na.rm = TRUE  )
N175_genic_D <- mean(N175_D_pdfr[N175_D_pdfr$Position > 44000 & N175_D_pdfr$Position < 83000, ]$D, na.rm = TRUE  )
N125_exgenic_D <- mean(N125_D_pdfr[N125_D_pdfr$Position < 44000 | N125_D_pdfr$Position > 83000, ]$D, na.rm = TRUE  )
N125_genic_D <- mean(N125_D_pdfr[N125_D_pdfr$Position > 44000 & N125_D_pdfr$Position < 83000, ]$D, na.rm = TRUE  )
N75_exgenic_D <- mean(N75_D_pdfr[N75_D_pdfr$Position < 44000 | N75_D_pdfr$Position > 83000, ]$D, na.rm = TRUE  )
N75_genic_D <- mean(N75_D_pdfr[N75_D_pdfr$Position > 44000 & N75_D_pdfr$Position < 83000, ]$D, na.rm = TRUE  )
S25_exgenic_D <- mean(S25_D_pdfr[S25_D_pdfr$Position < 44000 | S25_D_pdfr$Position > 83000, ]$D, na.rm = TRUE  )
S25_genic_D <- mean(S25_D_pdfr[S25_D_pdfr$Position > 44000 & S25_D_pdfr$Position < 83000, ]$D, na.rm = TRUE  )
S125_exgenic_D <- mean(S125_D_pdfr[S125_D_pdfr$Position < 44000 | S125_D_pdfr$Position > 83000, ]$D, na.rm = TRUE  )
S125_genic_D <- mean(S125_D_pdfr[S125_D_pdfr$Position > 44000 & S125_D_pdfr$Position < 83000, ]$D, na.rm = TRUE  )
S175_exgenic_D <- mean(S175_D_pdfr[S175_D_pdfr$Position < 44000 | S175_D_pdfr$Position > 83000, ]$D, na.rm = TRUE  )
S175_genic_D <- mean(S175_D_pdfr[S175_D_pdfr$Position > 44000 & S175_D_pdfr$Position < 83000, ]$D, na.rm = TRUE  )
S225_exgenic_D <- mean(S225_D_pdfr[S225_D_pdfr$Position < 44000 | S225_D_pdfr$Position > 83000, ]$D, na.rm = TRUE  )
S225_genic_D <- mean(S225_D_pdfr[S225_D_pdfr$Position > 44000 & S225_D_pdfr$Position < 83000, ]$D, na.rm = TRUE  )

#means and t tests
pdfr_genic_D <- c(N225_genic_D, N175_genic_D, N125_genic_D, N75_genic_D, S25_genic_D, S125_genic_D, S175_genic_D, S225_genic_D)
pdfr_exgenic_D <- c(N225_exgenic_D, N175_exgenic_D, N125_exgenic_D, N75_exgenic_D, S25_exgenic_D, S125_exgenic_D, S175_exgenic_D, S225_exgenic_D)
mean(pdfr_genic_D)
mean(pdfr_exgenic_D)
t.test(pdfr_genic_D, pdfr_exgenic_D)

t.test(pdfr_genic_D)
t.test(pdfr_exgenic_D)





### plotting for pdfr

FST_pdfr <- FST_m[which(FST_m$Scaffold == "scaffold87"), ]
plot(FST_pdfr$Window, FST_pdfr$N225_S75, cex = 0.2, pch = 19, ylim = c(0,1))
pdf("pdfr_fst_N225_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$N225_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$N225_UZ_ma, type = 'l', col = 'red')
rect(44493, 0.9, 46265, 1, col = 'black')
rect(46266, 0.9, 46349, 1, col = 'black')
rect(47065, 0.9, 47184, 1, col = 'black')
rect(52728, 0.9, 52847, 1, col = 'black')
rect(53062, 0.9, 53180, 1, col = 'black')
rect(54524, 0.9, 54707, 1, col = 'black')
rect(58394, 0.9, 58479, 1, col = 'black')
rect(60223, 0.9, 60292, 1, col = 'black')
rect(60605, 0.9, 60834, 1, col = 'black')
rect(61975, 0.9, 62125, 1, col = 'black')
rect(62228, 0.9, 62308, 1, col = 'black')
rect(66266, 0.9, 66411, 1, col = 'black')
rect(79859, 0.9, 80216, 1, col = 'green')
rect(82804, 0.9, 83222, 1, col = 'green')
par(op)
dev.off()

#repeat for other populations

library(dplyr)
#period plotting of FST
FST_period <- FST_m[which(FST_m$Scaffold == "scaffold532"), ]

UZ_period <- FST_period[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_period <- FST_period[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_period <- mutate_all(UZ_period, function(x) as.numeric(as.character(x)))
Other_comparisons_period <- mutate_all(Other_comparisons_period, function(x) as.numeric(as.character(x)))

UZ_period <- cbind(UZ_period, FST_period$Window)
UZ_period$color <- "grey"
UZ_period$color[66:137] <- "red"

UZ_period_ordered <- UZ_period[order(-UZ_period$N75_UZ),]
UZ_period$observation <- 1:nrow(UZ_period)
UZ_period_ordered$observation <- 1:nrow(UZ_period_ordered)
plot(UZ_period_ordered$observation, UZ_period_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_period_ordered$color)
UZ_period_ordered_N125 <- UZ_period[order(-UZ_period$N125_UZ),]
points(1:nrow(UZ_period_ordered_N125), UZ_period_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_period_ordered$color, pch = 19, cex = .5)
UZ_period_ordered_N175 <- UZ_period[order(-UZ_period$N175_UZ),]
points(1:nrow(UZ_period_ordered_N175), UZ_period_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_period_ordered$color, pch = 19, cex = .5)
UZ_period_ordered_N225 <- UZ_period[order(-UZ_period$N225_UZ),]
points(1:nrow(UZ_period_ordered_N225), UZ_period_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_period_ordered$color, pch = 19, cex = .5)
UZ_period_ordered_S25 <- UZ_period[order(-UZ_period$S25_UZ),]
points(1:nrow(UZ_period_ordered_S25), UZ_period_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_period_ordered$color, pch = 19, cex = .5)
UZ_period_ordered_S125 <- UZ_period[order(-UZ_period$S125_UZ),]
points(1:nrow(UZ_period_ordered_S125), UZ_period_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_period_ordered$color, pch = 19, cex = .5)
UZ_period_ordered_S175 <- UZ_period[order(-UZ_period$S175_UZ),]
points(1:nrow(UZ_period_ordered_S175), UZ_period_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_period_ordered$color, pch = 19, cex = .5)
UZ_period_ordered_S225 <- UZ_period[order(-UZ_period$S225_UZ),]
points(1:nrow(UZ_period_ordered_S225), UZ_period_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_period_ordered$color, pch = 19, cex = .5)

#pdfr plotting of FST
## pdfr
FST_pdfr <- FST_m[which(FST_m$Scaffold == "scaffold87"), ]
FST_pdfr <- FST_pdfr[sample(nrow(FST_pdfr), 306), ]
UZ_pdfr <- FST_pdfr[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_pdfr <- FST_pdfr[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_pdfr <- mutate_all(UZ_pdfr, function(x) as.numeric(as.character(x)))
Other_comparisons_pdfr <- mutate_all(Other_comparisons_pdfr, function(x) as.numeric(as.character(x)))

UZ_pdfr <- cbind(UZ_pdfr, FST_pdfr$Window)
UZ_pdfr$color <- "lightblue"
UZ_pdfr$color[26:57] <- "red"

UZ_pdfr_ordered <- UZ_pdfr[order(-UZ_pdfr$N75_UZ),]
UZ_pdfr$observation <- 1:nrow(UZ_pdfr)
UZ_pdfr_ordered$observation <- 1:nrow(UZ_pdfr_ordered)
plot(UZ_pdfr_ordered$observation, UZ_pdfr_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_pdfr_ordered$color)
UZ_pdfr_ordered_N125 <- UZ_pdfr[order(-UZ_pdfr$N125_UZ),]
points(1:nrow(UZ_pdfr_ordered_N125), UZ_pdfr_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_pdfr_ordered$color, pch = 19, cex = .5)
UZ_pdfr_ordered_N175 <- UZ_pdfr[order(-UZ_pdfr$N175_UZ),]
points(1:nrow(UZ_pdfr_ordered_N175), UZ_pdfr_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_pdfr_ordered$color, pch = 19, cex = .5)
UZ_pdfr_ordered_N225 <- UZ_pdfr[order(-UZ_pdfr$N225_UZ),]
points(1:nrow(UZ_pdfr_ordered_N225), UZ_pdfr_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_pdfr_ordered$color, pch = 19, cex = .5)
UZ_pdfr_ordered_S25 <- UZ_pdfr[order(-UZ_pdfr$S25_UZ),]
points(1:nrow(UZ_pdfr_ordered_S25), UZ_pdfr_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_pdfr_ordered$color, pch = 19, cex = .5)
UZ_pdfr_ordered_S125 <- UZ_pdfr[order(-UZ_pdfr$S125_UZ),]
points(1:nrow(UZ_pdfr_ordered_S125), UZ_pdfr_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_pdfr_ordered$color, pch = 19, cex = .5)
UZ_pdfr_ordered_S175 <- UZ_pdfr[order(-UZ_pdfr$S175_UZ),]
points(1:nrow(UZ_pdfr_ordered_S175), UZ_pdfr_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_pdfr_ordered$color, pch = 19, cex = .5)
UZ_pdfr_ordered_S225 <- UZ_pdfr[order(-UZ_pdfr$S225_UZ),]
points(1:nrow(UZ_pdfr_ordered_S225), UZ_pdfr_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_pdfr_ordered$color, pch = 19, cex = .5)




