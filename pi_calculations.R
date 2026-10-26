# script to process popoolation stats
# showing population N225 here, but repeat for other populations
setwd("~/Desktop/pi")
library(tidyr)
library(adegenet)
library(ape)
library(vegan)
library(RColorBrewer)
z_chr_scaffolds_order <- c('scaffold4435', 'scaffold751', 'scaffold35', 'scaffold572', 'scaffold740', 'scaffold1546', 'scaffold247', 'scaffold196', 'scaffold779', 'scaffold2358', 'scaffold600', 'scaffold216', 'scaffold390', 'scaffold364', 'scaffold1895', 'scaffold409', 'scaffold444', 'scaffold325', 'scaffold15', 'scaffold169', 'scaffold348', 'scaffold865', 'scaffold518', 'scaffold2225', 'scaffold546', 'scaffold645', 'scaffold1618', 'scaffold789', 'scaffold218', 'scaffold1660', 'scaffold1044', 'scaffold581', 'scaffold1237', 'scaffold1246', 'scaffold1756', 'scaffold1161', 'scaffold155', 'scaffold282', 'scaffold804', 'scaffold502', 'scaffold3599', 'scaffold1040', 'scaffold275', 'scaffold4243', 'scaffold532', 'scaffold459', 'scaffold1889', 'scaffold221', 'scaffold931', 'scaffold140', 'scaffold1038', 'scaffold598', 'scaffold48', 'scaffold1486', 'scaffold193', 'scaffold137', 'bab_gene_region_CC', 'scaffold87', 'scaffold50', 'scaffold662', 'scaffold728', 'scaffold229', 'scaffold330', 'scaffold673', 'scaffold812', 'scaffold1046')
#read input data
N225 <- read.csv("N225_pi1filtered.kb.txt", sep = '\t', col.names = c("Scaffold", "Window", "Depth", "X", "pi"))

N225$pop = "N225"

N225$ID <- paste(N225$Scaffold, N225$Window, sep = '_')
# filtering z and autosomes
N225_z <- N225[N225$Scaffold %in% z_chr_scaffolds_order, ]
N225_a<- N225[!N225$Scaffold %in% z_chr_scaffolds_order, ]

N225_a_r <- N225_a[, c(6, 9)]
#reordering columns
pi_all <- merge(N225_a_r, N175_a_r, by = 'ID')
pi_all_2 <- merge(pi_all, N125_a_r, by = 'ID')
pi_all_3 <- merge(pi_all_2, N75_a_r, by = 'ID')
pi_all_4 <- merge(pi_all_3, S25_a_r, by = 'ID')
pi_all_5 <- merge(pi_all_4, S125_a_r, by = 'ID')
pi_all_6 <- merge(pi_all_5, S175_a_r, by = 'ID')
pi_all_7 <- merge(pi_all_6, S225_a_r, by = 'ID')
pi_all_complete <- merge(pi_all_7, UZ_a_r, by = 'ID')
colnames(pi_all_complete) <- c("ID", "N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ")
#complete N225, repeat for other populations
pi_all_complete$N225 <- as.numeric(pi_all_complete$N225)

#bind data to pi complete frame
data_m <- pi_all_complete[, 2]
data_m <- as.data.frame(data_m)
data_m <- rbind(data_m, pi_all_complete[, 3])
data_m <- rbind(data_m, pi_all_complete[, 4])
data_m <- rbind(data_m, pi_all_complete[, 5])
data_m <- rbind(data_m, pi_all_complete[, 6])
data_m <- rbind(data_m, pi_all_complete[, 7])
data_m <- rbind(data_m, pi_all_complete[, 8])
data_m <- rbind(data_m, pi_all_complete[, 9])
#quick mean and ttests
mean(data_m[, 1], na.rm= TRUE)
t.test(data_m[, 1], na.rm=TRUE)
mean(pi_all_complete[, 10], na.rm= TRUE)

#complete cases only
pi_all_complete <- pi_all_complete[complete.cases(pi_all_complete), ]
data_m <- pi_all_complete[ ]
#pcoa 
dds.pcoa=pcoa(vegdist(t((pi_all_complete[,2:10])), method = "euclidean")/1000)

scores=dds.pcoa$vectors

percent <- dds.pcoa$values$Eigenvalues

percent /sum(percent)

#plotting pcoa
#palette <- brewer.pal(9, "PuOr")
png("~/Desktop/autosomal_pi_pcoa_2.png", height = 400, width = 400)
pdf("~/Desktop/autosomal_pi_pcoa_2.pdf", height = 6, width = 6)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=FALSE)
plot(scores[,1], scores[,2], col = palette, pch = 19, xlab = "PC1 (68.1%)", ylab = "PC2 (6.6%)", cex = 2, yaxt = 'n', xaxt = 'n')
abline(h = 0, lty = 'dashed')
abline(v = 0, lty = 'dashed')
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright", bg = 'white', c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ"), col = palette, pch = 19, inset=c(-.2,-.2))
dev.off()

#### Z ####

N225_z_r <- N225_z[, c(6, 9)]

pi_all <- merge(N225_z_r, N175_z_r, by = 'ID')
pi_all_2 <- merge(pi_all, N125_z_r, by = 'ID')
pi_all_3 <- merge(pi_all_2, N75_z_r, by = 'ID')
pi_all_4 <- merge(pi_all_3, S25_z_r, by = 'ID')
pi_all_5 <- merge(pi_all_4, S125_z_r, by = 'ID')
pi_all_6 <- merge(pi_all_5, S175_z_r, by = 'ID')
pi_all_7 <- merge(pi_all_6, S225_z_r, by = 'ID')
pi_all_complete <- merge(pi_all_7, UZ_z_r, by = 'ID')
colnames(pi_all_complete) <- c("ID", "N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ")
# repeat for other populations
pi_all_complete$N225 <- as.numeric(pi_all_complete$N225)


# bind relevant columns
data_m <- pi_all_complete[, 2]
data_m <- as.data.frame(data_m)
data_m <- rbind(data_m, pi_all_complete[, 3])
data_m <- rbind(data_m, pi_all_complete[, 4])
data_m <- rbind(data_m, pi_all_complete[, 5])
data_m <- rbind(data_m, pi_all_complete[, 6])
data_m <- rbind(data_m, pi_all_complete[, 7])
data_m <- rbind(data_m, pi_all_complete[, 8])
data_m <- rbind(data_m, pi_all_complete[, 9])
#quick means and ttests
mean(data_m[, 1], na.rm= TRUE)
t.test(data_m[, 1], na.rm= TRUE)
mean(pi_all_complete[, 10], na.rm= TRUE)

#pcoas on complete data (missing removed)
pi_all_complete <- pi_all_complete[complete.cases(pi_all_complete), ]

dds.pcoa=pcoa(vegdist(t((pi_all_complete[,2:10])), method = "euclidean")/1000)

scores=dds.pcoa$vectors

percent <- dds.pcoa$values$Eigenvalues

percent /sum(percent)
#plot pcoa 

palette <- c(brewer.pal(8, "PRGn"), "black")
png("~/Desktop/zlinked_pi_pcoa_2.png", height = 400, width = 400)
pdf("~/Desktop/zlinked_pi_pcoa_2.pdf", height = 6, width = 6)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=FALSE)
plot(scores[,1], -scores[,2], col = palette, pch = 19, xlab = "PC1 (67.6%)", ylab = "PC2 (6.2%)", cex = 2, xaxt = "n", yaxt = "n")
abline(h = 0, lty = 'dashed')
abline(v = 0, lty = 'dashed')
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#legend("topright", c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ"), col = palette, pch = 19, inset=c(-.2,-.2))
dev.off()


#### end Z ####
mean(mean(N225_z_r$pi, na.rm = TRUE),

mean(mean(N225_a_r$pi, na.rm = TRUE),
# ttests of multiple pops, make sure you've applied everything above on all pops you want to test
t.test(c(N225_a_r$pi, (N175_a_r$pi), (N125_a_r$pi), (N75_a_r$pi),
(S25_a_r$pi),
(S125_a_r$pi),
(S175_a_r$pi),
(S225_a_r$pi)))

#UZ population handled separately
mean(UZ_a_r$pi, na.rm = TRUE)
t.test(UZ_a_r$pi, na.rm = TRUE)
t.test(UZ_z_r$pi, na.rm = TRUE)


test_matrix <- cbind(N225_a$pi)
#plotting pi, repeat for other pops
plot(N225$observation, N225$pi, cex = 0.02, type = 'l', ylim = c(0.0, 0.5))
mean(N225$pi, na.rm = TRUE)

#filtering down to individual scaffolds
N225_532 <- N225[which(N225$Scaffold == 'scaffold532'), ]

ma <- function(x, n = 10){stats::filter(x, rep(1 / n, n), sides = 2)}
#setting 0s to NA
N225_532$pi[is.na(N225_532$pi)] <- 0
#moving averages
N225_532$pi_ma <- ma(N225_532$pi)
#period plot, repeat for other populations and scaffolds
pdf("period_pi_N225.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(N225_532$Window, N225_532$pi, pch = 19, cex = 0.25, ylim = c(0, 0.03))
points(N225_532$Window, N225_532$pi_ma, type = 'l', col = 'red')
rect(75057, 0.026, 75082, 0.029 , col = 'green')
rect(77479, 16, 77748, 17 , col = 'green')
rect(101450, 0.026, 101463, 0.029 , col = 'green')
rect(101464, 0.026, 101551, 0.029 , col = 'black')
rect(103643, 0.026, 103728, 0.029 , col = 'black')
rect(104601, 0.026, 104754, 0.029 , col = 'black')
rect(110704, 0.026, 110801, 0.029 , col = 'black')
rect(110887, 0.026, 110970, 0.029 , col = 'black')
rect(112874, 0.026, 113006, 0.029 , col = 'black')
rect(113155, 0.026, 113214, 0.029 , col = 'black')
rect(113215, 0.026, 113231, 0.029 , col = 'black')
rect(114358, 0.026, 114563, 0.029 , col = 'black')
rect(114924, 0.026, 115155, 0.029 , col = 'black')
rect(115522, 0.026, 115636, 0.029 , col = 'black')
rect(115734, 0.026, 115858, 0.029 , col = 'black')
rect(116369, 0.026, 116521, 0.029 , col = 'black')
rect(125722, 0.026, 125837, 0.029 , col = 'black')
rect(125922, 0.026, 126020, 0.029 , col = 'black')
rect(128758, 0.026, 128908, 0.029 , col = 'black')
rect(130219, 0.026, 130411, 0.029 , col = 'black')
rect(132239, 0.026, 132332, 0.029 , col = 'black')
rect(133776, 0.026, 133915, 0.029 , col = 'black')
rect(134917, 0.026, 135031, 0.029 , col = 'black')
rect(135932, 0.026, 136111, 0.029 , col = 'black')
rect(138505, 0.026, 138647, 0.029 , col = 'black')
rect(139391, 0.026, 139498, 0.029 , col = 'black')
rect(140310, 0.026, 140515, 0.029 , col = 'black')
par(op)
dev.off()

## pdfr repeated the above steps but swapping the scaffold

