# script to process popoolation stats
# filtered for an individual population here, but repeated for each population
setwd("~/Desktop/D")
library(tidyr)
library(dplyr)
library(ggplot2)
#open D files per pop and rename columns
N225 <- read.csv("N225_D1filtered.kb.txt", sep = '\t', col.names = c("Scaffold", "Window", "nSNPs", "X", "D"))

# filter
N225_complete <- N225[complete.cases(N225), ]
N225_complete <- separate(data = N225_complete, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
N225_complete$Scaffold<-gsub("v02","",as.character(N225_complete$Scaffold))
N225_complete$length<-gsub("v02","",as.character(N225_complete$length))
N225_complete$length<-replace_na(N225_complete$length,"414000")
N225_complete$length <- as.numeric(N225_complete$length)

z_chr_scaffolds_order <- c('scaffold4435', 'scaffold751', 'scaffold35', 'scaffold572', 'scaffold740', 'scaffold1546', 'scaffold247', 'scaffold196', 'scaffold779', 'scaffold2358', 'scaffold600', 'scaffold216', 'scaffold390', 'scaffold364', 'scaffold1895', 'scaffold409', 'scaffold444', 'scaffold325', 'scaffold15', 'scaffold169', 'scaffold348', 'scaffold865', 'scaffold518', 'scaffold2225', 'scaffold546', 'scaffold645', 'scaffold1618', 'scaffold789', 'scaffold218', 'scaffold1660', 'scaffold1044', 'scaffold581', 'scaffold1237', 'scaffold1246', 'scaffold1756', 'scaffold1161', 'scaffold155', 'scaffold282', 'scaffold804', 'scaffold502', 'scaffold3599', 'scaffold1040', 'scaffold275', 'scaffold4243', 'scaffold532', 'scaffold459', 'scaffold1889', 'scaffold221', 'scaffold931', 'scaffold140', 'scaffold1038', 'scaffold598', 'scaffold48', 'scaffold1486', 'scaffold193', 'scaffold137', 'scaffold87', 'scaffold50', 'scaffold662', 'scaffold728', 'scaffold229', 'scaffold330', 'scaffold673', 'scaffold812', 'scaffold1046')
#filter z vs autosomes
N225_D_z <- N225_complete[N225_complete$Scaffold %in% z_chr_scaffolds_order, ]
N225_D_a <- N225_complete[!N225_complete$Scaffold %in% z_chr_scaffolds_order, ]
#quick mean and t test check
mean(N225_complete$D)
t.test(N225_complete$D)

median(N225_complete$D)


N225_D_z <- N225_complete[N225_complete$Scaffold %in% z_chr_scaffolds_order, ]
N225_D_a <- N225_complete[!N225_complete$Scaffold %in% z_chr_scaffolds_order, ]
mean(N225_complete$D)
median(N225_complete$D)



mean(N225_D_z$D)
median(N225_D_z$D)


mean(N225_D_a$D)
median(N225_D_a$D)
#filter to pdfr scaffold
N225_D_z_pdfr <- N225_D_z[N225_D_z$Scaffold == "scaffold87", ]
#filter to period
N225_D_z_period <- N225_D_z[N225_D_z$Scaffold == "scaffold532", ]
# genic regions
N225_D_z_period_genic <- N225_D_z_period[(N225_D_z_period$Window < 149301 & N225_D_z_period$Window > 75057) , ]
#exgenic region
N225_D_z_period_exgenic <- N225_D_z_period[!(N225_D_z_period$Window < 149301 & N225_D_z_period$Window > 75057) , ]


N225_D_z_pdfr_exgenic <- N225_D_z_pdfr[!(N225_D_z_pdfr$Window < 83222 & N225_D_z_pdfr$Window > 44493) , ]

N225_D_z_pdfr_genic <- N225_D_z_pdfr[(N225_D_z_pdfr$Window < 83222 & N225_D_z_pdfr$Window > 44493) , ]
#wilcox tests
wilcox.test(N225_D_z_period_genic$D, N225_D_z_period_exgenic$D)
#t tests, relies on having applied the above code to the other populations
t.test(c(N225_D_z_period_genic$D, N175_D_z_period_genic$D, N125_D_z_period_genic$D, N75_D_z_period_genic$D, S25_D_z_period_genic$D, S125_D_z_period_genic$D, S175_D_z_period_genic$D, S225_D_z_period_genic$D))
t.test(c(N225_D_z_pdfr_genic$D, N175_D_z_pdfr_genic$D, N125_D_z_pdfr_genic$D, N75_D_z_pdfr_genic$D, S25_D_z_pdfr_genic$D, S125_D_z_pdfr_genic$D, S175_D_z_pdfr_genic$D, S225_D_z_pdfr_genic$D))

t.test(c(N225_D_z_period_genic$D, N175_D_z_period_genic$D, N125_D_z_period_genic$D, N75_D_z_period_genic$D, S25_D_z_period_genic$D, S125_D_z_period_genic$D, S175_D_z_period_genic$D, S225_D_z_period_genic$D, UZ_D_z_period_genic$D))
t.test(c(N225_D_z_pdfr_genic$D, N175_D_z_pdfr_genic$D, N125_D_z_pdfr_genic$D, N75_D_z_pdfr_genic$D, S25_D_z_pdfr_genic$D, S125_D_z_pdfr_genic$D, S175_D_z_pdfr_genic$D, S225_D_z_pdfr_genic$D, UZ_D_z_pdfr_genic$D))


wilcox.test(N225_D_z_pdfr_genic$D, N225_D_z_pdfr_exgenic$D)

N225_D_a$cat <- 1
N225_D_z$cat <- 2
N225_D <- rbind(N225_D_a, N225_D_z)
N225_D$pop <- 1

#combining and plotting statistics
comb <- rbind(N225_D, N175_D, N125_D, N75_D, S25_D, S125_D, S175_D, S225_D, UZ_D)
png("~/Desktop/td_boxplot.png", height = 400, width = 800)
pdf("~/Desktop/td_boxplot.pdf", height = 6, width = 12)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=FALSE)
boxplot(D~cat*pop, data = comb, ylab = "Tajima's D", names= c("N225", "", "N175", "", "N125", "", "N75", "", "S25", "", "S125", "", "S175", "", "S225", "", "UZ", ""), col = c("darkblue", "lightblue"), outline = FALSE, las = 2, notch = FALSE,
        at = c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23, 27, 28), xaxt = "n")
axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5, 22.5, 27.5), labels = c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ"))
abline(h = 0, lty = 'dashed')
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
legend("topright", bg = 'white',  inset=c(-.1,-.1), legend = c("Autosomal", "Z-linked"), col = c("darkblue", "lightblue"), pch = 15)
dev.off()
boxplot(N225_D_a$D)
boxplot(N225_D_z$D)

mean(comb$D)

comb_cline <- comb[!comb$pop == 9, ]
comb_cline_z <- comb_cline[comb_cline$cat == 2, ]
comb_cline_auto <- comb_cline[comb_cline$cat == 1, ]

mean(comb_cline_auto$D)
mean(comb_cline_z$D)

comb_uz <- comb[comb$pop == 9, ]
comb_uz_z <- comb_uz[comb_uz$cat == 2, ]
comb_uz_auto <- comb_uz[comb_uz$cat == 1, ]
mean(comb_uz_auto$D)
mean(comb_uz_z$D)




N225_complete <- N225[complete.cases(N225), ]
N225_complete <- separate(data = N225_complete, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
N225_complete$Scaffold<-gsub("v02","",as.character(N225_complete$Scaffold))
N225_complete$length<-gsub("v02","",as.character(N225_complete$length))
N225_complete$length<-replace_na(N225_complete$length,"414000")
N225_complete$length <- as.numeric(N225_complete$length)
N225_complete <- N225_complete[N225_complete$length >= 50000, ]

mean(N225_complete$D)

N225_complete_sorted <- N225_complete[order(-N225_complete$D), ]
N225_complete_sorted_r <- N225_complete[order(N225_complete$D), ]

z_linked <- c('scaffold4435', 'scaffold751', 'scaffold35', 'scaffold572', 'scaffold740', 'scaffold1546', 'scaffold247', 'scaffold196', 'scaffold779', 'scaffold2358', 'scaffold600', 'scaffold216', 'scaffold390', 'scaffold364', 'scaffold1895', 'scaffold409', 'scaffold444', 'scaffold325', 'scaffold15', 'scaffold169', 'scaffold348', 'scaffold865', 'scaffold518', 'scaffold2225', 'scaffold546', 'scaffold645', 'scaffold1618', 'scaffold789', 'scaffold218', 'scaffold1660', 'scaffold1044', 'scaffold581', 'scaffold1237', 'scaffold1246', 'scaffold1756', 'scaffold1161', 'scaffold155', 'scaffold282', 'scaffold804', 'scaffold502', 'scaffold3599', 'scaffold1040', 'scaffold275', 'scaffold4243', 'scaffold532', 'scaffold459', 'scaffold1889', 'scaffold221', 'scaffold931', 'scaffold140', 'scaffold1038', 'scaffold598', 'scaffold48', 'scaffold1486', 'scaffold193', 'scaffold137', 'scaffold87', 'scaffold50', 'scaffold662', 'scaffold728', 'scaffold229', 'scaffold330', 'scaffold673', 'scaffold812', 'scaffold1046')
N225_z <- N225_complete[N225_complete$Scaffold %in% z_linked, ]
N225_z_ordered <- N225_z %>% arrange(factor(Scaffold, levels = z_linked))
N225_z_ordered$observation <- 1:nrow(N225_z_ordered)
plot(N225_z_ordered$observation, N225_z_ordered$D, pch = 19, cex = 0.3)
ma <- function(x, n = 100){stats::filter(x, rep(1 / n, n), sides = 2)}
N225_z_ordered$D <- as.numeric(N225_z_ordered$D)
ma_N225_z <- ma(N225_z_ordered$D)
points(1:nrow(N225_z_ordered), ma_N225_z, type = 'l', col = 'red')
N225_z_ordered$moving_average <- ma_N225_z
plot(N225_z_ordered$observation, N225_z_ordered$D, pch = 19, cex = 0.3, xlim = c(12500, 13500))
points(1:nrow(N225_z_ordered), ma_N225_z, type = 'l', col = 'red')
plot(N225_z_ordered$observation, N225_z_ordered$D, pch = 19, cex = 0.3, xlim = c(9800, 10300))
points(1:nrow(N225_z_ordered), ma_N225_z, type = 'l', col = 'red')

mean(N225_z$D)
N225_z_ordered_by_ma <- N225_z_ordered[order(N225_z_ordered$moving_average), ]
N225_z_ordered_by_ma_f <- N225_z_ordered[order(-N225_z_ordered$moving_average), ]

hist(N225_z$D, breaks = 50)
hist(N225_z_ordered$moving_average, breaks = 50)

N225_140 <- N225_complete[N225_complete$Scaffold == 'scaffold140', ]
plot(N225_140$Window, N225_140$D, pch = 19, cex = 0.5, xlim = c(380000, 540000))


N225_216 <- N225_complete[N225_complete$Scaffold == 'scaffold216', ]
plot(N225_216$Window, N225_216$D, pch = 19)
N225_87 <- N225_complete[N225_complete$Scaffold == 'scaffold87', ]
plot(N225_87$Window, N225_87$D, pch = 19, cex = 0.5)
N225_532 <- N225_complete[N225_complete$Scaffold == 'scaffold532', ]
plot(N225_532$Window, N225_532$D, pch = 19, cex = 0.5)

N225_153 <- N225_complete[N225_complete$Scaffold == 'scaffold153', ]
plot(N225_153$Window, N225_153$D, pch = 19)
N225_152 <- N225_complete[N225_complete$Scaffold == 'scaffold152', ]
plot(N225_152$Window, N225_152$D, pch = 19)
N225_155 <- N225_complete[N225_complete$Scaffold == 'scaffold155', ]
plot(N225_155$Window, N225_155$D, pch = 19)

N225_153 <- N225_complete[N225_complete$Scaffold == 'scaffold153', ]
plot(N225_153$Window, N225_153$D, pch = 19)

N225_182 <- N225_complete[N225_complete$Scaffold == 'scaffold182', ]
plot(N225_182$Window, N225_182$D, pch = 19)
N225_178 <- N225_complete[N225_complete$Scaffold == 'scaffold178', ]
plot(N225_178$Window, N225_178$D, pch = 19)

# plotting D across period for N225 pop
pdf("period_D_N225.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(N225_532$Window, N225_532$D, pch = 19, cex = 0.25, ylim = c(-3, 3))
abline(h = 0, lty = 'dashed')
points(N225_532$Window, N225_532$D_ma, type = 'l', col = 'red')
rect(75057, 2.5, 75082, 3.0 , col = 'green')
rect(77479, 16, 77748, 17 , col = 'green')
rect(101450, 2.5, 101463, 3.0 , col = 'green')
rect(101464, 2.5, 101551, 3.0 , col = 'black')
rect(103643, 2.5, 103728, 3.0 , col = 'black')
rect(104601, 2.5, 104754, 3.0 , col = 'black')
rect(110704, 2.5, 110801, 3.0 , col = 'black')
rect(110887, 2.5, 110970, 3.0 , col = 'black')
rect(112874, 2.5, 113006, 3.0 , col = 'black')
rect(113155, 2.5, 113214, 3.0 , col = 'black')
rect(113215, 2.5, 113231, 3.0 , col = 'black')
rect(114358, 2.5, 114563, 3.0 , col = 'black')
rect(114924, 2.5, 115155, 3.0 , col = 'black')
rect(115522, 2.5, 115636, 3.0 , col = 'black')
rect(115734, 2.5, 115858, 3.0 , col = 'black')
rect(116369, 2.5, 116521, 3.0 , col = 'black')
rect(125722, 2.5, 125837, 3.0 , col = 'black')
rect(125922, 2.5, 126020, 3.0 , col = 'black')
rect(128758, 2.5, 128908, 3.0 , col = 'black')
rect(130219, 2.5, 130411, 3.0 , col = 'black')
rect(132239, 2.5, 132332, 3.0 , col = 'black')
rect(133776, 2.5, 133915, 3.0 , col = 'black')
rect(134917, 2.5, 135031, 3.0 , col = 'black')
rect(135932, 2.5, 136111, 3.0 , col = 'black')
rect(138505, 2.5, 138647, 3.0 , col = 'black')
rect(139391, 2.5, 139498, 3.0 , col = 'black')
rect(140310, 2.5, 140515, 3.0 , col = 'black')
par(op)
dev.off()


#plotting D for pdfr for pop N225


pdf("pdfr_D_N225.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(N225_87$Window, N225_87$D, pch = 19, cex = 0.25, ylim = c(-3, 3))
rect(46266, 2.5, 46349, 3.0, col = 'black')
rect(47065, 2.5, 47184, 3.0, col = 'black')
rect(52728, 2.5, 52847, 3.0, col = 'black')
rect(53062, 2.5, 53180, 3.0, col = 'black')
rect(54524, 2.5, 54707, 3.0, col = 'black')
rect(58394, 2.5, 58479, 3.0, col = 'black')
rect(60223, 2.5, 60292, 3.0, col = 'black')
rect(60605, 2.5, 60834, 3.0, col = 'black')
rect(61975, 2.5, 62125, 3.0, col = 'black')
rect(62228, 2.5, 62308, 3.0, col = 'black')
rect(66266, 2.5, 66411, 3.0, col = 'black')
rect(79859, 2.5, 80216, 3.0, col = 'green')
rect(82804, 2.5, 83222, 3.0, col = 'green')
abline(h = 0, lty = 'dashed')
points(N225_87$Window, N225_87$D_ma, type = 'l', col = 'red')
par(op)
dev.off()


plot(N225_532$Window[180:250], N225_532$D[180:250], pch = 19, cex = 0.5)


plot(N225_87$Window, N225_87$D, pch = 19, cex = 0.5)

N225_153 <- N225[which(N225$Scaffold == 'scaffold153|size1284462'), ]
plot(N225_153$Window, N225_153$D)







#N225
N225_complete <- N225[complete.cases(N225), ]
N225_complete <- separate(data = N225_complete, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
N225_complete$Scaffold<-gsub("v02","",as.character(N225_complete$Scaffold))
N225_complete$length<-gsub("v02","",as.character(N225_complete$length))
N225_complete$length<-replace_na(N225_complete$length,"414000")
N225_complete$length <- as.numeric(N225_complete$length)
N225_complete <- N225_complete[N225_complete$length >= 50000, ]

mean(N225_complete$D)

N225_complete_sorted <- N225_complete[order(-N225_complete$D), ]
N225_complete_sorted_r <- N225_complete[order(N225_complete$D), ]

z_linked <- c('scaffold4435', 'scaffold751', 'scaffold35', 'scaffold572', 'scaffold740', 'scaffold1546', 'scaffold247', 'scaffold196', 'scaffold779', 'scaffold2358', 'scaffold600', 'scaffold216', 'scaffold390', 'scaffold364', 'scaffold1895', 'scaffold409', 'scaffold444', 'scaffold325', 'scaffold15', 'scaffold169', 'scaffold348', 'scaffold865', 'scaffold518', 'scaffold2225', 'scaffold546', 'scaffold645', 'scaffold1618', 'scaffold789', 'scaffold218', 'scaffold1660', 'scaffold1044', 'scaffold581', 'scaffold1237', 'scaffold1246', 'scaffold1756', 'scaffold1161', 'scaffold155', 'scaffold282', 'scaffold804', 'scaffold502', 'scaffold3599', 'scaffold1040', 'scaffold275', 'scaffold4243', 'scaffold532', 'scaffold459', 'scaffold1889', 'scaffold221', 'scaffold931', 'scaffold140', 'scaffold1038', 'scaffold598', 'scaffold48', 'scaffold1486', 'scaffold193', 'scaffold137', 'scaffold87', 'scaffold50', 'scaffold662', 'scaffold728', 'scaffold229', 'scaffold330', 'scaffold673', 'scaffold812', 'scaffold1046')
N225_z <- N225_complete[N225_complete$Scaffold %in% z_linked, ]
N225_z_ordered <- N225_z %>% arrange(factor(Scaffold, levels = z_linked))
N225_z_ordered$observation <- 1:nrow(N225_z_ordered)
png("~/Desktop/N225_z_plot.png", height = 4, width = 6, units = "in", res = 1200)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(N225_z_ordered$observation, N225_z_ordered$D, pch = 19, cex = 0.3, ylim = c(-3, 3), col = alpha('black', 0.2))
ma <- function(x, n = 100){stats::filter(x, rep(1 / n, n), sides = 2)}
N225_z_ordered$D <- as.numeric(N225_z_ordered$D)
#moving average
ma_N225_z <- ma(N225_z_ordered$D)
points(1:nrow(N225_z_ordered), ma_N225_z, type = 'l', col = 'red')
dev.off()
png("~/Desktop/N225_z_plot_zoomed.png", height = 4, width = 6, units = "in", res = 1200)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(N225_z_ordered$observation, N225_z_ordered$D, pch = 19, cex = 0.3, ylim = c(-3, 3), col = alpha('black', 0.2), xlim = c(12000, 15000))
points(1:nrow(N225_z_ordered), ma_N225_z, type = 'l', col = 'red')
dev.off()
png("~/Desktop/N225_z_plot_zoomed_2.png", height = 4, width = 6, units = "in", res = 1200)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(N225_z_ordered$observation, N225_z_ordered$D, pch = 19, cex = 0.3, ylim = c(-3, 3), col = alpha('black', 0.6), xlim = c(13350, 13700))
points(1:nrow(N225_z_ordered), ma_N225_z, type = 'l', col = 'red')
dev.off()

#plotting moving averages
N225_z_ordered$moving_average <- ma_N225_z
plot(N225_z_ordered$observation, N225_z_ordered$D, pch = 19, cex = 0.3, xlim = c(12500, 13500))
points(1:nrow(N225_z_ordered), ma_N225_z, type = 'l', col = 'red')
plot(N225_z_ordered$observation, N225_z_ordered$D, pch = 19, cex = 0.3, xlim = c(9800, 10300))
points(1:nrow(N225_z_ordered), ma_N225_z, type = 'l', col = 'red')

mean(N225_z$D)
N225_z_ordered_by_ma <- N225_z_ordered[order(N225_z_ordered$moving_average), ]
N225_z_ordered_by_ma_f <- N225_z_ordered[order(-N225_z_ordered$moving_average), ]

hist(N225_z$D, breaks = 50)
hist(N225_z_ordered$moving_average, breaks = 50)
