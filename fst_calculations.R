# script to process popool stats
setwd("~/Desktop/FST")
library(tidyr)
library(vegan)
library(fossil)
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

A <- function(x) gsub(".*=","",x)
FST_m <- cbind(FST_main[1:5], apply(FST_main[6:60],2, A))
B <- function(x) as.numeric(as.character(x))
FST_m <- cbind(FST_main[1:5], apply(FST_main[6:60],2, A))
write.table(FST_m, "FST_by_pop", sep = '\t', quote = FALSE)


# start here for reanalysis
FST_m <- read.table("FST_by_pop", sep = '\t', header = TRUE)
mean(FST_m[, 16])
t.test(FST_m[, 16])

FST_m_z <- FST_m[FST_m$Scaffold %in% z_chr_scaffolds, ]
FST_m_auto <- FST_m[!FST_m$Scaffold %in% z_chr_scaffolds, ]
fst_out_auto <- matrix(nrow = 9, ncol = 9)
fst_out_z <- matrix(nrow = 9, ncol = 9)

# skip if not filtering
FST_m_z <- FST_m_z[FST_m_z$nSNPs > 20, ]
FST_m_z <- FST_m_z[FST_m_z$fraction_covered > 0.20, ]
FST_m_auto <- FST_m_auto[FST_m_auto$nSNPs > 20, ]
FST_m_auto <- FST_m_auto[FST_m_auto$fraction_covered > 0.20, ]

colnames(fst_out_auto) <- c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ")
row.names(fst_out_auto) <- c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ")

write.table(fst_out_auto, "~/Desktop/fst_out_autosomes.csv", quote = FALSE, sep = '\t')
write.table(fst_out_z, "~/Desktop/fst_out_z_linked.csv", quote = FALSE, sep = '\t')

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

fst_s_auto <- cbind((FST_m_auto[, 33]), (FST_m_auto[, 26]), (FST_m_auto[, 18]), (FST_m_auto[, 40]), (FST_m_auto[, 42]), (FST_m_auto[, 43]), (FST_m_auto[, 44]), (FST_m_auto[, 45]), (FST_m_auto[, 25]), (FST_m_auto[, 17]), (FST_m_auto[, 34]), (FST_m_auto[, 36]), (FST_m_auto[, 37]), (FST_m_auto[, 38]), (FST_m_auto[, 39]), (FST_m_auto[, 16]), (FST_m_auto[, 27]), (FST_m_auto[, 29]), (FST_m_auto[, 30]), (FST_m_auto[, 31]), (FST_m_auto[, 32]), (FST_m_auto[, 19]), (FST_m_auto[, 21]), (FST_m_auto[, 22]), (FST_m_auto[, 23]), (FST_m_auto[, 24]), (FST_m_auto[, 47]), (FST_m_auto[, 48]), (FST_m_auto[, 49]), (FST_m_auto[, 50]), (FST_m_auto[, 55]), (FST_m_auto[, 56]), (FST_m_auto[, 57]), (FST_m_auto[, 58]), (FST_m_auto[, 59]), (FST_m_auto[, 60]))
fst_noUZ_auto <- cbind(FST_m_auto[, 1:5], (FST_m_auto[, 33]), (FST_m_auto[, 26]), (FST_m_auto[, 18]), (FST_m_auto[, 40]), (FST_m_auto[, 42]), (FST_m_auto[, 43]), (FST_m_auto[, 44]), (FST_m_auto[, 25]), (FST_m_auto[, 17]), (FST_m_auto[, 34]), (FST_m_auto[, 36]), (FST_m_auto[, 37]), (FST_m_auto[, 38]), (FST_m_auto[, 16]), (FST_m_auto[, 27]), (FST_m_auto[, 29]), (FST_m_auto[, 30]), (FST_m_auto[, 31]), (FST_m_auto[, 19]), (FST_m_auto[, 21]), (FST_m_auto[, 22]), (FST_m_auto[, 23]), (FST_m_auto[, 47]), (FST_m_auto[, 48]), (FST_m_auto[, 49]), (FST_m_auto[, 55]), (FST_m_auto[, 56]), (FST_m_auto[, 58]))
fst_UZ_auto <- cbind(FST_m_auto[, 1:5], (FST_m_auto[, 45]), (FST_m_auto[, 39]), (FST_m_auto[, 32]), (FST_m_auto[, 24]), (FST_m_auto[, 50]), (FST_m_auto[, 57]), (FST_m_auto[, 59]), (FST_m_auto[, 60]))
fst_UZ_z <- cbind(FST_m_z[, 1:5], (FST_m_z[, 45]), (FST_m_z[, 39]), (FST_m_z[, 32]), (FST_m_z[, 24]), (FST_m_z[, 50]), (FST_m_z[, 57]), (FST_m_z[, 59]), (FST_m_z[, 60]))
mean(as.matrix(fst_UZ_z[, 6:13]))
test_mat <- as.matrix(fst_UZ_z[,6:13])
mean(as.matrix(fst_UZ_auto[, 6:13]))


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
#mean(fst_noUZ_auto_filtered$fraction_covered)
fst_noUZ_auto_filtered <- fst_noUZ_auto[fst_noUZ_auto$nSNPs >  20, ]
fst_noUZ_auto_filtered <- fst_noUZ_auto_filtered[fst_noUZ_auto_filtered$fraction_covered >  0.2, ]

hist(as.numeric(as.matrix(fst_noUZ_auto[, 7:34])))
mean(as.numeric(as.matrix(fst_noUZ_auto[, 7:34])))
mean(as.numeric(as.matrix(fst_noUZ_auto_filtered[, 7:34])))
pdf("auto_fst.pdf", height = 4, width = 4)
hist(as.numeric(as.matrix(fst_noUZ_auto_filtered[, 7:34])), xlim = c(0,1), col = "black", main = "Distribution of autosomal FST values", xlab = "FST")
dev.off()

fst_auto_g_pthree <- fst_noUZ_auto_filtered[fst_noUZ_auto_filtered[, 7:34] >= 0.25, ]
pdf("auto_fst_zoomed.pdf", height = 4, width = 6)
hist(as.numeric(as.matrix(fst_noUZ_auto_filtered[, 7:34])), xlim = c(0.4, 0.6), ylim = c(0,40), col = 'black', main = "", xlab = "", ylab = "", breaks = 50)
dev.off()

mean(as.matrix(fst_noUZ_auto_filtered[, 7:34]))
mean(as.matrix(fst_noUZ_z_filtered[, 7:34]))

t.test(as.matrix(fst_noUZ_auto_filtered[, 7:34]))
t.test(as.matrix(fst_noUZ_z_filtered[, 7:34]))


##z

fst_noUZ_z <- cbind(FST_m_z[, 1:5], (FST_m_z[, 33]), (FST_m_z[, 26]), (FST_m_z[, 18]), (FST_m_z[, 40]), (FST_m_z[, 42]), (FST_m_z[, 43]), (FST_m_z[, 44]), (FST_m_z[, 25]), (FST_m_z[, 17]), (FST_m_z[, 34]), (FST_m_z[, 36]), (FST_m_z[, 37]), (FST_m_z[, 38]), (FST_m_z[, 16]), (FST_m_z[, 27]), (FST_m_z[, 29]), (FST_m_z[, 30]), (FST_m_z[, 31]), (FST_m_z[, 19]), (FST_m_z[, 21]), (FST_m_z[, 22]), (FST_m_z[, 23]), (FST_m_z[, 47]), (FST_m_z[, 48]), (FST_m_z[, 49]), (FST_m_z[, 55]), (FST_m_z[, 56]), (FST_m_z[, 58]))

fst_noUZ_z <- separate(data = fst_noUZ_z, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
fst_noUZ_z$Scaffold<-gsub("v02","",as.character(fst_noUZ_z$Scaffold))
fst_noUZ_z$length<-gsub("v02","",as.character(fst_noUZ_z$length))
fst_noUZ_z$length<-replace_na(fst_noUZ_z$length,"414000")
colnames(fst_noUZ_z) <- c("Scaffold", "Length", "Window", "nSNPs", "fraction_covered", "coverage_ave",
                         "N75_N125", "N75_N175", "N75_N225", "N75_S25", "N75_S125", 
                         "N75_S175", "N75_S225", "N125_N175", "N125_N225", "N125_S25", 
                         "N125_S125", "N125_S175", "N125_S225", "N175_N225", "N175_S25", 
                          "N175_S125", "N175_S175", "N175_S225", "N225_S25", 
                          "N225_S125", "N225_S175", "N225_S225", 
                         "S25_S125", "S25_S175", "S25_S225", 
                        "S125_S175", "S125_S225",
                         "S175_S225")

fst_noUZ_z_filtered <- fst_noUZ_z[fst_noUZ_z$nSNPs >  20, ]
fst_noUZ_z_filtered <- fst_noUZ_z_filtered[fst_noUZ_z_filtered$fraction_covered >  0.2, ]

mean(fst_noUZ_z_filtered$fraction_covered)


hist(as.numeric(as.matrix(fst_noUZ_z[, 7:34])))
mean(as.numeric(as.matrix(fst_noUZ_z[, 7:34])))
mean(as.numeric(as.matrix(fst_noUZ_z_filtered[, 7:34])))

pdf("z_fst.pdf", height = 4, width = 4)
hist(as.numeric(as.matrix(fst_noUZ_z_filtered[, 7:34])), xlim = c(0,1), col = "black", main = "Distribution of z-linked marker FST values", xlab = "FST")
dev.off()
pdf("z_fst_zoomed.pdf", height = 4, width = 6)
hist(as.numeric(as.matrix(fst_noUZ_z_filtered[, 7:34])), xlim = c(0.3, .4), ylim = c(0,40), col = 'black', main = "Distribution of z-linked marker FST values", xlab = "FST")
dev.off()


plot(fst_noUZ_auto_filtered$nSNPs, fst_noUZ_auto_filtered$fraction_covered, pch = 19, cex = 0.3)
#hist(fst_noUZ_auto, xlim = c(0.5, 1), ylim = c(0, 100))
max(fst_noUZ_auto_filtered[, 7:34])

plot(1:nrow(fst_noUZ_auto_filtered), fst_noUZ_auto_filtered[, 7], pch = 19, cex = 0.2, ylim = c(0,1))
for(i in 8:ncol(fst_noUZ_auto_filtered)){
  points(1:nrow(fst_noUZ_auto_filtered), fst_noUZ_auto_filtered[, i], pch = 19, cex = 0.2)
  
}

which.max(fst_noUZ_auto[, 1])

nrow(fst_noUZ_auto)


fst_out_z <- matrix(nrow = 9, ncol = 9)

colnames(fst_out_z) <- c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ")
row.names(fst_out_z) <- c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225", "UZ")
diag(fst_out_z) <- 0
fst_out_z[1, 2] <- mean(FST_m_z[, 33])
fst_out_z[1, 3] <- mean(FST_m_z[, 26])
fst_out_z[1, 4] <- mean(FST_m_z[, 18])
fst_out_z[1, 5] <- mean(FST_m_z[, 40])
fst_out_z[1, 6] <- mean(FST_m_z[, 42])
fst_out_z[1, 7] <- mean(FST_m_z[, 43])
fst_out_z[1, 8] <- mean(FST_m_z[, 44])
fst_out_z[1, 9] <- mean(FST_m_z[, 45])
fst_out_z[2, 3] <- mean(FST_m_z[, 25])
fst_out_z[2, 4] <- mean(FST_m_z[, 17])
fst_out_z[2, 5] <- mean(FST_m_z[, 34])
fst_out_z[2, 6] <- mean(FST_m_z[, 36])
fst_out_z[2, 7] <- mean(FST_m_z[, 37])
fst_out_z[2, 8] <- mean(FST_m_z[, 38])
fst_out_z[2, 9] <- mean(FST_m_z[, 39])
fst_out_z[3, 4] <- mean(FST_m_z[, 16])
fst_out_z[3, 5] <- mean(FST_m_z[, 27])
fst_out_z[3, 6] <- mean(FST_m_z[, 29])
fst_out_z[3, 7] <- mean(FST_m_z[, 30])
fst_out_z[3, 8] <- mean(FST_m_z[, 31])
fst_out_z[3, 9] <- mean(FST_m_z[, 32])
fst_out_z[4, 5] <- mean(FST_m_z[, 19])
fst_out_z[4, 6] <- mean(FST_m_z[, 21])
fst_out_z[4, 7] <- mean(FST_m_z[, 22])
fst_out_z[4, 8] <- mean(FST_m_z[, 23])
fst_out_z[4, 9] <- mean(FST_m_z[, 24])
fst_out_z[5, 6] <- mean(FST_m_z[, 47])
fst_out_z[5, 7] <- mean(FST_m_z[, 48])
fst_out_z[5, 8] <- mean(FST_m_z[, 49])
fst_out_z[5, 9] <- mean(FST_m_z[, 50])
fst_out_z[6, 7] <- mean(FST_m_z[, 55])
fst_out_z[6, 8] <- mean(FST_m_z[, 56])
fst_out_z[6, 9] <- mean(FST_m_z[, 57])
fst_out_z[7, 8] <- mean(FST_m_z[, 58])
fst_out_z[7, 9] <- mean(FST_m_z[, 59])
fst_out_z[8, 9] <- mean(FST_m_z[, 60])

fst_out_z<- Matrix::forceSymmetric(fst_out_z)
fst_out_z <- as.matrix(fst_out_z)
mean(fst_out_z[1:8, 1:8])
mean(fst_out_z[9, 1:8])

t.test(fst_out_z[1:8, 1:8])
t.test(fst_out_z[9, 1:8])







fst_out_auto<- Matrix::forceSymmetric(fst_out_auto)
fst_out_auto <- as.matrix(fst_out_auto)
mean(fst_out_auto[1:8, 1:8])
mean(fst_out_auto[9, 1:8])
# Z

write.table(fst_out_auto, "~/Desktop/fst_auto_table.txt", col.names = TRUE)

write.table(fst_out_z, "~/Desktop/fst_z_table.txt", col.names = TRUE)

write.table(fst_dist_auto, "~/Desktop/dist_auto_table.txt", col.names = TRUE)

write.table(fst_dist_z, "~/Desktop/dist_z_table.txt", col.names = TRUE)

fst_dist_z <- fst_out_z / (1-fst_out_z)

distances <- read.table("../distances.txt", header = TRUE)
row.names(distances) <-(distances[, 1])
distances <- distances[, -1]

dist_z <- mantel(fst_dist_z[1:8, 1:8], distances)
dist_z_test <- mantel(fst_dist_z[1:7, 1:7], distances[1:7, 1:7])

fst_dist_auto <- fst_out_auto / (1-fst_out_auto)
dist_auto <- mantel(fst_dist_auto[1:8, 1:8], distances)
dist_auto_test <- mantel(fst_dist_auto[1:7, 1:7], distances[1:7, 1:7])


fst_dist_auto[1:7, 1:7]
heatmap(fst_dist_z[1:8, 1:8], Colv = "Rowv", Rowv = NA)
heatmap(fst_dist_auto[1:8, 1:8], Colv = "Rowv", Rowv = NA)


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


N225 <- separate(data = N225, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
N225$Scaffold<-gsub("v02","",as.character(N225$Scaffold))
N225$length<-gsub("v02","",as.character(N225$length))
N225$length<-replace_na(N225$length,"414000")


#u <- 2.9e-9

FST_m <- separate(data = FST_m, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
FST_m$Scaffold<-gsub("v02","",as.character(FST_m$Scaffold))
FST_m$length<-gsub("v02","",as.character(FST_m$length))
FST_m$length<-replace_na(FST_m$length,"414000")

ma <- function(x, n = 10){stats::filter(x, rep(1 / n, n), sides = 2)}

FST_period <- FST_m[which(FST_m$Scaffold == "scaffold532"), ]
plot(FST_period$Window, FST_period$N225_S75, cex = 0.2, pch = 19, ylim = c(0,1))

FST_period$N225_UZ_ma <- ma(FST_period$N225_UZ)
FST_period$N175_UZ_ma <- ma(FST_period$N175_UZ)
FST_period$N125_UZ_ma <- ma(FST_period$N125_UZ)
FST_period$N75_UZ_ma <- ma(FST_period$N75_UZ)
FST_period$S25_UZ_ma <- ma(FST_period$S25_UZ)
FST_period$S125_UZ_ma <- ma(FST_period$S125_UZ)
FST_period$S175_UZ_ma <- ma(FST_period$S175_UZ)
FST_period$S225_UZ_ma <- ma(FST_period$S225_UZ)

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

#period
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


test_vector <- abs(as.numeric(ord_pdfr_z_all[, 36 ]) - as.numeric(ord_pdfr_z_all[, 37]))

plot(1:length(test_vector), test_vector, pch = 19, cex = 0.5)

pdfr_frame <- as.matrix(cbind(ord_pdfr_z_all[, 2], test_vector))
pdfr_frame <- as.data.frame(pdfr_frame)
pdfr_frame$V1 <- as.numeric(as.character(pdfr_frame$V1))
pdfr_frame$test_vector <- as.numeric(as.character(pdfr_frame$test_vector))

plot(pdfr_frame$V1, pdfr_frame$test_vector, pch = 19, cex = 0.5)

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



FST_test <- FST_m[which(FST_m$Scaffold == "bab_gene_region_CC"), ]
plot(FST_test$Window, FST_test$N225_UZ, cex = 0.2, pch = 19, ylim = c(0,1))


mean(pdfr_fst[c(1:4, 6, 8:10)])
mean(period_fst[c(1:4, 6, 8:10)])
mean(pdfr_out_fst[c(1:4, 6, 8:10)])
mean(period_out_fst[c(1:4, 6, 8:10)])

t.test(pdfr_fst[c(1:4, 6, 8:10)], pdfr_out_fst[c(1:4, 6, 8:10)])
t.test(period_fst[c(1:4, 6, 8:10)], period_out_fst[c(1:4, 6, 8:10)])

t.test(pdfr_fst[c(1:4, 6, 8:10)], period_fst[c(1:4, 6, 8:10)])
t.test(pdfr_out_fst[c(1:4, 6, 8:10)], period_out_fst[c(1:4, 6, 8:10)])


#### pi

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


N175_pi <- read.table("~/Desktop/pi/N175_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
N175_pi <- separate(data = N175_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
N175_pi$Scaffold<-gsub("v02","",as.character(N175_pi$Scaffold))
N175_pi$length<-gsub("v02","",as.character(N175_pi$length))
N175_pi$length<-replace_na(N175_pi$length,"414000")
head(N175_pi)

N175_pi_period <- N175_pi[which(N175_pi$Scaffold == "scaffold532"), ]
N175_pi_pdfr <- N175_pi[which(N175_pi$Scaffold == "scaffold87"), ]

plot(N175_pi_period$Position, N175_pi_period$Pi, pch = 19, cex = .3)
plot(N175_pi_pdfr$Position, N175_pi_pdfr$Pi, pch = 19, cex = .3)

N125_pi <- read.table("~/Desktop/pi/N125_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
N125_pi <- separate(data = N125_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
N125_pi$Scaffold<-gsub("v02","",as.character(N125_pi$Scaffold))
N125_pi$length<-gsub("v02","",as.character(N125_pi$length))
N125_pi$length<-replace_na(N125_pi$length,"414000")
head(N125_pi)

N125_pi_period <- N125_pi[which(N125_pi$Scaffold == "scaffold532"), ]
N125_pi_pdfr <- N125_pi[which(N125_pi$Scaffold == "scaffold87"), ]

plot(N125_pi_period$Position, N125_pi_period$Pi, pch = 19, cex = .3)
plot(N125_pi_pdfr$Position, N125_pi_pdfr$Pi, pch = 19, cex = .3)

N75_pi <- read.table("~/Desktop/pi/N75_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
N75_pi <- separate(data = N75_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
N75_pi$Scaffold<-gsub("v02","",as.character(N75_pi$Scaffold))
N75_pi$length<-gsub("v02","",as.character(N75_pi$length))
N75_pi$length<-replace_na(N75_pi$length,"414000")
head(N75_pi)

N75_pi_period <- N75_pi[which(N75_pi$Scaffold == "scaffold532"), ]
N75_pi_pdfr <- N75_pi[which(N75_pi$Scaffold == "scaffold87"), ]

plot(N75_pi_period$Position, N75_pi_period$Pi, pch = 19, cex = .3)
plot(N75_pi_pdfr$Position, N75_pi_pdfr$Pi, pch = 19, cex = .3)

N25_pi <- read.table("~/Desktop/pi/N25_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
N25_pi <- separate(data = N25_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
N25_pi$Scaffold<-gsub("v02","",as.character(N25_pi$Scaffold))
N25_pi$length<-gsub("v02","",as.character(N25_pi$length))
N25_pi$length<-replace_na(N25_pi$length,"414000")
head(N25_pi)

N25_pi_period <- N25_pi[which(N25_pi$Scaffold == "scaffold532"), ]
N25_pi_pdfr <- N25_pi[which(N25_pi$Scaffold == "scaffold87"), ]

plot(N25_pi_period$Position, N25_pi_period$Pi, pch = 19, cex = .3)
plot(N25_pi_pdfr$Position, N25_pi_pdfr$Pi, pch = 19, cex = .3)

S25_pi <- read.table("~/Desktop/pi/S25_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
S25_pi <- separate(data = S25_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
S25_pi$Scaffold<-gsub("v02","",as.character(S25_pi$Scaffold))
S25_pi$length<-gsub("v02","",as.character(S25_pi$length))
S25_pi$length<-replace_na(S25_pi$length,"414000")
head(S25_pi)

S25_pi_period <- S25_pi[which(S25_pi$Scaffold == "scaffold532"), ]
S25_pi_pdfr <- S25_pi[which(S25_pi$Scaffold == "scaffold87"), ]

plot(S25_pi_period$Position, S25_pi_period$Pi, pch = 19, cex = .3)
plot(S25_pi_pdfr$Position, S25_pi_pdfr$Pi, pch = 19, cex = .3)


S75_pi <- read.table("~/Desktop/pi/S75_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
S75_pi <- separate(data = S75_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
S75_pi$Scaffold<-gsub("v02","",as.character(S75_pi$Scaffold))
S75_pi$length<-gsub("v02","",as.character(S75_pi$length))
S75_pi$length<-replace_na(S75_pi$length,"414000")
head(S75_pi)

S75_pi_period <- S75_pi[which(S75_pi$Scaffold == "scaffold532"), ]
S75_pi_pdfr <- S75_pi[which(S75_pi$Scaffold == "scaffold87"), ]

plot(S75_pi_period$Position, S75_pi_period$Pi, pch = 19, cex = .3)
plot(S75_pi_pdfr$Position, S75_pi_pdfr$Pi, pch = 19, cex = .3)

S125_pi <- read.table("~/Desktop/pi/S125_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
S125_pi <- separate(data = S125_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
S125_pi$Scaffold<-gsub("v02","",as.character(S125_pi$Scaffold))
S125_pi$length<-gsub("v02","",as.character(S125_pi$length))
S125_pi$length<-replace_na(S125_pi$length,"414000")
head(S125_pi)

S125_pi_period <- S125_pi[which(S125_pi$Scaffold == "scaffold532"), ]
S125_pi_pdfr <- S125_pi[which(S125_pi$Scaffold == "scaffold87"), ]

plot(S125_pi_period$Position, S125_pi_period$Pi, pch = 19, cex = .3)
plot(S125_pi_pdfr$Position, S125_pi_pdfr$Pi, pch = 19, cex = .3)

S175_pi <- read.table("~/Desktop/pi/S175_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
S175_pi <- separate(data = S175_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
S175_pi$Scaffold<-gsub("v02","",as.character(S175_pi$Scaffold))
S175_pi$length<-gsub("v02","",as.character(S175_pi$length))
S175_pi$length<-replace_na(S175_pi$length,"414000")
head(S175_pi)

S175_pi_period <- S175_pi[which(S175_pi$Scaffold == "scaffold532"), ]
S175_pi_pdfr <- S175_pi[which(S175_pi$Scaffold == "scaffold87"), ]

plot(S175_pi_period$Position, S175_pi_period$Pi, pch = 19, cex = .3)
plot(S175_pi_pdfr$Position, S175_pi_pdfr$Pi, pch = 19, cex = .3)

S225_pi <- read.table("~/Desktop/pi/S225_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
S225_pi <- separate(data = S225_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
S225_pi$Scaffold<-gsub("v02","",as.character(S225_pi$Scaffold))
S225_pi$length<-gsub("v02","",as.character(S225_pi$length))
S225_pi$length<-replace_na(S225_pi$length,"414000")
head(S225_pi)

S225_pi_period <- S225_pi[which(S225_pi$Scaffold == "scaffold532"), ]
S225_pi_pdfr <- S225_pi[which(S225_pi$Scaffold == "scaffold87"), ]

plot(S225_pi_period$Position, S225_pi_period$Pi, pch = 19, cex = .3)
plot(S225_pi_pdfr$Position, S225_pi_pdfr$Pi, pch = 19, cex = .3)

UZ1_pi <- read.table("~/Desktop/pi/UZ1_pi1filtered.kb.txt", col.names = c("Scaffold", "Position", "NSNPs", "Coverage", "Pi"))
UZ1_pi <- separate(data = UZ1_pi, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
UZ1_pi$Scaffold<-gsub("v02","",as.character(UZ1_pi$Scaffold))
UZ1_pi$length<-gsub("v02","",as.character(UZ1_pi$length))
UZ1_pi$length<-replace_na(UZ1_pi$length,"414000")
head(UZ1_pi)

UZ1_pi_period <- UZ1_pi[which(UZ1_pi$Scaffold == "scaffold532"), ]
UZ1_pi_pdfr <- UZ1_pi[which(UZ1_pi$Scaffold == "scaffold87"), ]

plot(UZ1_pi_period$Position, UZ1_pi_period$Pi, pch = 19, cex = .3)
plot(UZ1_pi_pdfr$Position, UZ1_pi_pdfr$Pi, pch = 19, cex = .3)


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


UZ1_exgenic_pi <- mean(UZ1_pi_period[UZ1_pi_period$Position < 75000 | UZ1_pi_period$Position > 145000, ]$Pi, na.rm = TRUE  )
UZ1_genic_pi <- mean(UZ1_pi_period[UZ1_pi_period$Position > 75000 & UZ1_pi_period$Position < 145000, ]$Pi, na.rm = TRUE  )
N225_exgenic_pi <- mean(N225_pi_period[N225_pi_period$Position < 75000 | N225_pi_period$Position > 145000, ]$Pi, na.rm = TRUE  )
N225_genic_pi <- mean(N225_pi_period[N225_pi_period$Position > 75000 & N225_pi_period$Position < 145000, ]$Pi, na.rm = TRUE  )
N175_exgenic_pi <- mean(N175_pi_period[N175_pi_period$Position < 75000 | N175_pi_period$Position > 145000, ]$Pi, na.rm = TRUE  )
N175_genic_pi <- mean(N175_pi_period[N175_pi_period$Position > 75000 & N175_pi_period$Position < 145000, ]$Pi, na.rm = TRUE  )
N125_exgenic_pi <- mean(N125_pi_period[N125_pi_period$Position < 75000 | N125_pi_period$Position > 145000, ]$Pi, na.rm = TRUE  )
N125_genic_pi <- mean(N125_pi_period[N125_pi_period$Position > 75000 & N125_pi_period$Position < 145000, ]$Pi, na.rm = TRUE  )
N75_exgenic_pi <- mean(N75_pi_period[N75_pi_period$Position < 75000 | N75_pi_period$Position > 145000, ]$Pi, na.rm = TRUE  )
N75_genic_pi <- mean(N75_pi_period[N75_pi_period$Position > 75000 & N75_pi_period$Position < 145000, ]$Pi, na.rm = TRUE  )
S25_exgenic_pi <- mean(S25_pi_period[S25_pi_period$Position < 75000 | S25_pi_period$Position > 145000, ]$Pi, na.rm = TRUE  )
S25_genic_pi <- mean(S25_pi_period[S25_pi_period$Position > 75000 & S25_pi_period$Position < 145000, ]$Pi, na.rm = TRUE  )
S125_exgenic_pi <- mean(S125_pi_period[S125_pi_period$Position < 75000 | S125_pi_period$Position > 145000, ]$Pi, na.rm = TRUE  )
S125_genic_pi <- mean(S125_pi_period[S125_pi_period$Position > 75000 & S125_pi_period$Position < 145000, ]$Pi, na.rm = TRUE  )
S175_exgenic_pi <- mean(S175_pi_period[S175_pi_period$Position < 75000 | S175_pi_period$Position > 145000, ]$Pi, na.rm = TRUE  )
S175_genic_pi <- mean(S175_pi_period[S175_pi_period$Position > 75000 & S175_pi_period$Position < 145000, ]$Pi, na.rm = TRUE  )
S225_exgenic_pi <- mean(S225_pi_period[S225_pi_period$Position < 75000 | S225_pi_period$Position > 145000, ]$Pi, na.rm = TRUE  )
S225_genic_pi <- mean(S225_pi_period[S225_pi_period$Position > 75000 & S225_pi_period$Position < 145000, ]$Pi, na.rm = TRUE  )

period_genic_pi <- c(N225_genic_pi, N175_genic_pi, N125_genic_pi, N75_genic_pi, S25_genic_pi, S125_genic_pi, S175_genic_pi, S225_genic_pi)
period_exgenic_pi <- c(N225_exgenic_pi, N175_exgenic_pi, N125_exgenic_pi, N75_exgenic_pi, S25_exgenic_pi, S125_exgenic_pi, S175_exgenic_pi, S225_exgenic_pi)
mean(period_genic_pi)
mean(period_exgenic_pi)
t.test(period_genic_pi, period_exgenic_pi)

t.test(period_genic_pi)
t.test(period_exgenic_pi)




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

N175_D_pdfr <- N175_D[N175_D$Scaffold == 'scaffold87|size517236', ]
plot(N175_D_pdfr$Position, N175_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
N125_D_pdfr <- N125_D[N125_D$Scaffold == 'scaffold87|size517236', ]
plot(N125_D_pdfr$Position, N125_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
N75_D_pdfr <- N75_D[N75_D$Scaffold == 'scaffold87|size517236', ]
plot(N75_D_pdfr$Position, N75_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
N25_D_pdfr <- N25_D[N25_D$Scaffold == 'scaffold87|size517236', ]
plot(N25_D_pdfr$Position, N25_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S25_D_pdfr <- S25_D[S25_D$Scaffold == 'scaffold87|size517236', ]
plot(S25_D_pdfr$Position, S25_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S75_D_pdfr <- S75_D[S75_D$Scaffold == 'scaffold87|size517236', ]
plot(S75_D_pdfr$Position, S75_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S125_D_pdfr <- S125_D[S125_D$Scaffold == 'scaffold87|size517236', ]
plot(S125_D_pdfr$Position, S125_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S175_D_pdfr <- S175_D[S175_D$Scaffold == 'scaffold87|size517236', ]
plot(S175_D_pdfr$Position, S175_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S225_D_pdfr <- S225_D[S225_D$Scaffold == 'scaffold87|size517236', ]
plot(S225_D_pdfr$Position, S225_D_pdfr$D, pch = 19, cex = 0.7, ylim = c(-3, 3))


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

N175_D_period <- N175_D[N175_D$Scaffold == 'scaffold532|size371787', ]
plot(N175_D_period$Position, N175_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
N125_D_period <- N125_D[N125_D$Scaffold == 'scaffold532|size371787', ]
plot(N125_D_period$Position, N125_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
N75_D_period <- N75_D[N75_D$Scaffold == 'scaffold532|size371787', ]
plot(N75_D_period$Position, N75_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
N25_D_period <- N25_D[N25_D$Scaffold == 'scaffold532|size371787', ]
plot(N25_D_period$Position, N25_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S25_D_period <- S25_D[S25_D$Scaffold == 'scaffold532|size371787', ]
plot(S25_D_period$Position, S25_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S75_D_period <- S75_D[S75_D$Scaffold == 'scaffold532|size371787', ]
plot(S75_D_period$Position, S75_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S125_D_period <- S125_D[S125_D$Scaffold == 'scaffold532|size371787', ]
plot(S125_D_period$Position, S125_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S175_D_period <- S175_D[S175_D$Scaffold == 'scaffold532|size371787', ]
plot(S175_D_period$Position, S175_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3))
S225_D_period <- S225_D[S225_D$Scaffold == 'scaffold532|size371787', ]
plot(S225_D_period$Position, S225_D_period$D, pch = 19, cex = 0.7, ylim = c(-3, 3))





# distinguishing taj D

#UZ1_exgenic_D <- mean(UZ1_D_period[UZ1_D_period$Position < 75000 | UZ1_D_period$Position > 145000, ]$D, na.rm = TRUE  )
#UZ1_genic_D <- mean(UZ1_D_period[UZ1_D_period$Position > 75000 & UZ1_D_period$Position < 145000, ]$D, na.rm = TRUE  )
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

pdfr_genic_D <- c(N225_genic_D, N175_genic_D, N125_genic_D, N75_genic_D, S25_genic_D, S125_genic_D, S175_genic_D, S225_genic_D)
pdfr_exgenic_D <- c(N225_exgenic_D, N175_exgenic_D, N125_exgenic_D, N75_exgenic_D, S25_exgenic_D, S125_exgenic_D, S175_exgenic_D, S225_exgenic_D)
mean(pdfr_genic_D)
mean(pdfr_exgenic_D)
t.test(pdfr_genic_D, pdfr_exgenic_D)

t.test(pdfr_genic_D)
t.test(pdfr_exgenic_D)





### pdfr plot

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
pdf("pdfr_fst_N175_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$N175_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$N175_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("pdfr_fst_N125_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$N125_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$N125_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("pdfr_fst_N75_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$N75_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$N75_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("pdfr_fst_N25_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$N25_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$N25_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("pdfr_fst_S25_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$S25_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$S25_UZ_ma, type = 'l', col = 'red')
points(FST_pdfr$Window, FST_pdfr$S25_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("pdfr_fst_S75_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$S75_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$S75_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("pdfr_fst_S125_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$S125_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$S125_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("pdfr_fst_S175_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$S175_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$S175_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()
pdf("pdfr_fst_S225_UZ.pdf", height = 1, width = 2.33)
op <- par(mar = rep(0, 4))
plot(FST_pdfr$Window, FST_pdfr$S225_UZ, cex = 0.25, pch = 19, ylim = c(0,1))
points(FST_pdfr$Window, FST_pdfr$S225_UZ_ma, type = 'l', col = 'red')
par(op)
dev.off()


library(dplyr)

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









## scaffold178
FST_scaffold178 <- FST_m[which(FST_m$Scaffold == "scaffold178"), ]
FST_scaffold178 <- FST_scaffold178[sample(nrow(FST_scaffold178), 306), ]

UZ_scaffold178 <- FST_scaffold178[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold178 <- FST_scaffold178[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold178 <- mutate_all(UZ_scaffold178, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold178 <- mutate_all(Other_comparisons_scaffold178, function(x) as.numeric(as.character(x)))

UZ_scaffold178 <- cbind(UZ_scaffold178, FST_scaffold178$Window)
UZ_scaffold178$color <- "black"

UZ_scaffold178_ordered <- UZ_scaffold178[order(-UZ_scaffold178$N75_UZ),]
UZ_scaffold178$observation <- 1:nrow(UZ_scaffold178)
UZ_scaffold178_ordered$observation <- 1:nrow(UZ_scaffold178_ordered)
plot(UZ_scaffold178_ordered$observation, UZ_scaffold178_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold178_ordered$color)
UZ_scaffold178_ordered_N125 <- UZ_scaffold178[order(-UZ_scaffold178$N125_UZ),]
points(1:nrow(UZ_scaffold178_ordered_N125), UZ_scaffold178_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold178_ordered$color, pch = 19, cex = .5)
UZ_scaffold178_ordered_N175 <- UZ_scaffold178[order(-UZ_scaffold178$N175_UZ),]
points(1:nrow(UZ_scaffold178_ordered_N175), UZ_scaffold178_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold178_ordered$color, pch = 19, cex = .5)
UZ_scaffold178_ordered_N225 <- UZ_scaffold178[order(-UZ_scaffold178$N225_UZ),]
points(1:nrow(UZ_scaffold178_ordered_N225), UZ_scaffold178_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold178_ordered$color, pch = 19, cex = .5)
UZ_scaffold178_ordered_S25 <- UZ_scaffold178[order(-UZ_scaffold178$S25_UZ),]
points(1:nrow(UZ_scaffold178_ordered_S25), UZ_scaffold178_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold178_ordered$color, pch = 19, cex = .5)
UZ_scaffold178_ordered_S125 <- UZ_scaffold178[order(-UZ_scaffold178$S125_UZ),]
points(1:nrow(UZ_scaffold178_ordered_S125), UZ_scaffold178_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold178_ordered$color, pch = 19, cex = .5)
UZ_scaffold178_ordered_S175 <- UZ_scaffold178[order(-UZ_scaffold178$S175_UZ),]
points(1:nrow(UZ_scaffold178_ordered_S175), UZ_scaffold178_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold178_ordered$color, pch = 19, cex = .5)
UZ_scaffold178_ordered_S225 <- UZ_scaffold178[order(-UZ_scaffold178$S225_UZ),]
points(1:nrow(UZ_scaffold178_ordered_S225), UZ_scaffold178_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold178_ordered$color, pch = 19, cex = .5)




## scaffold1
FST_scaffold1 <- FST_m[which(FST_m$Scaffold == "scaffold1"), ]
FST_scaffold1 <- FST_scaffold1[sample(nrow(FST_scaffold1), 306), ]

UZ_scaffold1 <- FST_scaffold1[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold1 <- FST_scaffold1[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold1 <- mutate_all(UZ_scaffold1, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold1 <- mutate_all(Other_comparisons_scaffold1, function(x) as.numeric(as.character(x)))

UZ_scaffold1 <- cbind(UZ_scaffold1, FST_scaffold1$Window)
UZ_scaffold1$color <- "black"

UZ_scaffold1_ordered <- UZ_scaffold1[order(-UZ_scaffold1$N75_UZ),]
UZ_scaffold1$observation <- 1:nrow(UZ_scaffold1)
UZ_scaffold1_ordered$observation <- 1:nrow(UZ_scaffold1_ordered)
#plot(UZ_scaffold1_ordered$observation, UZ_scaffold1_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold1_ordered$color)
UZ_scaffold1_ordered_N125 <- UZ_scaffold1[order(-UZ_scaffold1$N125_UZ),]
points(1:nrow(UZ_scaffold1_ordered_N125), UZ_scaffold1_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold1_ordered$color, pch = 19, cex = .5)
UZ_scaffold1_ordered_N175 <- UZ_scaffold1[order(-UZ_scaffold1$N175_UZ),]
points(1:nrow(UZ_scaffold1_ordered_N175), UZ_scaffold1_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold1_ordered$color, pch = 19, cex = .5)
UZ_scaffold1_ordered_N225 <- UZ_scaffold1[order(-UZ_scaffold1$N225_UZ),]
points(1:nrow(UZ_scaffold1_ordered_N225), UZ_scaffold1_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold1_ordered$color, pch = 19, cex = .5)
UZ_scaffold1_ordered_S25 <- UZ_scaffold1[order(-UZ_scaffold1$S25_UZ),]
points(1:nrow(UZ_scaffold1_ordered_S25), UZ_scaffold1_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold1_ordered$color, pch = 19, cex = .5)
UZ_scaffold1_ordered_S125 <- UZ_scaffold1[order(-UZ_scaffold1$S125_UZ),]
points(1:nrow(UZ_scaffold1_ordered_S125), UZ_scaffold1_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold1_ordered$color, pch = 19, cex = .5)
UZ_scaffold1_ordered_S175 <- UZ_scaffold1[order(-UZ_scaffold1$S175_UZ),]
points(1:nrow(UZ_scaffold1_ordered_S175), UZ_scaffold1_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold1_ordered$color, pch = 19, cex = .5)
UZ_scaffold1_ordered_S225 <- UZ_scaffold1[order(-UZ_scaffold1$S225_UZ),]
points(1:nrow(UZ_scaffold1_ordered_S225), UZ_scaffold1_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold1_ordered$color, pch = 19, cex = .5)

## scaffold2
FST_scaffold2 <- FST_m[which(FST_m$Scaffold == "scaffold2"), ]
FST_scaffold2 <- FST_scaffold2[sample(nrow(FST_scaffold2), 306), ]

UZ_scaffold2 <- FST_scaffold2[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold2 <- FST_scaffold2[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold2 <- mutate_all(UZ_scaffold2, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold2 <- mutate_all(Other_comparisons_scaffold2, function(x) as.numeric(as.character(x)))

UZ_scaffold2 <- cbind(UZ_scaffold2, FST_scaffold2$Window)
UZ_scaffold2$color <- "black"

UZ_scaffold2_ordered <- UZ_scaffold2[order(-UZ_scaffold2$N75_UZ),]
UZ_scaffold2$observation <- 1:nrow(UZ_scaffold2)
UZ_scaffold2_ordered$observation <- 1:nrow(UZ_scaffold2_ordered)
#plot(UZ_scaffold2_ordered$observation, UZ_scaffold2_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold2_ordered$color)
UZ_scaffold2_ordered_N125 <- UZ_scaffold2[order(-UZ_scaffold2$N125_UZ),]
points(1:nrow(UZ_scaffold2_ordered_N125), UZ_scaffold2_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold2_ordered$color, pch = 19, cex = .5)
UZ_scaffold2_ordered_N175 <- UZ_scaffold2[order(-UZ_scaffold2$N175_UZ),]
points(1:nrow(UZ_scaffold2_ordered_N175), UZ_scaffold2_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold2_ordered$color, pch = 19, cex = .5)
UZ_scaffold2_ordered_N225 <- UZ_scaffold2[order(-UZ_scaffold2$N225_UZ),]
points(1:nrow(UZ_scaffold2_ordered_N225), UZ_scaffold2_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold2_ordered$color, pch = 19, cex = .5)
UZ_scaffold2_ordered_S25 <- UZ_scaffold2[order(-UZ_scaffold2$S25_UZ),]
points(1:nrow(UZ_scaffold2_ordered_S25), UZ_scaffold2_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold2_ordered$color, pch = 19, cex = .5)
UZ_scaffold2_ordered_S125 <- UZ_scaffold2[order(-UZ_scaffold2$S125_UZ),]
points(1:nrow(UZ_scaffold2_ordered_S125), UZ_scaffold2_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold2_ordered$color, pch = 19, cex = .5)
UZ_scaffold2_ordered_S175 <- UZ_scaffold2[order(-UZ_scaffold2$S175_UZ),]
points(1:nrow(UZ_scaffold2_ordered_S175), UZ_scaffold2_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold2_ordered$color, pch = 19, cex = .5)
UZ_scaffold2_ordered_S225 <- UZ_scaffold2[order(-UZ_scaffold2$S225_UZ),]
points(1:nrow(UZ_scaffold2_ordered_S225), UZ_scaffold2_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold2_ordered$color, pch = 19, cex = .5)


## scaffold3
FST_scaffold3 <- FST_m[which(FST_m$Scaffold == "scaffold3"), ]
FST_scaffold3 <- FST_scaffold3[sample(nrow(FST_scaffold3), 306), ]

UZ_scaffold3 <- FST_scaffold3[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold3 <- FST_scaffold3[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold3 <- mutate_all(UZ_scaffold3, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold3 <- mutate_all(Other_comparisons_scaffold3, function(x) as.numeric(as.character(x)))

UZ_scaffold3 <- cbind(UZ_scaffold3, FST_scaffold3$Window)
UZ_scaffold3$color <- "black"

UZ_scaffold3_ordered <- UZ_scaffold3[order(-UZ_scaffold3$N75_UZ),]
UZ_scaffold3$observation <- 1:nrow(UZ_scaffold3)
UZ_scaffold3_ordered$observation <- 1:nrow(UZ_scaffold3_ordered)
#plot(UZ_scaffold3_ordered$observation, UZ_scaffold3_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold3_ordered$color)
UZ_scaffold3_ordered_N125 <- UZ_scaffold3[order(-UZ_scaffold3$N125_UZ),]
points(1:nrow(UZ_scaffold3_ordered_N125), UZ_scaffold3_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold3_ordered$color, pch = 19, cex = .5)
UZ_scaffold3_ordered_N175 <- UZ_scaffold3[order(-UZ_scaffold3$N175_UZ),]
points(1:nrow(UZ_scaffold3_ordered_N175), UZ_scaffold3_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold3_ordered$color, pch = 19, cex = .5)
UZ_scaffold3_ordered_N225 <- UZ_scaffold3[order(-UZ_scaffold3$N225_UZ),]
points(1:nrow(UZ_scaffold3_ordered_N225), UZ_scaffold3_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold3_ordered$color, pch = 19, cex = .5)
UZ_scaffold3_ordered_S25 <- UZ_scaffold3[order(-UZ_scaffold3$S25_UZ),]
points(1:nrow(UZ_scaffold3_ordered_S25), UZ_scaffold3_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold3_ordered$color, pch = 19, cex = .5)
UZ_scaffold3_ordered_S125 <- UZ_scaffold3[order(-UZ_scaffold3$S125_UZ),]
points(1:nrow(UZ_scaffold3_ordered_S125), UZ_scaffold3_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold3_ordered$color, pch = 19, cex = .5)
UZ_scaffold3_ordered_S175 <- UZ_scaffold3[order(-UZ_scaffold3$S175_UZ),]
points(1:nrow(UZ_scaffold3_ordered_S175), UZ_scaffold3_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold3_ordered$color, pch = 19, cex = .5)
UZ_scaffold3_ordered_S225 <- UZ_scaffold3[order(-UZ_scaffold3$S225_UZ),]
points(1:nrow(UZ_scaffold3_ordered_S225), UZ_scaffold3_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold3_ordered$color, pch = 19, cex = .5)




## scaffold4
FST_scaffold4 <- FST_m[which(FST_m$Scaffold == "scaffold4"), ]
FST_scaffold4 <- FST_scaffold4[sample(nrow(FST_scaffold4), 306), ]

UZ_scaffold4 <- FST_scaffold4[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold4 <- FST_scaffold4[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold4 <- mutate_all(UZ_scaffold4, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold4 <- mutate_all(Other_comparisons_scaffold4, function(x) as.numeric(as.character(x)))

UZ_scaffold4 <- cbind(UZ_scaffold4, FST_scaffold4$Window)
UZ_scaffold4$color <- "black"

UZ_scaffold4_ordered <- UZ_scaffold4[order(-UZ_scaffold4$N75_UZ),]
UZ_scaffold4$observation <- 1:nrow(UZ_scaffold4)
UZ_scaffold4_ordered$observation <- 1:nrow(UZ_scaffold4_ordered)
#plot(UZ_scaffold4_ordered$observation, UZ_scaffold4_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold4_ordered$color)
UZ_scaffold4_ordered_N125 <- UZ_scaffold4[order(-UZ_scaffold4$N125_UZ),]
points(1:nrow(UZ_scaffold4_ordered_N125), UZ_scaffold4_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold4_ordered$color, pch = 19, cex = .5)
UZ_scaffold4_ordered_N175 <- UZ_scaffold4[order(-UZ_scaffold4$N175_UZ),]
points(1:nrow(UZ_scaffold4_ordered_N175), UZ_scaffold4_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold4_ordered$color, pch = 19, cex = .5)
UZ_scaffold4_ordered_N225 <- UZ_scaffold4[order(-UZ_scaffold4$N225_UZ),]
points(1:nrow(UZ_scaffold4_ordered_N225), UZ_scaffold4_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold4_ordered$color, pch = 19, cex = .5)
UZ_scaffold4_ordered_S25 <- UZ_scaffold4[order(-UZ_scaffold4$S25_UZ),]
points(1:nrow(UZ_scaffold4_ordered_S25), UZ_scaffold4_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold4_ordered$color, pch = 19, cex = .5)
UZ_scaffold4_ordered_S125 <- UZ_scaffold4[order(-UZ_scaffold4$S125_UZ),]
points(1:nrow(UZ_scaffold4_ordered_S125), UZ_scaffold4_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold4_ordered$color, pch = 19, cex = .5)
UZ_scaffold4_ordered_S175 <- UZ_scaffold4[order(-UZ_scaffold4$S175_UZ),]
points(1:nrow(UZ_scaffold4_ordered_S175), UZ_scaffold4_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold4_ordered$color, pch = 19, cex = .5)
UZ_scaffold4_ordered_S225 <- UZ_scaffold4[order(-UZ_scaffold4$S225_UZ),]
points(1:nrow(UZ_scaffold4_ordered_S225), UZ_scaffold4_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold4_ordered$color, pch = 19, cex = .5)


## scaffold5
FST_scaffold5 <- FST_m[which(FST_m$Scaffold == "scaffold5"), ]
FST_scaffold5 <- FST_scaffold5[sample(nrow(FST_scaffold5), 306), ]

UZ_scaffold5 <- FST_scaffold5[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold5 <- FST_scaffold5[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold5 <- mutate_all(UZ_scaffold5, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold5 <- mutate_all(Other_comparisons_scaffold5, function(x) as.numeric(as.character(x)))

UZ_scaffold5 <- cbind(UZ_scaffold5, FST_scaffold5$Window)
UZ_scaffold5$color <- "black"

UZ_scaffold5_ordered <- UZ_scaffold5[order(-UZ_scaffold5$N75_UZ),]
UZ_scaffold5$observation <- 1:nrow(UZ_scaffold5)
UZ_scaffold5_ordered$observation <- 1:nrow(UZ_scaffold5_ordered)
#plot(UZ_scaffold5_ordered$observation, UZ_scaffold5_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold5_ordered$color)
UZ_scaffold5_ordered_N125 <- UZ_scaffold5[order(-UZ_scaffold5$N125_UZ),]
points(1:nrow(UZ_scaffold5_ordered_N125), UZ_scaffold5_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold5_ordered$color, pch = 19, cex = .5)
UZ_scaffold5_ordered_N175 <- UZ_scaffold5[order(-UZ_scaffold5$N175_UZ),]
points(1:nrow(UZ_scaffold5_ordered_N175), UZ_scaffold5_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold5_ordered$color, pch = 19, cex = .5)
UZ_scaffold5_ordered_N225 <- UZ_scaffold5[order(-UZ_scaffold5$N225_UZ),]
points(1:nrow(UZ_scaffold5_ordered_N225), UZ_scaffold5_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold5_ordered$color, pch = 19, cex = .5)
UZ_scaffold5_ordered_S25 <- UZ_scaffold5[order(-UZ_scaffold5$S25_UZ),]
points(1:nrow(UZ_scaffold5_ordered_S25), UZ_scaffold5_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold5_ordered$color, pch = 19, cex = .5)
UZ_scaffold5_ordered_S125 <- UZ_scaffold5[order(-UZ_scaffold5$S125_UZ),]
points(1:nrow(UZ_scaffold5_ordered_S125), UZ_scaffold5_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold5_ordered$color, pch = 19, cex = .5)
UZ_scaffold5_ordered_S175 <- UZ_scaffold5[order(-UZ_scaffold5$S175_UZ),]
points(1:nrow(UZ_scaffold5_ordered_S175), UZ_scaffold5_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold5_ordered$color, pch = 19, cex = .5)
UZ_scaffold5_ordered_S225 <- UZ_scaffold5[order(-UZ_scaffold5$S225_UZ),]
points(1:nrow(UZ_scaffold5_ordered_S225), UZ_scaffold5_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold5_ordered$color, pch = 19, cex = .5)


## scaffold6
FST_scaffold6 <- FST_m[which(FST_m$Scaffold == "scaffold6"), ]
FST_scaffold6 <- FST_scaffold6[sample(nrow(FST_scaffold6), 306), ]

UZ_scaffold6 <- FST_scaffold6[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold6 <- FST_scaffold6[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold6 <- mutate_all(UZ_scaffold6, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold6 <- mutate_all(Other_comparisons_scaffold6, function(x) as.numeric(as.character(x)))

UZ_scaffold6 <- cbind(UZ_scaffold6, FST_scaffold6$Window)
UZ_scaffold6$color <- "black"

UZ_scaffold6_ordered <- UZ_scaffold6[order(-UZ_scaffold6$N75_UZ),]
UZ_scaffold6$observation <- 1:nrow(UZ_scaffold6)
UZ_scaffold6_ordered$observation <- 1:nrow(UZ_scaffold6_ordered)
#plot(UZ_scaffold6_ordered$observation, UZ_scaffold6_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold6_ordered$color)
UZ_scaffold6_ordered_N125 <- UZ_scaffold6[order(-UZ_scaffold6$N125_UZ),]
points(1:nrow(UZ_scaffold6_ordered_N125), UZ_scaffold6_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold6_ordered$color, pch = 19, cex = .5)
UZ_scaffold6_ordered_N175 <- UZ_scaffold6[order(-UZ_scaffold6$N175_UZ),]
points(1:nrow(UZ_scaffold6_ordered_N175), UZ_scaffold6_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold6_ordered$color, pch = 19, cex = .5)
UZ_scaffold6_ordered_N225 <- UZ_scaffold6[order(-UZ_scaffold6$N225_UZ),]
points(1:nrow(UZ_scaffold6_ordered_N225), UZ_scaffold6_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold6_ordered$color, pch = 19, cex = .5)
UZ_scaffold6_ordered_S25 <- UZ_scaffold6[order(-UZ_scaffold6$S25_UZ),]
points(1:nrow(UZ_scaffold6_ordered_S25), UZ_scaffold6_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold6_ordered$color, pch = 19, cex = .5)
UZ_scaffold6_ordered_S125 <- UZ_scaffold6[order(-UZ_scaffold6$S125_UZ),]
points(1:nrow(UZ_scaffold6_ordered_S125), UZ_scaffold6_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold6_ordered$color, pch = 19, cex = .5)
UZ_scaffold6_ordered_S175 <- UZ_scaffold6[order(-UZ_scaffold6$S175_UZ),]
points(1:nrow(UZ_scaffold6_ordered_S175), UZ_scaffold6_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold6_ordered$color, pch = 19, cex = .5)
UZ_scaffold6_ordered_S225 <- UZ_scaffold6[order(-UZ_scaffold6$S225_UZ),]
points(1:nrow(UZ_scaffold6_ordered_S225), UZ_scaffold6_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold6_ordered$color, pch = 19, cex = .5)



## scaffold7
FST_scaffold7 <- FST_m[which(FST_m$Scaffold == "scaffold7"), ]
FST_scaffold7 <- FST_scaffold7[sample(nrow(FST_scaffold7), 306), ]

UZ_scaffold7 <- FST_scaffold7[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold7 <- FST_scaffold7[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold7 <- mutate_all(UZ_scaffold7, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold7 <- mutate_all(Other_comparisons_scaffold7, function(x) as.numeric(as.character(x)))

UZ_scaffold7 <- cbind(UZ_scaffold7, FST_scaffold7$Window)
UZ_scaffold7$color <- "black"

UZ_scaffold7_ordered <- UZ_scaffold7[order(-UZ_scaffold7$N75_UZ),]
UZ_scaffold7$observation <- 1:nrow(UZ_scaffold7)
UZ_scaffold7_ordered$observation <- 1:nrow(UZ_scaffold7_ordered)
#plot(UZ_scaffold7_ordered$observation, UZ_scaffold7_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold7_ordered$color)
UZ_scaffold7_ordered_N125 <- UZ_scaffold7[order(-UZ_scaffold7$N125_UZ),]
points(1:nrow(UZ_scaffold7_ordered_N125), UZ_scaffold7_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold7_ordered$color, pch = 19, cex = .5)
UZ_scaffold7_ordered_N175 <- UZ_scaffold7[order(-UZ_scaffold7$N175_UZ),]
points(1:nrow(UZ_scaffold7_ordered_N175), UZ_scaffold7_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold7_ordered$color, pch = 19, cex = .5)
UZ_scaffold7_ordered_N225 <- UZ_scaffold7[order(-UZ_scaffold7$N225_UZ),]
points(1:nrow(UZ_scaffold7_ordered_N225), UZ_scaffold7_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold7_ordered$color, pch = 19, cex = .5)
UZ_scaffold7_ordered_S25 <- UZ_scaffold7[order(-UZ_scaffold7$S25_UZ),]
points(1:nrow(UZ_scaffold7_ordered_S25), UZ_scaffold7_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold7_ordered$color, pch = 19, cex = .5)
UZ_scaffold7_ordered_S125 <- UZ_scaffold7[order(-UZ_scaffold7$S125_UZ),]
points(1:nrow(UZ_scaffold7_ordered_S125), UZ_scaffold7_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold7_ordered$color, pch = 19, cex = .5)
UZ_scaffold7_ordered_S175 <- UZ_scaffold7[order(-UZ_scaffold7$S175_UZ),]
points(1:nrow(UZ_scaffold7_ordered_S175), UZ_scaffold7_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold7_ordered$color, pch = 19, cex = .5)
UZ_scaffold7_ordered_S225 <- UZ_scaffold7[order(-UZ_scaffold7$S225_UZ),]
points(1:nrow(UZ_scaffold7_ordered_S225), UZ_scaffold7_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold7_ordered$color, pch = 19, cex = .5)



## scaffold8
FST_scaffold8 <- FST_m[which(FST_m$Scaffold == "scaffold8"), ]
FST_scaffold8 <- FST_scaffold8[sample(nrow(FST_scaffold8), 306), ]

UZ_scaffold8 <- FST_scaffold8[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold8 <- FST_scaffold8[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold8 <- mutate_all(UZ_scaffold8, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold8 <- mutate_all(Other_comparisons_scaffold8, function(x) as.numeric(as.character(x)))

UZ_scaffold8 <- cbind(UZ_scaffold8, FST_scaffold8$Window)
UZ_scaffold8$color <- "black"

UZ_scaffold8_ordered <- UZ_scaffold8[order(-UZ_scaffold8$N75_UZ),]
UZ_scaffold8$observation <- 1:nrow(UZ_scaffold8)
UZ_scaffold8_ordered$observation <- 1:nrow(UZ_scaffold8_ordered)
#plot(UZ_scaffold8_ordered$observation, UZ_scaffold8_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold8_ordered$color)
UZ_scaffold8_ordered_N125 <- UZ_scaffold8[order(-UZ_scaffold8$N125_UZ),]
points(1:nrow(UZ_scaffold8_ordered_N125), UZ_scaffold8_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold8_ordered$color, pch = 19, cex = .5)
UZ_scaffold8_ordered_N175 <- UZ_scaffold8[order(-UZ_scaffold8$N175_UZ),]
points(1:nrow(UZ_scaffold8_ordered_N175), UZ_scaffold8_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold8_ordered$color, pch = 19, cex = .5)
UZ_scaffold8_ordered_N225 <- UZ_scaffold8[order(-UZ_scaffold8$N225_UZ),]
points(1:nrow(UZ_scaffold8_ordered_N225), UZ_scaffold8_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold8_ordered$color, pch = 19, cex = .5)
UZ_scaffold8_ordered_S25 <- UZ_scaffold8[order(-UZ_scaffold8$S25_UZ),]
points(1:nrow(UZ_scaffold8_ordered_S25), UZ_scaffold8_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold8_ordered$color, pch = 19, cex = .5)
UZ_scaffold8_ordered_S125 <- UZ_scaffold8[order(-UZ_scaffold8$S125_UZ),]
points(1:nrow(UZ_scaffold8_ordered_S125), UZ_scaffold8_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold8_ordered$color, pch = 19, cex = .5)
UZ_scaffold8_ordered_S175 <- UZ_scaffold8[order(-UZ_scaffold8$S175_UZ),]
points(1:nrow(UZ_scaffold8_ordered_S175), UZ_scaffold8_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold8_ordered$color, pch = 19, cex = .5)
UZ_scaffold8_ordered_S225 <- UZ_scaffold8[order(-UZ_scaffold8$S225_UZ),]
points(1:nrow(UZ_scaffold8_ordered_S225), UZ_scaffold8_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold8_ordered$color, pch = 19, cex = .5)



## scaffold9
FST_scaffold9 <- FST_m[which(FST_m$Scaffold == "scaffold9"), ]
FST_scaffold9 <- FST_scaffold9[sample(nrow(FST_scaffold9), 306), ]

UZ_scaffold9 <- FST_scaffold9[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold9 <- FST_scaffold9[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold9 <- mutate_all(UZ_scaffold9, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold9 <- mutate_all(Other_comparisons_scaffold9, function(x) as.numeric(as.character(x)))

UZ_scaffold9 <- cbind(UZ_scaffold9, FST_scaffold9$Window)
UZ_scaffold9$color <- "black"

UZ_scaffold9_ordered <- UZ_scaffold9[order(-UZ_scaffold9$N75_UZ),]
UZ_scaffold9$observation <- 1:nrow(UZ_scaffold9)
UZ_scaffold9_ordered$observation <- 1:nrow(UZ_scaffold9_ordered)
#plot(UZ_scaffold9_ordered$observation, UZ_scaffold9_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold9_ordered$color)
UZ_scaffold9_ordered_N125 <- UZ_scaffold9[order(-UZ_scaffold9$N125_UZ),]
points(1:nrow(UZ_scaffold9_ordered_N125), UZ_scaffold9_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold9_ordered$color, pch = 19, cex = .5)
UZ_scaffold9_ordered_N175 <- UZ_scaffold9[order(-UZ_scaffold9$N175_UZ),]
points(1:nrow(UZ_scaffold9_ordered_N175), UZ_scaffold9_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold9_ordered$color, pch = 19, cex = .5)
UZ_scaffold9_ordered_N225 <- UZ_scaffold9[order(-UZ_scaffold9$N225_UZ),]
points(1:nrow(UZ_scaffold9_ordered_N225), UZ_scaffold9_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold9_ordered$color, pch = 19, cex = .5)
UZ_scaffold9_ordered_S25 <- UZ_scaffold9[order(-UZ_scaffold9$S25_UZ),]
points(1:nrow(UZ_scaffold9_ordered_S25), UZ_scaffold9_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold9_ordered$color, pch = 19, cex = .5)
UZ_scaffold9_ordered_S125 <- UZ_scaffold9[order(-UZ_scaffold9$S125_UZ),]
points(1:nrow(UZ_scaffold9_ordered_S125), UZ_scaffold9_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold9_ordered$color, pch = 19, cex = .5)
UZ_scaffold9_ordered_S175 <- UZ_scaffold9[order(-UZ_scaffold9$S175_UZ),]
points(1:nrow(UZ_scaffold9_ordered_S175), UZ_scaffold9_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold9_ordered$color, pch = 19, cex = .5)
UZ_scaffold9_ordered_S225 <- UZ_scaffold9[order(-UZ_scaffold9$S225_UZ),]
points(1:nrow(UZ_scaffold9_ordered_S225), UZ_scaffold9_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold9_ordered$color, pch = 19, cex = .5)



## scaffold10
FST_scaffold10 <- FST_m[which(FST_m$Scaffold == "scaffold10"), ]
FST_scaffold10 <- FST_scaffold10[sample(nrow(FST_scaffold10), 306), ]

UZ_scaffold10 <- FST_scaffold10[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold10 <- FST_scaffold10[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold10 <- mutate_all(UZ_scaffold10, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold10 <- mutate_all(Other_comparisons_scaffold10, function(x) as.numeric(as.character(x)))

UZ_scaffold10 <- cbind(UZ_scaffold10, FST_scaffold10$Window)
UZ_scaffold10$color <- "black"

UZ_scaffold10_ordered <- UZ_scaffold10[order(-UZ_scaffold10$N75_UZ),]
UZ_scaffold10$observation <- 1:nrow(UZ_scaffold10)
UZ_scaffold10_ordered$observation <- 1:nrow(UZ_scaffold10_ordered)
#plot(UZ_scaffold10_ordered$observation, UZ_scaffold10_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold10_ordered$color)
UZ_scaffold10_ordered_N125 <- UZ_scaffold10[order(-UZ_scaffold10$N125_UZ),]
points(1:nrow(UZ_scaffold10_ordered_N125), UZ_scaffold10_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold10_ordered$color, pch = 19, cex = .5)
UZ_scaffold10_ordered_N175 <- UZ_scaffold10[order(-UZ_scaffold10$N175_UZ),]
points(1:nrow(UZ_scaffold10_ordered_N175), UZ_scaffold10_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold10_ordered$color, pch = 19, cex = .5)
UZ_scaffold10_ordered_N225 <- UZ_scaffold10[order(-UZ_scaffold10$N225_UZ),]
points(1:nrow(UZ_scaffold10_ordered_N225), UZ_scaffold10_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold10_ordered$color, pch = 19, cex = .5)
UZ_scaffold10_ordered_S25 <- UZ_scaffold10[order(-UZ_scaffold10$S25_UZ),]
points(1:nrow(UZ_scaffold10_ordered_S25), UZ_scaffold10_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold10_ordered$color, pch = 19, cex = .5)
UZ_scaffold10_ordered_S125 <- UZ_scaffold10[order(-UZ_scaffold10$S125_UZ),]
points(1:nrow(UZ_scaffold10_ordered_S125), UZ_scaffold10_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold10_ordered$color, pch = 19, cex = .5)
UZ_scaffold10_ordered_S175 <- UZ_scaffold10[order(-UZ_scaffold10$S175_UZ),]
points(1:nrow(UZ_scaffold10_ordered_S175), UZ_scaffold10_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold10_ordered$color, pch = 19, cex = .5)
UZ_scaffold10_ordered_S225 <- UZ_scaffold10[order(-UZ_scaffold10$S225_UZ),]
points(1:nrow(UZ_scaffold10_ordered_S225), UZ_scaffold10_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold10_ordered$color, pch = 19, cex = .5)

## scaffold11
FST_scaffold11 <- FST_m[which(FST_m$Scaffold == "scaffold11"), ]
FST_scaffold11 <- FST_scaffold11[sample(nrow(FST_scaffold11), 306), ]

UZ_scaffold11 <- FST_scaffold11[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold11 <- FST_scaffold11[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold11 <- mutate_all(UZ_scaffold11, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold11 <- mutate_all(Other_comparisons_scaffold11, function(x) as.numeric(as.character(x)))

UZ_scaffold11 <- cbind(UZ_scaffold11, FST_scaffold11$Window)
UZ_scaffold11$color <- "black"

UZ_scaffold11_ordered <- UZ_scaffold11[order(-UZ_scaffold11$N75_UZ),]
UZ_scaffold11$observation <- 1:nrow(UZ_scaffold11)
UZ_scaffold11_ordered$observation <- 1:nrow(UZ_scaffold11_ordered)
#plot(UZ_scaffold11_ordered$observation, UZ_scaffold11_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold11_ordered$color)
UZ_scaffold11_ordered_N125 <- UZ_scaffold11[order(-UZ_scaffold11$N125_UZ),]
points(1:nrow(UZ_scaffold11_ordered_N125), UZ_scaffold11_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold11_ordered$color, pch = 19, cex = .5)
UZ_scaffold11_ordered_N175 <- UZ_scaffold11[order(-UZ_scaffold11$N175_UZ),]
points(1:nrow(UZ_scaffold11_ordered_N175), UZ_scaffold11_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold11_ordered$color, pch = 19, cex = .5)
UZ_scaffold11_ordered_N225 <- UZ_scaffold11[order(-UZ_scaffold11$N225_UZ),]
points(1:nrow(UZ_scaffold11_ordered_N225), UZ_scaffold11_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold11_ordered$color, pch = 19, cex = .5)
UZ_scaffold11_ordered_S25 <- UZ_scaffold11[order(-UZ_scaffold11$S25_UZ),]
points(1:nrow(UZ_scaffold11_ordered_S25), UZ_scaffold11_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold11_ordered$color, pch = 19, cex = .5)
UZ_scaffold11_ordered_S125 <- UZ_scaffold11[order(-UZ_scaffold11$S125_UZ),]
points(1:nrow(UZ_scaffold11_ordered_S125), UZ_scaffold11_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold11_ordered$color, pch = 19, cex = .5)
UZ_scaffold11_ordered_S175 <- UZ_scaffold11[order(-UZ_scaffold11$S175_UZ),]
points(1:nrow(UZ_scaffold11_ordered_S175), UZ_scaffold11_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold11_ordered$color, pch = 19, cex = .5)
UZ_scaffold11_ordered_S225 <- UZ_scaffold11[order(-UZ_scaffold11$S225_UZ),]
points(1:nrow(UZ_scaffold11_ordered_S225), UZ_scaffold11_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold11_ordered$color, pch = 19, cex = .5)

## scaffold12
FST_scaffold12 <- FST_m[which(FST_m$Scaffold == "scaffold12"), ]
FST_scaffold12 <- FST_scaffold12[sample(nrow(FST_scaffold12), 306), ]

UZ_scaffold12 <- FST_scaffold12[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold12 <- FST_scaffold12[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold12 <- mutate_all(UZ_scaffold12, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold12 <- mutate_all(Other_comparisons_scaffold12, function(x) as.numeric(as.character(x)))

UZ_scaffold12 <- cbind(UZ_scaffold12, FST_scaffold12$Window)
UZ_scaffold12$color <- "black"

UZ_scaffold12_ordered <- UZ_scaffold12[order(-UZ_scaffold12$N75_UZ),]
UZ_scaffold12$observation <- 1:nrow(UZ_scaffold12)
UZ_scaffold12_ordered$observation <- 1:nrow(UZ_scaffold12_ordered)
#plot(UZ_scaffold12_ordered$observation, UZ_scaffold12_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold12_ordered$color)
UZ_scaffold12_ordered_N125 <- UZ_scaffold12[order(-UZ_scaffold12$N125_UZ),]
points(1:nrow(UZ_scaffold12_ordered_N125), UZ_scaffold12_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold12_ordered$color, pch = 19, cex = .5)
UZ_scaffold12_ordered_N175 <- UZ_scaffold12[order(-UZ_scaffold12$N175_UZ),]
points(1:nrow(UZ_scaffold12_ordered_N175), UZ_scaffold12_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold12_ordered$color, pch = 19, cex = .5)
UZ_scaffold12_ordered_N225 <- UZ_scaffold12[order(-UZ_scaffold12$N225_UZ),]
points(1:nrow(UZ_scaffold12_ordered_N225), UZ_scaffold12_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold12_ordered$color, pch = 19, cex = .5)
UZ_scaffold12_ordered_S25 <- UZ_scaffold12[order(-UZ_scaffold12$S25_UZ),]
points(1:nrow(UZ_scaffold12_ordered_S25), UZ_scaffold12_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold12_ordered$color, pch = 19, cex = .5)
UZ_scaffold12_ordered_S125 <- UZ_scaffold12[order(-UZ_scaffold12$S125_UZ),]
points(1:nrow(UZ_scaffold12_ordered_S125), UZ_scaffold12_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold12_ordered$color, pch = 19, cex = .5)
UZ_scaffold12_ordered_S175 <- UZ_scaffold12[order(-UZ_scaffold12$S175_UZ),]
points(1:nrow(UZ_scaffold12_ordered_S175), UZ_scaffold12_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold12_ordered$color, pch = 19, cex = .5)
UZ_scaffold12_ordered_S225 <- UZ_scaffold12[order(-UZ_scaffold12$S225_UZ),]
points(1:nrow(UZ_scaffold12_ordered_S225), UZ_scaffold12_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold12_ordered$color, pch = 19, cex = .5)


## scaffold13
FST_scaffold13 <- FST_m[which(FST_m$Scaffold == "scaffold13"), ]
FST_scaffold13 <- FST_scaffold13[sample(nrow(FST_scaffold13), 306), ]

UZ_scaffold13 <- FST_scaffold13[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold13 <- FST_scaffold13[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold13 <- mutate_all(UZ_scaffold13, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold13 <- mutate_all(Other_comparisons_scaffold13, function(x) as.numeric(as.character(x)))

UZ_scaffold13 <- cbind(UZ_scaffold13, FST_scaffold13$Window)
UZ_scaffold13$color <- "black"

UZ_scaffold13_ordered <- UZ_scaffold13[order(-UZ_scaffold13$N75_UZ),]
UZ_scaffold13$observation <- 1:nrow(UZ_scaffold13)
UZ_scaffold13_ordered$observation <- 1:nrow(UZ_scaffold13_ordered)
#plot(UZ_scaffold13_ordered$observation, UZ_scaffold13_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold13_ordered$color)
UZ_scaffold13_ordered_N125 <- UZ_scaffold13[order(-UZ_scaffold13$N125_UZ),]
points(1:nrow(UZ_scaffold13_ordered_N125), UZ_scaffold13_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold13_ordered$color, pch = 19, cex = .5)
UZ_scaffold13_ordered_N175 <- UZ_scaffold13[order(-UZ_scaffold13$N175_UZ),]
points(1:nrow(UZ_scaffold13_ordered_N175), UZ_scaffold13_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold13_ordered$color, pch = 19, cex = .5)
UZ_scaffold13_ordered_N225 <- UZ_scaffold13[order(-UZ_scaffold13$N225_UZ),]
points(1:nrow(UZ_scaffold13_ordered_N225), UZ_scaffold13_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold13_ordered$color, pch = 19, cex = .5)
UZ_scaffold13_ordered_S25 <- UZ_scaffold13[order(-UZ_scaffold13$S25_UZ),]
points(1:nrow(UZ_scaffold13_ordered_S25), UZ_scaffold13_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold13_ordered$color, pch = 19, cex = .5)
UZ_scaffold13_ordered_S125 <- UZ_scaffold13[order(-UZ_scaffold13$S125_UZ),]
points(1:nrow(UZ_scaffold13_ordered_S125), UZ_scaffold13_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold13_ordered$color, pch = 19, cex = .5)
UZ_scaffold13_ordered_S175 <- UZ_scaffold13[order(-UZ_scaffold13$S175_UZ),]
points(1:nrow(UZ_scaffold13_ordered_S175), UZ_scaffold13_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold13_ordered$color, pch = 19, cex = .5)
UZ_scaffold13_ordered_S225 <- UZ_scaffold13[order(-UZ_scaffold13$S225_UZ),]
points(1:nrow(UZ_scaffold13_ordered_S225), UZ_scaffold13_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold13_ordered$color, pch = 19, cex = .5)




## scaffold14
FST_scaffold14 <- FST_m[which(FST_m$Scaffold == "scaffold14"), ]
FST_scaffold14 <- FST_scaffold14[sample(nrow(FST_scaffold14), 306), ]

UZ_scaffold14 <- FST_scaffold14[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold14 <- FST_scaffold14[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold14 <- mutate_all(UZ_scaffold14, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold14 <- mutate_all(Other_comparisons_scaffold14, function(x) as.numeric(as.character(x)))

UZ_scaffold14 <- cbind(UZ_scaffold14, FST_scaffold14$Window)
UZ_scaffold14$color <- "black"

UZ_scaffold14_ordered <- UZ_scaffold14[order(-UZ_scaffold14$N75_UZ),]
UZ_scaffold14$observation <- 1:nrow(UZ_scaffold14)
UZ_scaffold14_ordered$observation <- 1:nrow(UZ_scaffold14_ordered)
#plot(UZ_scaffold14_ordered$observation, UZ_scaffold14_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold14_ordered$color)
UZ_scaffold14_ordered_N125 <- UZ_scaffold14[order(-UZ_scaffold14$N125_UZ),]
points(1:nrow(UZ_scaffold14_ordered_N125), UZ_scaffold14_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold14_ordered$color, pch = 19, cex = .5)
UZ_scaffold14_ordered_N175 <- UZ_scaffold14[order(-UZ_scaffold14$N175_UZ),]
points(1:nrow(UZ_scaffold14_ordered_N175), UZ_scaffold14_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold14_ordered$color, pch = 19, cex = .5)
UZ_scaffold14_ordered_N225 <- UZ_scaffold14[order(-UZ_scaffold14$N225_UZ),]
points(1:nrow(UZ_scaffold14_ordered_N225), UZ_scaffold14_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold14_ordered$color, pch = 19, cex = .5)
UZ_scaffold14_ordered_S25 <- UZ_scaffold14[order(-UZ_scaffold14$S25_UZ),]
points(1:nrow(UZ_scaffold14_ordered_S25), UZ_scaffold14_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold14_ordered$color, pch = 19, cex = .5)
UZ_scaffold14_ordered_S125 <- UZ_scaffold14[order(-UZ_scaffold14$S125_UZ),]
points(1:nrow(UZ_scaffold14_ordered_S125), UZ_scaffold14_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold14_ordered$color, pch = 19, cex = .5)
UZ_scaffold14_ordered_S175 <- UZ_scaffold14[order(-UZ_scaffold14$S175_UZ),]
points(1:nrow(UZ_scaffold14_ordered_S175), UZ_scaffold14_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold14_ordered$color, pch = 19, cex = .5)
UZ_scaffold14_ordered_S225 <- UZ_scaffold14[order(-UZ_scaffold14$S225_UZ),]
points(1:nrow(UZ_scaffold14_ordered_S225), UZ_scaffold14_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold14_ordered$color, pch = 19, cex = .5)


## scaffold15
FST_scaffold15 <- FST_m[which(FST_m$Scaffold == "scaffold15"), ]
FST_scaffold15 <- FST_scaffold15[sample(nrow(FST_scaffold15), 306), ]

UZ_scaffold15 <- FST_scaffold15[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold15 <- FST_scaffold15[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold15 <- mutate_all(UZ_scaffold15, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold15 <- mutate_all(Other_comparisons_scaffold15, function(x) as.numeric(as.character(x)))

UZ_scaffold15 <- cbind(UZ_scaffold15, FST_scaffold15$Window)
UZ_scaffold15$color <- "black"

UZ_scaffold15_ordered <- UZ_scaffold15[order(-UZ_scaffold15$N75_UZ),]
UZ_scaffold15$observation <- 1:nrow(UZ_scaffold15)
UZ_scaffold15_ordered$observation <- 1:nrow(UZ_scaffold15_ordered)
#plot(UZ_scaffold15_ordered$observation, UZ_scaffold15_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold15_ordered$color)
UZ_scaffold15_ordered_N125 <- UZ_scaffold15[order(-UZ_scaffold15$N125_UZ),]
points(1:nrow(UZ_scaffold15_ordered_N125), UZ_scaffold15_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold15_ordered$color, pch = 19, cex = .5)
UZ_scaffold15_ordered_N175 <- UZ_scaffold15[order(-UZ_scaffold15$N175_UZ),]
points(1:nrow(UZ_scaffold15_ordered_N175), UZ_scaffold15_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold15_ordered$color, pch = 19, cex = .5)
UZ_scaffold15_ordered_N225 <- UZ_scaffold15[order(-UZ_scaffold15$N225_UZ),]
points(1:nrow(UZ_scaffold15_ordered_N225), UZ_scaffold15_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold15_ordered$color, pch = 19, cex = .5)
UZ_scaffold15_ordered_S25 <- UZ_scaffold15[order(-UZ_scaffold15$S25_UZ),]
points(1:nrow(UZ_scaffold15_ordered_S25), UZ_scaffold15_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold15_ordered$color, pch = 19, cex = .5)
UZ_scaffold15_ordered_S125 <- UZ_scaffold15[order(-UZ_scaffold15$S125_UZ),]
points(1:nrow(UZ_scaffold15_ordered_S125), UZ_scaffold15_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold15_ordered$color, pch = 19, cex = .5)
UZ_scaffold15_ordered_S175 <- UZ_scaffold15[order(-UZ_scaffold15$S175_UZ),]
points(1:nrow(UZ_scaffold15_ordered_S175), UZ_scaffold15_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold15_ordered$color, pch = 19, cex = .5)
UZ_scaffold15_ordered_S225 <- UZ_scaffold15[order(-UZ_scaffold15$S225_UZ),]
points(1:nrow(UZ_scaffold15_ordered_S225), UZ_scaffold15_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold15_ordered$color, pch = 19, cex = .5)


## scaffold16
FST_scaffold16 <- FST_m[which(FST_m$Scaffold == "scaffold16"), ]
FST_scaffold16 <- FST_scaffold16[sample(nrow(FST_scaffold16), 306), ]

UZ_scaffold16 <- FST_scaffold16[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold16 <- FST_scaffold16[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold16 <- mutate_all(UZ_scaffold16, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold16 <- mutate_all(Other_comparisons_scaffold16, function(x) as.numeric(as.character(x)))

UZ_scaffold16 <- cbind(UZ_scaffold16, FST_scaffold16$Window)
UZ_scaffold16$color <- "black"

UZ_scaffold16_ordered <- UZ_scaffold16[order(-UZ_scaffold16$N75_UZ),]
UZ_scaffold16$observation <- 1:nrow(UZ_scaffold16)
UZ_scaffold16_ordered$observation <- 1:nrow(UZ_scaffold16_ordered)
#plot(UZ_scaffold16_ordered$observation, UZ_scaffold16_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold16_ordered$color)
UZ_scaffold16_ordered_N125 <- UZ_scaffold16[order(-UZ_scaffold16$N125_UZ),]
points(1:nrow(UZ_scaffold16_ordered_N125), UZ_scaffold16_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold16_ordered$color, pch = 19, cex = .5)
UZ_scaffold16_ordered_N175 <- UZ_scaffold16[order(-UZ_scaffold16$N175_UZ),]
points(1:nrow(UZ_scaffold16_ordered_N175), UZ_scaffold16_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold16_ordered$color, pch = 19, cex = .5)
UZ_scaffold16_ordered_N225 <- UZ_scaffold16[order(-UZ_scaffold16$N225_UZ),]
points(1:nrow(UZ_scaffold16_ordered_N225), UZ_scaffold16_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold16_ordered$color, pch = 19, cex = .5)
UZ_scaffold16_ordered_S25 <- UZ_scaffold16[order(-UZ_scaffold16$S25_UZ),]
points(1:nrow(UZ_scaffold16_ordered_S25), UZ_scaffold16_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold16_ordered$color, pch = 19, cex = .5)
UZ_scaffold16_ordered_S125 <- UZ_scaffold16[order(-UZ_scaffold16$S125_UZ),]
points(1:nrow(UZ_scaffold16_ordered_S125), UZ_scaffold16_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold16_ordered$color, pch = 19, cex = .5)
UZ_scaffold16_ordered_S175 <- UZ_scaffold16[order(-UZ_scaffold16$S175_UZ),]
points(1:nrow(UZ_scaffold16_ordered_S175), UZ_scaffold16_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold16_ordered$color, pch = 19, cex = .5)
UZ_scaffold16_ordered_S225 <- UZ_scaffold16[order(-UZ_scaffold16$S225_UZ),]
points(1:nrow(UZ_scaffold16_ordered_S225), UZ_scaffold16_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold16_ordered$color, pch = 19, cex = .5)



## scaffold17
FST_scaffold17 <- FST_m[which(FST_m$Scaffold == "scaffold17"), ]
FST_scaffold17 <- FST_scaffold17[sample(nrow(FST_scaffold17), 306), ]

UZ_scaffold17 <- FST_scaffold17[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold17 <- FST_scaffold17[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold17 <- mutate_all(UZ_scaffold17, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold17 <- mutate_all(Other_comparisons_scaffold17, function(x) as.numeric(as.character(x)))

UZ_scaffold17 <- cbind(UZ_scaffold17, FST_scaffold17$Window)
UZ_scaffold17$color <- "black"

UZ_scaffold17_ordered <- UZ_scaffold17[order(-UZ_scaffold17$N75_UZ),]
UZ_scaffold17$observation <- 1:nrow(UZ_scaffold17)
UZ_scaffold17_ordered$observation <- 1:nrow(UZ_scaffold17_ordered)
#plot(UZ_scaffold17_ordered$observation, UZ_scaffold17_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold17_ordered$color)
UZ_scaffold17_ordered_N125 <- UZ_scaffold17[order(-UZ_scaffold17$N125_UZ),]
points(1:nrow(UZ_scaffold17_ordered_N125), UZ_scaffold17_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold17_ordered$color, pch = 19, cex = .5)
UZ_scaffold17_ordered_N175 <- UZ_scaffold17[order(-UZ_scaffold17$N175_UZ),]
points(1:nrow(UZ_scaffold17_ordered_N175), UZ_scaffold17_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold17_ordered$color, pch = 19, cex = .5)
UZ_scaffold17_ordered_N225 <- UZ_scaffold17[order(-UZ_scaffold17$N225_UZ),]
points(1:nrow(UZ_scaffold17_ordered_N225), UZ_scaffold17_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold17_ordered$color, pch = 19, cex = .5)
UZ_scaffold17_ordered_S25 <- UZ_scaffold17[order(-UZ_scaffold17$S25_UZ),]
points(1:nrow(UZ_scaffold17_ordered_S25), UZ_scaffold17_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold17_ordered$color, pch = 19, cex = .5)
UZ_scaffold17_ordered_S125 <- UZ_scaffold17[order(-UZ_scaffold17$S125_UZ),]
points(1:nrow(UZ_scaffold17_ordered_S125), UZ_scaffold17_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold17_ordered$color, pch = 19, cex = .5)
UZ_scaffold17_ordered_S175 <- UZ_scaffold17[order(-UZ_scaffold17$S175_UZ),]
points(1:nrow(UZ_scaffold17_ordered_S175), UZ_scaffold17_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold17_ordered$color, pch = 19, cex = .5)
UZ_scaffold17_ordered_S225 <- UZ_scaffold17[order(-UZ_scaffold17$S225_UZ),]
points(1:nrow(UZ_scaffold17_ordered_S225), UZ_scaffold17_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold17_ordered$color, pch = 19, cex = .5)



## scaffold18
FST_scaffold18 <- FST_m[which(FST_m$Scaffold == "scaffold18"), ]
FST_scaffold18 <- FST_scaffold18[sample(nrow(FST_scaffold18), 306), ]

UZ_scaffold18 <- FST_scaffold18[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold18 <- FST_scaffold18[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold18 <- mutate_all(UZ_scaffold18, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold18 <- mutate_all(Other_comparisons_scaffold18, function(x) as.numeric(as.character(x)))

UZ_scaffold18 <- cbind(UZ_scaffold18, FST_scaffold18$Window)
UZ_scaffold18$color <- "black"

UZ_scaffold18_ordered <- UZ_scaffold18[order(-UZ_scaffold18$N75_UZ),]
UZ_scaffold18$observation <- 1:nrow(UZ_scaffold18)
UZ_scaffold18_ordered$observation <- 1:nrow(UZ_scaffold18_ordered)
#plot(UZ_scaffold18_ordered$observation, UZ_scaffold18_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold18_ordered$color)
UZ_scaffold18_ordered_N125 <- UZ_scaffold18[order(-UZ_scaffold18$N125_UZ),]
points(1:nrow(UZ_scaffold18_ordered_N125), UZ_scaffold18_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold18_ordered$color, pch = 19, cex = .5)
UZ_scaffold18_ordered_N175 <- UZ_scaffold18[order(-UZ_scaffold18$N175_UZ),]
points(1:nrow(UZ_scaffold18_ordered_N175), UZ_scaffold18_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold18_ordered$color, pch = 19, cex = .5)
UZ_scaffold18_ordered_N225 <- UZ_scaffold18[order(-UZ_scaffold18$N225_UZ),]
points(1:nrow(UZ_scaffold18_ordered_N225), UZ_scaffold18_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold18_ordered$color, pch = 19, cex = .5)
UZ_scaffold18_ordered_S25 <- UZ_scaffold18[order(-UZ_scaffold18$S25_UZ),]
points(1:nrow(UZ_scaffold18_ordered_S25), UZ_scaffold18_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold18_ordered$color, pch = 19, cex = .5)
UZ_scaffold18_ordered_S125 <- UZ_scaffold18[order(-UZ_scaffold18$S125_UZ),]
points(1:nrow(UZ_scaffold18_ordered_S125), UZ_scaffold18_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold18_ordered$color, pch = 19, cex = .5)
UZ_scaffold18_ordered_S175 <- UZ_scaffold18[order(-UZ_scaffold18$S175_UZ),]
points(1:nrow(UZ_scaffold18_ordered_S175), UZ_scaffold18_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold18_ordered$color, pch = 19, cex = .5)
UZ_scaffold18_ordered_S225 <- UZ_scaffold18[order(-UZ_scaffold18$S225_UZ),]
points(1:nrow(UZ_scaffold18_ordered_S225), UZ_scaffold18_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold18_ordered$color, pch = 19, cex = .5)



## scaffold19
FST_scaffold19 <- FST_m[which(FST_m$Scaffold == "scaffold19"), ]
FST_scaffold19 <- FST_scaffold19[sample(nrow(FST_scaffold19), 306), ]

UZ_scaffold19 <- FST_scaffold19[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold19 <- FST_scaffold19[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold19 <- mutate_all(UZ_scaffold19, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold19 <- mutate_all(Other_comparisons_scaffold19, function(x) as.numeric(as.character(x)))

UZ_scaffold19 <- cbind(UZ_scaffold19, FST_scaffold19$Window)
UZ_scaffold19$color <- "black"

UZ_scaffold19_ordered <- UZ_scaffold19[order(-UZ_scaffold19$N75_UZ),]
UZ_scaffold19$observation <- 1:nrow(UZ_scaffold19)
UZ_scaffold19_ordered$observation <- 1:nrow(UZ_scaffold19_ordered)
#plot(UZ_scaffold19_ordered$observation, UZ_scaffold19_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold19_ordered$color)
UZ_scaffold19_ordered_N125 <- UZ_scaffold19[order(-UZ_scaffold19$N125_UZ),]
points(1:nrow(UZ_scaffold19_ordered_N125), UZ_scaffold19_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold19_ordered$color, pch = 19, cex = .5)
UZ_scaffold19_ordered_N175 <- UZ_scaffold19[order(-UZ_scaffold19$N175_UZ),]
points(1:nrow(UZ_scaffold19_ordered_N175), UZ_scaffold19_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold19_ordered$color, pch = 19, cex = .5)
UZ_scaffold19_ordered_N225 <- UZ_scaffold19[order(-UZ_scaffold19$N225_UZ),]
points(1:nrow(UZ_scaffold19_ordered_N225), UZ_scaffold19_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold19_ordered$color, pch = 19, cex = .5)
UZ_scaffold19_ordered_S25 <- UZ_scaffold19[order(-UZ_scaffold19$S25_UZ),]
points(1:nrow(UZ_scaffold19_ordered_S25), UZ_scaffold19_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold19_ordered$color, pch = 19, cex = .5)
UZ_scaffold19_ordered_S125 <- UZ_scaffold19[order(-UZ_scaffold19$S125_UZ),]
points(1:nrow(UZ_scaffold19_ordered_S125), UZ_scaffold19_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold19_ordered$color, pch = 19, cex = .5)
UZ_scaffold19_ordered_S175 <- UZ_scaffold19[order(-UZ_scaffold19$S175_UZ),]
points(1:nrow(UZ_scaffold19_ordered_S175), UZ_scaffold19_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold19_ordered$color, pch = 19, cex = .5)
UZ_scaffold19_ordered_S225 <- UZ_scaffold19[order(-UZ_scaffold19$S225_UZ),]
points(1:nrow(UZ_scaffold19_ordered_S225), UZ_scaffold19_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold19_ordered$color, pch = 19, cex = .5)



## scaffold110
FST_scaffold110 <- FST_m[which(FST_m$Scaffold == "scaffold110"), ]
FST_scaffold110 <- FST_scaffold110[sample(nrow(FST_scaffold110), 306), ]

UZ_scaffold110 <- FST_scaffold110[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_scaffold110 <- FST_scaffold110[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_scaffold110 <- mutate_all(UZ_scaffold110, function(x) as.numeric(as.character(x)))
Other_comparisons_scaffold110 <- mutate_all(Other_comparisons_scaffold110, function(x) as.numeric(as.character(x)))

UZ_scaffold110 <- cbind(UZ_scaffold110, FST_scaffold110$Window)
UZ_scaffold110$color <- "black"

UZ_scaffold110_ordered <- UZ_scaffold110[order(-UZ_scaffold110$N75_UZ),]
UZ_scaffold110$observation <- 1:nrow(UZ_scaffold110)
UZ_scaffold110_ordered$observation <- 1:nrow(UZ_scaffold110_ordered)
#plot(UZ_scaffold110_ordered$observation, UZ_scaffold110_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_scaffold110_ordered$color)
UZ_scaffold110_ordered_N125 <- UZ_scaffold110[order(-UZ_scaffold110$N125_UZ),]
points(1:nrow(UZ_scaffold110_ordered_N125), UZ_scaffold110_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_scaffold110_ordered$color, pch = 19, cex = .5)
UZ_scaffold110_ordered_N175 <- UZ_scaffold110[order(-UZ_scaffold110$N175_UZ),]
points(1:nrow(UZ_scaffold110_ordered_N175), UZ_scaffold110_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_scaffold110_ordered$color, pch = 19, cex = .5)
UZ_scaffold110_ordered_N225 <- UZ_scaffold110[order(-UZ_scaffold110$N225_UZ),]
points(1:nrow(UZ_scaffold110_ordered_N225), UZ_scaffold110_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_scaffold110_ordered$color, pch = 19, cex = .5)
UZ_scaffold110_ordered_S25 <- UZ_scaffold110[order(-UZ_scaffold110$S25_UZ),]
points(1:nrow(UZ_scaffold110_ordered_S25), UZ_scaffold110_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_scaffold110_ordered$color, pch = 19, cex = .5)
UZ_scaffold110_ordered_S125 <- UZ_scaffold110[order(-UZ_scaffold110$S125_UZ),]
points(1:nrow(UZ_scaffold110_ordered_S125), UZ_scaffold110_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_scaffold110_ordered$color, pch = 19, cex = .5)
UZ_scaffold110_ordered_S175 <- UZ_scaffold110[order(-UZ_scaffold110$S175_UZ),]
points(1:nrow(UZ_scaffold110_ordered_S175), UZ_scaffold110_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_scaffold110_ordered$color, pch = 19, cex = .5)
UZ_scaffold110_ordered_S225 <- UZ_scaffold110[order(-UZ_scaffold110$S225_UZ),]
points(1:nrow(UZ_scaffold110_ordered_S225), UZ_scaffold110_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_scaffold110_ordered$color, pch = 19, cex = .5)



## Scaffold20
FST_Scaffold20 <- FST_m[which(FST_m$Scaffold == "scaffold20"), ]
FST_Scaffold20 <- FST_Scaffold20[sample(nrow(FST_Scaffold20), 306), ]

UZ_Scaffold20 <- FST_Scaffold20[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold20 <- FST_Scaffold20[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold20 <- mutate_all(UZ_Scaffold20, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold20 <- mutate_all(Other_comparisons_Scaffold20, function(x) as.numeric(as.character(x)))

UZ_Scaffold20 <- cbind(UZ_Scaffold20, FST_Scaffold20$Window)
UZ_Scaffold20$color <- "black"

UZ_Scaffold20_ordered <- UZ_Scaffold20[order(-UZ_Scaffold20$N75_UZ),]
UZ_Scaffold20$observation <- 1:nrow(UZ_Scaffold20)
UZ_Scaffold20_ordered$observation <- 1:nrow(UZ_Scaffold20_ordered)
#plot(UZ_Scaffold20_ordered$observation, UZ_Scaffold20_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold20_ordered$color)
UZ_Scaffold20_ordered_N125 <- UZ_Scaffold20[order(-UZ_Scaffold20$N125_UZ),]
points(1:nrow(UZ_Scaffold20_ordered_N125), UZ_Scaffold20_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold20_ordered$color, pch = 19, cex = .5)
UZ_Scaffold20_ordered_N175 <- UZ_Scaffold20[order(-UZ_Scaffold20$N175_UZ),]
points(1:nrow(UZ_Scaffold20_ordered_N175), UZ_Scaffold20_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold20_ordered$color, pch = 19, cex = .5)
UZ_Scaffold20_ordered_N225 <- UZ_Scaffold20[order(-UZ_Scaffold20$N225_UZ),]
points(1:nrow(UZ_Scaffold20_ordered_N225), UZ_Scaffold20_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold20_ordered$color, pch = 19, cex = .5)
UZ_Scaffold20_ordered_S25 <- UZ_Scaffold20[order(-UZ_Scaffold20$S25_UZ),]
points(1:nrow(UZ_Scaffold20_ordered_S25), UZ_Scaffold20_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold20_ordered$color, pch = 19, cex = .5)
UZ_Scaffold20_ordered_S125 <- UZ_Scaffold20[order(-UZ_Scaffold20$S125_UZ),]
points(1:nrow(UZ_Scaffold20_ordered_S125), UZ_Scaffold20_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold20_ordered$color, pch = 19, cex = .5)
UZ_Scaffold20_ordered_S175 <- UZ_Scaffold20[order(-UZ_Scaffold20$S175_UZ),]
points(1:nrow(UZ_Scaffold20_ordered_S175), UZ_Scaffold20_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold20_ordered$color, pch = 19, cex = .5)
UZ_Scaffold20_ordered_S225 <- UZ_Scaffold20[order(-UZ_Scaffold20$S225_UZ),]
points(1:nrow(UZ_Scaffold20_ordered_S225), UZ_Scaffold20_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold20_ordered$color, pch = 19, cex = .5)

## Scaffold21
FST_Scaffold21 <- FST_m[which(FST_m$Scaffold == "scaffold21"), ]
FST_Scaffold21 <- FST_Scaffold21[sample(nrow(FST_Scaffold21), 306), ]

UZ_Scaffold21 <- FST_Scaffold21[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold21 <- FST_Scaffold21[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold21 <- mutate_all(UZ_Scaffold21, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold21 <- mutate_all(Other_comparisons_Scaffold21, function(x) as.numeric(as.character(x)))

UZ_Scaffold21 <- cbind(UZ_Scaffold21, FST_Scaffold21$Window)
UZ_Scaffold21$color <- "black"

UZ_Scaffold21_ordered <- UZ_Scaffold21[order(-UZ_Scaffold21$N75_UZ),]
UZ_Scaffold21$observation <- 1:nrow(UZ_Scaffold21)
UZ_Scaffold21_ordered$observation <- 1:nrow(UZ_Scaffold21_ordered)
#plot(UZ_Scaffold21_ordered$observation, UZ_Scaffold21_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold21_ordered$color)
UZ_Scaffold21_ordered_N125 <- UZ_Scaffold21[order(-UZ_Scaffold21$N125_UZ),]
points(1:nrow(UZ_Scaffold21_ordered_N125), UZ_Scaffold21_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold21_ordered$color, pch = 19, cex = .5)
UZ_Scaffold21_ordered_N175 <- UZ_Scaffold21[order(-UZ_Scaffold21$N175_UZ),]
points(1:nrow(UZ_Scaffold21_ordered_N175), UZ_Scaffold21_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold21_ordered$color, pch = 19, cex = .5)
UZ_Scaffold21_ordered_N225 <- UZ_Scaffold21[order(-UZ_Scaffold21$N225_UZ),]
points(1:nrow(UZ_Scaffold21_ordered_N225), UZ_Scaffold21_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold21_ordered$color, pch = 19, cex = .5)
UZ_Scaffold21_ordered_S25 <- UZ_Scaffold21[order(-UZ_Scaffold21$S25_UZ),]
points(1:nrow(UZ_Scaffold21_ordered_S25), UZ_Scaffold21_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold21_ordered$color, pch = 19, cex = .5)
UZ_Scaffold21_ordered_S125 <- UZ_Scaffold21[order(-UZ_Scaffold21$S125_UZ),]
points(1:nrow(UZ_Scaffold21_ordered_S125), UZ_Scaffold21_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold21_ordered$color, pch = 19, cex = .5)
UZ_Scaffold21_ordered_S175 <- UZ_Scaffold21[order(-UZ_Scaffold21$S175_UZ),]
points(1:nrow(UZ_Scaffold21_ordered_S175), UZ_Scaffold21_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold21_ordered$color, pch = 19, cex = .5)
UZ_Scaffold21_ordered_S225 <- UZ_Scaffold21[order(-UZ_Scaffold21$S225_UZ),]
points(1:nrow(UZ_Scaffold21_ordered_S225), UZ_Scaffold21_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold21_ordered$color, pch = 19, cex = .5)

## Scaffold22
FST_Scaffold22 <- FST_m[which(FST_m$Scaffold == "scaffold22"), ]
FST_Scaffold22 <- FST_Scaffold22[sample(nrow(FST_Scaffold22), 306), ]

UZ_Scaffold22 <- FST_Scaffold22[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold22 <- FST_Scaffold22[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold22 <- mutate_all(UZ_Scaffold22, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold22 <- mutate_all(Other_comparisons_Scaffold22, function(x) as.numeric(as.character(x)))

UZ_Scaffold22 <- cbind(UZ_Scaffold22, FST_Scaffold22$Window)
UZ_Scaffold22$color <- "black"

UZ_Scaffold22_ordered <- UZ_Scaffold22[order(-UZ_Scaffold22$N75_UZ),]
UZ_Scaffold22$observation <- 1:nrow(UZ_Scaffold22)
UZ_Scaffold22_ordered$observation <- 1:nrow(UZ_Scaffold22_ordered)
#plot(UZ_Scaffold22_ordered$observation, UZ_Scaffold22_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold22_ordered$color)
UZ_Scaffold22_ordered_N125 <- UZ_Scaffold22[order(-UZ_Scaffold22$N125_UZ),]
points(1:nrow(UZ_Scaffold22_ordered_N125), UZ_Scaffold22_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold22_ordered$color, pch = 19, cex = .5)
UZ_Scaffold22_ordered_N175 <- UZ_Scaffold22[order(-UZ_Scaffold22$N175_UZ),]
points(1:nrow(UZ_Scaffold22_ordered_N175), UZ_Scaffold22_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold22_ordered$color, pch = 19, cex = .5)
UZ_Scaffold22_ordered_N225 <- UZ_Scaffold22[order(-UZ_Scaffold22$N225_UZ),]
points(1:nrow(UZ_Scaffold22_ordered_N225), UZ_Scaffold22_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold22_ordered$color, pch = 19, cex = .5)
UZ_Scaffold22_ordered_S25 <- UZ_Scaffold22[order(-UZ_Scaffold22$S25_UZ),]
points(1:nrow(UZ_Scaffold22_ordered_S25), UZ_Scaffold22_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold22_ordered$color, pch = 19, cex = .5)
UZ_Scaffold22_ordered_S125 <- UZ_Scaffold22[order(-UZ_Scaffold22$S125_UZ),]
points(1:nrow(UZ_Scaffold22_ordered_S125), UZ_Scaffold22_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold22_ordered$color, pch = 19, cex = .5)
UZ_Scaffold22_ordered_S175 <- UZ_Scaffold22[order(-UZ_Scaffold22$S175_UZ),]
points(1:nrow(UZ_Scaffold22_ordered_S175), UZ_Scaffold22_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold22_ordered$color, pch = 19, cex = .5)
UZ_Scaffold22_ordered_S225 <- UZ_Scaffold22[order(-UZ_Scaffold22$S225_UZ),]
points(1:nrow(UZ_Scaffold22_ordered_S225), UZ_Scaffold22_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold22_ordered$color, pch = 19, cex = .5)


## Scaffold23
FST_Scaffold23 <- FST_m[which(FST_m$Scaffold == "scaffold23"), ]
FST_Scaffold23 <- FST_Scaffold23[sample(nrow(FST_Scaffold23), 306), ]

UZ_Scaffold23 <- FST_Scaffold23[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold23 <- FST_Scaffold23[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold23 <- mutate_all(UZ_Scaffold23, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold23 <- mutate_all(Other_comparisons_Scaffold23, function(x) as.numeric(as.character(x)))

UZ_Scaffold23 <- cbind(UZ_Scaffold23, FST_Scaffold23$Window)
UZ_Scaffold23$color <- "black"

UZ_Scaffold23_ordered <- UZ_Scaffold23[order(-UZ_Scaffold23$N75_UZ),]
UZ_Scaffold23$observation <- 1:nrow(UZ_Scaffold23)
UZ_Scaffold23_ordered$observation <- 1:nrow(UZ_Scaffold23_ordered)
#plot(UZ_Scaffold23_ordered$observation, UZ_Scaffold23_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold23_ordered$color)
UZ_Scaffold23_ordered_N125 <- UZ_Scaffold23[order(-UZ_Scaffold23$N125_UZ),]
points(1:nrow(UZ_Scaffold23_ordered_N125), UZ_Scaffold23_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold23_ordered$color, pch = 19, cex = .5)
UZ_Scaffold23_ordered_N175 <- UZ_Scaffold23[order(-UZ_Scaffold23$N175_UZ),]
points(1:nrow(UZ_Scaffold23_ordered_N175), UZ_Scaffold23_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold23_ordered$color, pch = 19, cex = .5)
UZ_Scaffold23_ordered_N225 <- UZ_Scaffold23[order(-UZ_Scaffold23$N225_UZ),]
points(1:nrow(UZ_Scaffold23_ordered_N225), UZ_Scaffold23_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold23_ordered$color, pch = 19, cex = .5)
UZ_Scaffold23_ordered_S25 <- UZ_Scaffold23[order(-UZ_Scaffold23$S25_UZ),]
points(1:nrow(UZ_Scaffold23_ordered_S25), UZ_Scaffold23_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold23_ordered$color, pch = 19, cex = .5)
UZ_Scaffold23_ordered_S125 <- UZ_Scaffold23[order(-UZ_Scaffold23$S125_UZ),]
points(1:nrow(UZ_Scaffold23_ordered_S125), UZ_Scaffold23_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold23_ordered$color, pch = 19, cex = .5)
UZ_Scaffold23_ordered_S175 <- UZ_Scaffold23[order(-UZ_Scaffold23$S175_UZ),]
points(1:nrow(UZ_Scaffold23_ordered_S175), UZ_Scaffold23_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold23_ordered$color, pch = 19, cex = .5)
UZ_Scaffold23_ordered_S225 <- UZ_Scaffold23[order(-UZ_Scaffold23$S225_UZ),]
points(1:nrow(UZ_Scaffold23_ordered_S225), UZ_Scaffold23_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold23_ordered$color, pch = 19, cex = .5)




## Scaffold24
FST_Scaffold24 <- FST_m[which(FST_m$Scaffold == "scaffold24"), ]
FST_Scaffold24 <- FST_Scaffold24[sample(nrow(FST_Scaffold24), 306), ]

UZ_Scaffold24 <- FST_Scaffold24[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold24 <- FST_Scaffold24[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold24 <- mutate_all(UZ_Scaffold24, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold24 <- mutate_all(Other_comparisons_Scaffold24, function(x) as.numeric(as.character(x)))

UZ_Scaffold24 <- cbind(UZ_Scaffold24, FST_Scaffold24$Window)
UZ_Scaffold24$color <- "black"

UZ_Scaffold24_ordered <- UZ_Scaffold24[order(-UZ_Scaffold24$N75_UZ),]
UZ_Scaffold24$observation <- 1:nrow(UZ_Scaffold24)
UZ_Scaffold24_ordered$observation <- 1:nrow(UZ_Scaffold24_ordered)
#plot(UZ_Scaffold24_ordered$observation, UZ_Scaffold24_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold24_ordered$color)
UZ_Scaffold24_ordered_N125 <- UZ_Scaffold24[order(-UZ_Scaffold24$N125_UZ),]
points(1:nrow(UZ_Scaffold24_ordered_N125), UZ_Scaffold24_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold24_ordered$color, pch = 19, cex = .5)
UZ_Scaffold24_ordered_N175 <- UZ_Scaffold24[order(-UZ_Scaffold24$N175_UZ),]
points(1:nrow(UZ_Scaffold24_ordered_N175), UZ_Scaffold24_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold24_ordered$color, pch = 19, cex = .5)
UZ_Scaffold24_ordered_N225 <- UZ_Scaffold24[order(-UZ_Scaffold24$N225_UZ),]
points(1:nrow(UZ_Scaffold24_ordered_N225), UZ_Scaffold24_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold24_ordered$color, pch = 19, cex = .5)
UZ_Scaffold24_ordered_S25 <- UZ_Scaffold24[order(-UZ_Scaffold24$S25_UZ),]
points(1:nrow(UZ_Scaffold24_ordered_S25), UZ_Scaffold24_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold24_ordered$color, pch = 19, cex = .5)
UZ_Scaffold24_ordered_S125 <- UZ_Scaffold24[order(-UZ_Scaffold24$S125_UZ),]
points(1:nrow(UZ_Scaffold24_ordered_S125), UZ_Scaffold24_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold24_ordered$color, pch = 19, cex = .5)
UZ_Scaffold24_ordered_S175 <- UZ_Scaffold24[order(-UZ_Scaffold24$S175_UZ),]
points(1:nrow(UZ_Scaffold24_ordered_S175), UZ_Scaffold24_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold24_ordered$color, pch = 19, cex = .5)
UZ_Scaffold24_ordered_S225 <- UZ_Scaffold24[order(-UZ_Scaffold24$S225_UZ),]
points(1:nrow(UZ_Scaffold24_ordered_S225), UZ_Scaffold24_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold24_ordered$color, pch = 19, cex = .5)


## Scaffold25
FST_Scaffold25 <- FST_m[which(FST_m$Scaffold == "scaffold25"), ]
FST_Scaffold25 <- FST_Scaffold25[sample(nrow(FST_Scaffold25), 306), ]

UZ_Scaffold25 <- FST_Scaffold25[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold25 <- FST_Scaffold25[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold25 <- mutate_all(UZ_Scaffold25, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold25 <- mutate_all(Other_comparisons_Scaffold25, function(x) as.numeric(as.character(x)))

UZ_Scaffold25 <- cbind(UZ_Scaffold25, FST_Scaffold25$Window)
UZ_Scaffold25$color <- "black"

UZ_Scaffold25_ordered <- UZ_Scaffold25[order(-UZ_Scaffold25$N75_UZ),]
UZ_Scaffold25$observation <- 1:nrow(UZ_Scaffold25)
UZ_Scaffold25_ordered$observation <- 1:nrow(UZ_Scaffold25_ordered)
#plot(UZ_Scaffold25_ordered$observation, UZ_Scaffold25_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold25_ordered$color)
UZ_Scaffold25_ordered_N125 <- UZ_Scaffold25[order(-UZ_Scaffold25$N125_UZ),]
points(1:nrow(UZ_Scaffold25_ordered_N125), UZ_Scaffold25_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold25_ordered$color, pch = 19, cex = .5)
UZ_Scaffold25_ordered_N175 <- UZ_Scaffold25[order(-UZ_Scaffold25$N175_UZ),]
points(1:nrow(UZ_Scaffold25_ordered_N175), UZ_Scaffold25_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold25_ordered$color, pch = 19, cex = .5)
UZ_Scaffold25_ordered_N225 <- UZ_Scaffold25[order(-UZ_Scaffold25$N225_UZ),]
points(1:nrow(UZ_Scaffold25_ordered_N225), UZ_Scaffold25_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold25_ordered$color, pch = 19, cex = .5)
UZ_Scaffold25_ordered_S25 <- UZ_Scaffold25[order(-UZ_Scaffold25$S25_UZ),]
points(1:nrow(UZ_Scaffold25_ordered_S25), UZ_Scaffold25_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold25_ordered$color, pch = 19, cex = .5)
UZ_Scaffold25_ordered_S125 <- UZ_Scaffold25[order(-UZ_Scaffold25$S125_UZ),]
points(1:nrow(UZ_Scaffold25_ordered_S125), UZ_Scaffold25_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold25_ordered$color, pch = 19, cex = .5)
UZ_Scaffold25_ordered_S175 <- UZ_Scaffold25[order(-UZ_Scaffold25$S175_UZ),]
points(1:nrow(UZ_Scaffold25_ordered_S175), UZ_Scaffold25_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold25_ordered$color, pch = 19, cex = .5)
UZ_Scaffold25_ordered_S225 <- UZ_Scaffold25[order(-UZ_Scaffold25$S225_UZ),]
points(1:nrow(UZ_Scaffold25_ordered_S225), UZ_Scaffold25_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold25_ordered$color, pch = 19, cex = .5)


## Scaffold26
FST_Scaffold26 <- FST_m[which(FST_m$Scaffold == "scaffold26"), ]
FST_Scaffold26 <- FST_Scaffold26[sample(nrow(FST_Scaffold26), 306), ]

UZ_Scaffold26 <- FST_Scaffold26[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold26 <- FST_Scaffold26[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold26 <- mutate_all(UZ_Scaffold26, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold26 <- mutate_all(Other_comparisons_Scaffold26, function(x) as.numeric(as.character(x)))

UZ_Scaffold26 <- cbind(UZ_Scaffold26, FST_Scaffold26$Window)
UZ_Scaffold26$color <- "black"

UZ_Scaffold26_ordered <- UZ_Scaffold26[order(-UZ_Scaffold26$N75_UZ),]
UZ_Scaffold26$observation <- 1:nrow(UZ_Scaffold26)
UZ_Scaffold26_ordered$observation <- 1:nrow(UZ_Scaffold26_ordered)
#plot(UZ_Scaffold26_ordered$observation, UZ_Scaffold26_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold26_ordered$color)
UZ_Scaffold26_ordered_N125 <- UZ_Scaffold26[order(-UZ_Scaffold26$N125_UZ),]
points(1:nrow(UZ_Scaffold26_ordered_N125), UZ_Scaffold26_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold26_ordered$color, pch = 19, cex = .5)
UZ_Scaffold26_ordered_N175 <- UZ_Scaffold26[order(-UZ_Scaffold26$N175_UZ),]
points(1:nrow(UZ_Scaffold26_ordered_N175), UZ_Scaffold26_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold26_ordered$color, pch = 19, cex = .5)
UZ_Scaffold26_ordered_N225 <- UZ_Scaffold26[order(-UZ_Scaffold26$N225_UZ),]
points(1:nrow(UZ_Scaffold26_ordered_N225), UZ_Scaffold26_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold26_ordered$color, pch = 19, cex = .5)
UZ_Scaffold26_ordered_S25 <- UZ_Scaffold26[order(-UZ_Scaffold26$S25_UZ),]
points(1:nrow(UZ_Scaffold26_ordered_S25), UZ_Scaffold26_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold26_ordered$color, pch = 19, cex = .5)
UZ_Scaffold26_ordered_S125 <- UZ_Scaffold26[order(-UZ_Scaffold26$S125_UZ),]
points(1:nrow(UZ_Scaffold26_ordered_S125), UZ_Scaffold26_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold26_ordered$color, pch = 19, cex = .5)
UZ_Scaffold26_ordered_S175 <- UZ_Scaffold26[order(-UZ_Scaffold26$S175_UZ),]
points(1:nrow(UZ_Scaffold26_ordered_S175), UZ_Scaffold26_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold26_ordered$color, pch = 19, cex = .5)
UZ_Scaffold26_ordered_S225 <- UZ_Scaffold26[order(-UZ_Scaffold26$S225_UZ),]
points(1:nrow(UZ_Scaffold26_ordered_S225), UZ_Scaffold26_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold26_ordered$color, pch = 19, cex = .5)



## Scaffold27
FST_Scaffold27 <- FST_m[which(FST_m$Scaffold == "scaffold27"), ]
FST_Scaffold27 <- FST_Scaffold27[sample(nrow(FST_Scaffold27), 306), ]

UZ_Scaffold27 <- FST_Scaffold27[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold27 <- FST_Scaffold27[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold27 <- mutate_all(UZ_Scaffold27, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold27 <- mutate_all(Other_comparisons_Scaffold27, function(x) as.numeric(as.character(x)))

UZ_Scaffold27 <- cbind(UZ_Scaffold27, FST_Scaffold27$Window)
UZ_Scaffold27$color <- "black"

UZ_Scaffold27_ordered <- UZ_Scaffold27[order(-UZ_Scaffold27$N75_UZ),]
UZ_Scaffold27$observation <- 1:nrow(UZ_Scaffold27)
UZ_Scaffold27_ordered$observation <- 1:nrow(UZ_Scaffold27_ordered)
#plot(UZ_Scaffold27_ordered$observation, UZ_Scaffold27_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold27_ordered$color)
UZ_Scaffold27_ordered_N125 <- UZ_Scaffold27[order(-UZ_Scaffold27$N125_UZ),]
points(1:nrow(UZ_Scaffold27_ordered_N125), UZ_Scaffold27_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold27_ordered$color, pch = 19, cex = .5)
UZ_Scaffold27_ordered_N175 <- UZ_Scaffold27[order(-UZ_Scaffold27$N175_UZ),]
points(1:nrow(UZ_Scaffold27_ordered_N175), UZ_Scaffold27_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold27_ordered$color, pch = 19, cex = .5)
UZ_Scaffold27_ordered_N225 <- UZ_Scaffold27[order(-UZ_Scaffold27$N225_UZ),]
points(1:nrow(UZ_Scaffold27_ordered_N225), UZ_Scaffold27_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold27_ordered$color, pch = 19, cex = .5)
UZ_Scaffold27_ordered_S25 <- UZ_Scaffold27[order(-UZ_Scaffold27$S25_UZ),]
points(1:nrow(UZ_Scaffold27_ordered_S25), UZ_Scaffold27_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold27_ordered$color, pch = 19, cex = .5)
UZ_Scaffold27_ordered_S125 <- UZ_Scaffold27[order(-UZ_Scaffold27$S125_UZ),]
points(1:nrow(UZ_Scaffold27_ordered_S125), UZ_Scaffold27_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold27_ordered$color, pch = 19, cex = .5)
UZ_Scaffold27_ordered_S175 <- UZ_Scaffold27[order(-UZ_Scaffold27$S175_UZ),]
points(1:nrow(UZ_Scaffold27_ordered_S175), UZ_Scaffold27_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold27_ordered$color, pch = 19, cex = .5)
UZ_Scaffold27_ordered_S225 <- UZ_Scaffold27[order(-UZ_Scaffold27$S225_UZ),]
points(1:nrow(UZ_Scaffold27_ordered_S225), UZ_Scaffold27_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold27_ordered$color, pch = 19, cex = .5)



## Scaffold28
FST_Scaffold28 <- FST_m[which(FST_m$Scaffold == "scaffold28"), ]
FST_Scaffold28 <- FST_Scaffold28[sample(nrow(FST_Scaffold28), 306), ]

UZ_Scaffold28 <- FST_Scaffold28[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold28 <- FST_Scaffold28[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold28 <- mutate_all(UZ_Scaffold28, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold28 <- mutate_all(Other_comparisons_Scaffold28, function(x) as.numeric(as.character(x)))

UZ_Scaffold28 <- cbind(UZ_Scaffold28, FST_Scaffold28$Window)
UZ_Scaffold28$color <- "black"

UZ_Scaffold28_ordered <- UZ_Scaffold28[order(-UZ_Scaffold28$N75_UZ),]
UZ_Scaffold28$observation <- 1:nrow(UZ_Scaffold28)
UZ_Scaffold28_ordered$observation <- 1:nrow(UZ_Scaffold28_ordered)
#plot(UZ_Scaffold28_ordered$observation, UZ_Scaffold28_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold28_ordered$color)
UZ_Scaffold28_ordered_N125 <- UZ_Scaffold28[order(-UZ_Scaffold28$N125_UZ),]
points(1:nrow(UZ_Scaffold28_ordered_N125), UZ_Scaffold28_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold28_ordered$color, pch = 19, cex = .5)
UZ_Scaffold28_ordered_N175 <- UZ_Scaffold28[order(-UZ_Scaffold28$N175_UZ),]
points(1:nrow(UZ_Scaffold28_ordered_N175), UZ_Scaffold28_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold28_ordered$color, pch = 19, cex = .5)
UZ_Scaffold28_ordered_N225 <- UZ_Scaffold28[order(-UZ_Scaffold28$N225_UZ),]
points(1:nrow(UZ_Scaffold28_ordered_N225), UZ_Scaffold28_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold28_ordered$color, pch = 19, cex = .5)
UZ_Scaffold28_ordered_S25 <- UZ_Scaffold28[order(-UZ_Scaffold28$S25_UZ),]
points(1:nrow(UZ_Scaffold28_ordered_S25), UZ_Scaffold28_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold28_ordered$color, pch = 19, cex = .5)
UZ_Scaffold28_ordered_S125 <- UZ_Scaffold28[order(-UZ_Scaffold28$S125_UZ),]
points(1:nrow(UZ_Scaffold28_ordered_S125), UZ_Scaffold28_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold28_ordered$color, pch = 19, cex = .5)
UZ_Scaffold28_ordered_S175 <- UZ_Scaffold28[order(-UZ_Scaffold28$S175_UZ),]
points(1:nrow(UZ_Scaffold28_ordered_S175), UZ_Scaffold28_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold28_ordered$color, pch = 19, cex = .5)
UZ_Scaffold28_ordered_S225 <- UZ_Scaffold28[order(-UZ_Scaffold28$S225_UZ),]
points(1:nrow(UZ_Scaffold28_ordered_S225), UZ_Scaffold28_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold28_ordered$color, pch = 19, cex = .5)



## Scaffold29
FST_Scaffold29 <- FST_m[which(FST_m$Scaffold == "scaffold29"), ]
FST_Scaffold29 <- FST_Scaffold29[sample(nrow(FST_Scaffold29), 306), ]

UZ_Scaffold29 <- FST_Scaffold29[, c(25, 33, 40, 46, 51, 58, 60, 61)]
Other_comparisons_Scaffold29 <- FST_Scaffold29[, c(17:20, 22:24, 26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59)]

UZ_Scaffold29 <- mutate_all(UZ_Scaffold29, function(x) as.numeric(as.character(x)))
Other_comparisons_Scaffold29 <- mutate_all(Other_comparisons_Scaffold29, function(x) as.numeric(as.character(x)))

UZ_Scaffold29 <- cbind(UZ_Scaffold29, FST_Scaffold29$Window)
UZ_Scaffold29$color <- "black"

UZ_Scaffold29_ordered <- UZ_Scaffold29[order(-UZ_Scaffold29$N75_UZ),]
UZ_Scaffold29$observation <- 1:nrow(UZ_Scaffold29)
UZ_Scaffold29_ordered$observation <- 1:nrow(UZ_Scaffold29_ordered)
#plot(UZ_Scaffold29_ordered$observation, UZ_Scaffold29_ordered$N75_UZ, ylim = c(0, 1), cex = 0.5, pch = 19, col= UZ_Scaffold29_ordered$color)
UZ_Scaffold29_ordered_N125 <- UZ_Scaffold29[order(-UZ_Scaffold29$N125_UZ),]
points(1:nrow(UZ_Scaffold29_ordered_N125), UZ_Scaffold29_ordered_N125$N125_UZ, ylim = c(0, 1), col= UZ_Scaffold29_ordered$color, pch = 19, cex = .5)
UZ_Scaffold29_ordered_N175 <- UZ_Scaffold29[order(-UZ_Scaffold29$N175_UZ),]
points(1:nrow(UZ_Scaffold29_ordered_N175), UZ_Scaffold29_ordered_N175$N175_UZ, ylim = c(0, 1), col= UZ_Scaffold29_ordered$color, pch = 19, cex = .5)
UZ_Scaffold29_ordered_N225 <- UZ_Scaffold29[order(-UZ_Scaffold29$N225_UZ),]
points(1:nrow(UZ_Scaffold29_ordered_N225), UZ_Scaffold29_ordered_N225$N225_UZ, ylim = c(0, 1), col= UZ_Scaffold29_ordered$color, pch = 19, cex = .5)
UZ_Scaffold29_ordered_S25 <- UZ_Scaffold29[order(-UZ_Scaffold29$S25_UZ),]
points(1:nrow(UZ_Scaffold29_ordered_S25), UZ_Scaffold29_ordered_S25$S25_UZ, ylim = c(0, 1), col= UZ_Scaffold29_ordered$color, pch = 19, cex = .5)
UZ_Scaffold29_ordered_S125 <- UZ_Scaffold29[order(-UZ_Scaffold29$S125_UZ),]
points(1:nrow(UZ_Scaffold29_ordered_S125), UZ_Scaffold29_ordered_S125$S125_UZ, ylim = c(0, 1), col= UZ_Scaffold29_ordered$color, pch = 19, cex = .5)
UZ_Scaffold29_ordered_S175 <- UZ_Scaffold29[order(-UZ_Scaffold29$S175_UZ),]
points(1:nrow(UZ_Scaffold29_ordered_S175), UZ_Scaffold29_ordered_S175$S175_UZ, ylim = c(0, 1), col= UZ_Scaffold29_ordered$color, pch = 19, cex = .5)
UZ_Scaffold29_ordered_S225 <- UZ_Scaffold29[order(-UZ_Scaffold29$S225_UZ),]
points(1:nrow(UZ_Scaffold29_ordered_S225), UZ_Scaffold29_ordered_S225$S225_UZ, ylim = c(0, 1), col= UZ_Scaffold29_ordered$color, pch = 19, cex = .5)

FST_m_z <- FST_m[FST_m$Scaffold %in% z_chr_scaffolds, ]
FST_m_auto <- FST_m[!FST_m$Scaffold %in% z_chr_scaffolds, ]


FST_m_auto <- FST_m_auto[(FST_m_auto$nSNPs >= 10 & FST_m_auto$fraction_covered >= 0.25), ]
FST_m_auto <- separate(data = FST_m_auto, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
FST_m_auto$Scaffold<-gsub("v02","",as.character(FST_m_auto$Scaffold))
FST_m_auto$length<-gsub("v02","",as.character(FST_m_auto$length))
FST_m_auto$length<-replace_na(FST_m_auto$length,"414000")

FST_m_auto <- FST_m_auto[as.numeric(FST_m_auto$length) >= 50000, ]
par(xpd=FALSE)


ma <- function(x, n = 50){stats::filter(x, rep(1 / n, n), sides = 2)}
fst_auto_N175_N225 <- FST_m_auto[, c(1:5,34)]
fst_auto_N175_N225$MA <- ma(fst_auto_N175_N225$N175_N225)
fst_auto_N175_N225$observation <- 1:nrow(fst_auto_N175_N225)

plot(1:nrow(fst_auto_N175_N225), fst_auto_N175_N225$N175_N225, pch = 19, cex = 0.3)
max(fst_auto_N175_N225$MA, na.rm = TRUE)
points(fst_auto_N175_N225$MA, type = 'l', col = 'red')
par(xpd=FALSE)
plot(1:nrow(fst_auto_N175_N225), fst_auto_N175_N225$N175_N225, pch = 19, cex = 0.3, xlim = c(224000, 226000))
points(fst_auto_N175_N225$MA, type = 'l', col = 'red')


fst_auto_N225_S225 <- FST_m_auto[, c(1:5, 45)]
fst_auto_N225_S225$MA <- ma(fst_auto_N225_S225$N225_S225)
fst_auto_N225_S225$observation <- 1:nrow(fst_auto_N225_S225)

plot(1:nrow(fst_auto_N225_S225), fst_auto_N225_S225$N225_S225, pch = 19, cex = 0.3, ylim = c(0,1))
max(fst_auto_N225_S225$N225_S225, na.rm = TRUE)
max(fst_auto_N225_S225$MA, na.rm = TRUE)
fst_auto_N225_S225_MA <- fst_auto_N225_S225[order(-fst_auto_N225_S225$MA), ]
points(fst_auto_N225_S225$MA, type = 'l', col = 'red')
plot(1:nrow(fst_auto_N225_S225), fst_auto_N225_S225$N225_S225, pch = 19, cex = 0.3, xlim = c(224000, 226000))
points(fst_auto_N225_S225$MA, type = 'l', col = 'red')

plot(1:nrow(fst_auto_N225_S225), fst_auto_N225_S225$N225_S225, pch = 19, cex = 0.3, xlim = c(170000, 180000))
points(fst_auto_N225_S225$MA, type = 'l', col = 'red')







# 
# fst_max_auto_N175_N225 <- FST_m_auto[, c(1:5,33)][order(-FST_m_auto[, 33]), ]
# fst_max_auto_N125_N225 <- FST_m_auto[, c(1:5,26)][order(-FST_m_auto[, 26]), ]
# fst_max_auto_N75_N225 <- FST_m_auto[, c(1:5,18)][order(-FST_m_auto[, 18]), ]
# fst_max_auto[1, 5] <- mean(FST_m_auto[, 40])
# fst_max_auto[1, 6] <- mean(FST_m_auto[, 42])
# fst_max_auto[1, 7] <- mean(FST_m_auto[, 43])
# fst_max_auto[1, 8] <- mean(FST_m_auto[, 44])
# fst_max_auto[1, 9] <- mean(FST_m_auto[, 45])
# fst_max_auto[2, 3] <- mean(FST_m_auto[, 25])
# fst_max_auto[2, 4] <- mean(FST_m_auto[, 17])
# fst_max_auto[2, 5] <- mean(FST_m_auto[, 34])
# fst_max_auto[2, 6] <- mean(FST_m_auto[, 36])
# fst_max_auto[2, 7] <- mean(FST_m_auto[, 37])
# fst_max_auto[2, 8] <- mean(FST_m_auto[, 38])
# fst_max_auto[2, 9] <- mean(FST_m_auto[, 39])
# fst_max_auto[3, 4] <- mean(FST_m_auto[, 16])
# fst_max_auto[3, 5] <- mean(FST_m_auto[, 27])
# fst_max_auto[3, 6] <- mean(FST_m_auto[, 29])
# fst_max_auto[3, 7] <- mean(FST_m_auto[, 30])
# fst_max_auto[3, 8] <- mean(FST_m_auto[, 31])
# fst_max_auto[3, 9] <- mean(FST_m_auto[, 32])
# fst_max_auto[4, 5] <- mean(FST_m_auto[, 19])
# fst_max_auto[4, 6] <- mean(FST_m_auto[, 21])
# fst_max_auto[4, 7] <- mean(FST_m_auto[, 22])
# fst_max_auto[4, 8] <- mean(FST_m_auto[, 23])
# fst_max_auto[4, 9] <- mean(FST_m_auto[, 24])
# fst_max_auto[5, 6] <- mean(FST_m_auto[, 47])
# fst_max_auto[5, 7] <- mean(FST_m_auto[, 48])
# fst_max_auto[5, 8] <- mean(FST_m_auto[, 49])
# fst_max_auto[5, 9] <- mean(FST_m_auto[, 50])
# fst_max_auto[6, 7] <- mean(FST_m_auto[, 55])
# fst_max_auto[6, 8] <- mean(FST_m_auto[, 56])
# fst_max_auto[6, 9] <- mean(FST_m_auto[, 57])
# fst_max_auto[7, 8] <- mean(FST_m_auto[, 58])
# fst_max_auto[7, 9] <- mean(FST_m_auto[, 59])
# fst_max_auto[8, 9] <- mean(FST_m_auto[, 60])
# 
# 
# 
# FST_pdfr <- FST_m[which(FST_m$Scaffold == "scaffold87|size517236"), ]
# 
# 
# model_1 <- lm(fst_max_auto_N175_N225$N175_N225~fst_max_auto_N175_N225$fraction_covered)
# plot(fst_max_auto_N175_N225$fraction_covered, fst_max_auto_N175_N225$N175_N225, cex = 0.5, pch = 19)
# 
# abline(model_1, col = 'red')
# 
# 
# plot(1:nrow(FST_m_z), FST_m_z$N225_S225, pch = 19, cex = 0.3)






FST_m_z <- FST_m[FST_m$Scaffold %in% z_chr_scaffolds, ]
FST_m_auto <- FST_m[!FST_m$Scaffold %in% z_chr_scaffolds, ]



FST_m_z <- FST_m_z[(FST_m_z$nSNPs >= 10 & FST_m_z$fraction_covered >= 0.25), ]
FST_m_z <- separate(data = FST_m_z, col = Scaffold, into = c("Scaffold", "length"), sep = "\\|size")
FST_m_z$Scaffold<-gsub("v02","",as.character(FST_m_z$Scaffold))
FST_m_z$length<-gsub("v02","",as.character(FST_m_z$length))
FST_m_z$length<-replace_na(FST_m_z$length,"414000")

z_linked <- c('scaffold4435', 'scaffold751', 'scaffold35', 'scaffold572', 'scaffold740', 'scaffold1546', 'scaffold247', 'scaffold196', 'scaffold779', 'scaffold2358', 'scaffold600', 'scaffold216', 'scaffold390', 'scaffold364', 'scaffold1895', 'scaffold409', 'scaffold444', 'scaffold325', 'scaffold15', 'scaffold169', 'scaffold348', 'scaffold865', 'scaffold518', 'scaffold2225', 'scaffold546', 'scaffold645', 'scaffold1618', 'scaffold789', 'scaffold218', 'scaffold1660', 'scaffold1044', 'scaffold581', 'scaffold1237', 'scaffold1246', 'scaffold1756', 'scaffold1161', 'scaffold155', 'scaffold282', 'scaffold804', 'scaffold502', 'scaffold3599', 'scaffold1040', 'scaffold275', 'scaffold4243', 'scaffold532', 'scaffold459', 'scaffold1889', 'scaffold221', 'scaffold931', 'scaffold140', 'scaffold1038', 'scaffold598', 'scaffold48', 'scaffold1486', 'scaffold193', 'scaffold137', 'scaffold87', 'scaffold50', 'scaffold662', 'scaffold728', 'scaffold229', 'scaffold330', 'scaffold673', 'scaffold812', 'scaffold1046')
FST_m_z <- FST_m_z %>% arrange(factor(Scaffold, levels = z_linked))


FST_m_z <- FST_m_z[as.numeric(FST_m_z$length) >= 50000, ]
par(xpd=FALSE)


fst_z_N225_S225 <- FST_m_z[, c(1:5, 45)]
fst_z_N225_S225$MA <- ma(fst_z_N225_S225$N225_S225)
fst_z_N225_S225$observation <- 1:nrow(fst_z_N225_S225)

plot(1:nrow(fst_z_N225_S225), fst_z_N225_S225$N225_S225, pch = 19, cex = 0.3, ylim = c(0,1))
max(fst_z_N225_S225$N225_S225, na.rm = TRUE)
max(fst_z_N225_S225$MA, na.rm = TRUE)
fst_z_N225_S225_MA <- fst_z_N225_S225[order(-fst_z_N225_S225$MA), ]
points(fst_z_N225_S225$MA, type = 'l', col = 'red')


plot(1:nrow(fst_z_N225_S225), fst_z_N225_S225$N225_S225, pch = 19, cex = 0.3, xlim = c(1400, 1600))
points(fst_z_N225_S225$MA, type = 'l', col = 'red')


fst_z_N75_S25 <- FST_m_z[, c(1:5, 20)]
fst_z_N75_S25$MA <- ma(fst_z_N75_S25$N75_S25)
fst_z_N75_S25$observation <- 1:nrow(fst_z_N75_S25)

plot(1:nrow(fst_z_N75_S25), fst_z_N75_S25$N75_S25, pch = 19, cex = 0.3, ylim = c(0,1))
max(fst_z_N75_S25$N75_S25, na.rm = TRUE)
max(fst_z_N75_S25$MA, na.rm = TRUE)
fst_z_N75_S25_MA <- fst_z_N75_S25[order(-fst_z_N75_S25$MA), ]
points(fst_z_N75_S25$MA, type = 'l', col = 'red')


plot(1:nrow(fst_z_N75_S25), fst_z_N75_S25$N75_S25, pch = 19, cex = 0.3, xlim = c(1400, 1600))
points(fst_z_N75_S25$MA, type = 'l', col = 'red')


fst_m_z_removals <- FST_m_z[, c(7:16, 21, 25, 29, 33, 36, 40, 42, 46:47, 51:55, 58, 60, 61)]
fst_m_z_kept <- FST_m_z[, c(1:6,17:20,22:24,26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59), ]

max(fst_m_z_kept[, 7:34])
#median(fst_m_z_kept[, 7:34])

fst_m_auto_removals <- FST_m_auto[, c(7:16, 21, 25, 29, 33, 36, 40, 42, 46:47, 51:55, 58, 60, 61)]
fst_m_auto_kept <- FST_m_auto[, c(1:6,17:20,22:24,26:28, 30:32, 34:35, 37:39, 41, 43:45, 48:50, 56:57, 59), ]
mat <- fst_m_auto_kept[, 7:34]

which(mat == max(mat), arr.ind=TRUE)
max(mat)

fst_max <- fst_m_auto_kept[122173, ]

fst_max_mat <- matrix(nrow = 8, ncol = 8)
fst_max_mat <- as.data.frame(matrix(nrow = 8, ncol = 8))

colnames(fst_max_mat) <- c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225")
row.names(fst_max_mat) <- c("N225", "N175", "N125", "N75", "S25", "S125", "S175", "S225")

diag(fst_max_mat) <- 0
fst_max_mat[1,2] <- fst_max[1,20]
fst_max_mat[1,3] <- fst_max[1,15]
fst_max_mat[1,4] <- fst_max[1,9]
fst_max_mat[1,5] <- fst_max[1,25]
fst_max_mat[1,6] <- fst_max[1,26]
fst_max_mat[1,7] <- fst_max[1,27]
fst_max_mat[1,8] <- fst_max[1,28]
fst_max_mat[2,3] <- fst_max[1,14]
fst_max_mat[2,4] <- fst_max[1,8]
fst_max_mat[2,5] <- fst_max[1,21]
fst_max_mat[2,6] <- fst_max[1,22]
fst_max_mat[2,7] <- fst_max[1,23]
fst_max_mat[2,8] <- fst_max[1,24]
fst_max_mat[3,4] <- fst_max[1,7]
fst_max_mat[3,5] <- fst_max[1,16]
fst_max_mat[3,6] <- fst_max[1,17]
fst_max_mat[3,7] <- fst_max[1,18]
fst_max_mat[3,8] <- fst_max[1,19]
fst_max_mat[4,5] <- fst_max[1,10]
fst_max_mat[4,6] <- fst_max[1,11]
fst_max_mat[4,7] <- fst_max[1,12]
fst_max_mat[4,8] <- fst_max[1,13]
fst_max_mat[5,6] <- fst_max[1,29]
fst_max_mat[5,7] <- fst_max[1,30]
fst_max_mat[5,8] <- fst_max[1,31]
fst_max_mat[6,7] <- fst_max[1,32]
fst_max_mat[6,8] <- fst_max[1,33]
fst_max_mat[7,8] <- fst_max[1,34]

fst_max_mat<- Matrix::forceSymmetric(as.matrix(fst_max_mat))
fst_max_mat<- as.matrix(fst_max_mat)

heatmap(fst_max_mat, Colv = "Rowv", Rowv = NA)
a<-c(3,4,6,2,3)
which(a==max(a))


scaff_254 <- fst_m_auto_kept[fst_m_auto_kept$Scaffold == 'scaffold254', ]
plot(scaff_254$Window, scaff_254$N175_S175, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N125_N175, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N75_N175, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N175_N225, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N175_S25, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N175_S125, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N175_S175, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N175_S225, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N225_S25, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N225_S125, type = 'l', ylim = c(0,1))
plot(scaff_254$Window, scaff_254$N225_S175, type = 'l', ylim = c(0,1))


medians <- c()
for(i in 7:34){
medians[i-6] <- median(fst_m_auto_kept[, i])
}
median(medians)
means <- c()
for(i in 7:34){
means[i-6] <- mean(fst_m_auto_kept[, i])
}
mean(means)

medians_z <- c()
for(i in 7:34){
medians_z[i-6] <- median(fst_m_z_kept[, i])
}
median(medians_z)
means_z <- c()
for(i in 7:34){
means_z[i-6] <- mean(fst_m_z_kept[, i])
}
mean(means_z)

#FST_recalc

FST_pdfr_UZ <- FST_pdfr[, c("Scaffold", "length", "Window", "nSNPs", "fraction_covered", "coverage_ave", "N225_UZ", "N175_UZ", "N125_UZ", "N75_UZ", "S25_UZ", "S125_UZ", "S175_UZ", "S225_UZ")]
FST_pdfr_UZ <- FST_pdfr_UZ[order(FST_pdfr_UZ$Window), ]
drop_cols <- c("N225_UZ", "N175_UZ", "N125_UZ", "N75_UZ", "S25_UZ", "S125_UZ", "S175_UZ", "S225_UZ", "N25_N225", "N225_S75", "N25_N175", "N175_S75", "N25_N125", "N125_S75", "N25_N75", "N75_S75", "N25_S25", "N25_S75", "N25_S125", "N25_S175", "N25_S225", "S75_S125", "S75_S125", "S75_S175", "S75_S225", "N25_UZ", "S75_UZ")

FST_pdfr_cline <- FST_pdfr[, !(names(FST_pdfr) %in% drop_cols)]


FST_exgenic_pdfr_UZ <- FST_pdfr_UZ[FST_pdfr_UZ$Window < 44000 | FST_pdfr_UZ$Window > 83000, ]
FST_genic_pdfr_UZ <- FST_pdfr_UZ[FST_pdfr_UZ$Window > 44000 & FST_pdfr_UZ$Window < 83000, ]


t.test(FST_exgenic_pdfr_UZ[, 7:14])
t.test(FST_genic_pdfr_UZ[, 7:14])




FST_period_UZ <- FST_period[, c("Scaffold", "length", "Window", "nSNPs", "fraction_covered", "coverage_ave", "N225_UZ", "N175_UZ", "N125_UZ", "N75_UZ", "S25_UZ", "S125_UZ", "S175_UZ", "S225_UZ")]
FST_period_UZ <- FST_period_UZ[order(FST_period_UZ$Window), ]
drop_cols <- c("N225_UZ", "N175_UZ", "N125_UZ", "N75_UZ", "S25_UZ", "S125_UZ", "S175_UZ", "S225_UZ", "N25_N225", "N225_S75", "N25_N175", "N175_S75", "N25_N125", "N125_S75", "N25_N75", "N75_S75", "N25_S25", "N25_S75", "N25_S125", "N25_S175", "N25_S225", "S75_S125", "S75_S125", "S75_S175", "S75_S225", "N25_UZ", "S75_UZ")

FST_period_cline <- FST_period[, !(names(FST_period) %in% drop_cols)]


FST_exgenic_period_UZ <- FST_period_UZ[FST_period_UZ$Window < 75000 | FST_period_UZ$Window > 145000, ]
FST_genic_period_UZ <- FST_period_UZ[FST_period_UZ$Window > 75000 & FST_period_UZ$Window < 145000, ]


t.test(FST_exgenic_period_UZ[, 7:14])
t.test(FST_genic_period_UZ[, 7:14])


#Z

FST_Z_UZ <- FST_m_z[, c("Scaffold", "Window", "nSNPs", "fraction_covered", "coverage_ave", "N225_UZ", "N175_UZ", "N125_UZ", "N75_UZ", "S25_UZ", "S125_UZ", "S175_UZ", "S225_UZ")]
FST_Z_UZ <- FST_m_z[order(FST_m_z$Window), ]
drop_cols <- c("N225_UZ", "N175_UZ", "N125_UZ", "N75_UZ", "S25_UZ", "S125_UZ", "S175_UZ", "S225_UZ", "N25_N225", "N225_S75", "N25_N175", "N175_S75", "N25_N125", "N125_S75", "N25_N75", "N75_S75", "N25_S25", "N25_S75", "N25_S125", "N25_S175", "N25_S225", "S75_S125", "S75_S125", "S75_S175", "S75_S225", "N25_UZ", "S75_UZ")

FST_Z_cline <- FST_m_z[, !(names(FST_m_z) %in% drop_cols)]

t.test(FST_Z_UZ[, 6:13])
t.test(FST_Z_cline[, 6:34])


#auto

FST_a_UZ <- FST_m_auto[, c("Scaffold", "Window", "nSNPs", "fraction_covered", "coverage_ave", "N225_UZ", "N175_UZ", "N125_UZ", "N75_UZ", "S25_UZ", "S125_UZ", "S175_UZ", "S225_UZ")]
FST_a_UZ <- FST_m_auto[order(FST_m_auto$Window), ]
drop_cols <- c("N225_UZ", "N175_UZ", "N125_UZ", "N75_UZ", "S25_UZ", "S125_UZ", "S175_UZ", "S225_UZ", "N25_N225", "N225_S75", "N25_N175", "N175_S75", "N25_N125", "N125_S75", "N25_N75", "N75_S75", "N25_S25", "N25_S75", "N25_S125", "N25_S175", "N25_S225", "S75_S125", "S75_S125", "S75_S175", "S75_S225", "N25_UZ", "S75_UZ")

FST_a_cline <- FST_m_auto[, !(names(FST_m_auto) %in% drop_cols)]

t.test(FST_a_UZ[, 6:13])
t.test(FST_a_cline[, 6:34])
