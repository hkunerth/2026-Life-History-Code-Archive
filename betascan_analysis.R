setwd("~/Desktop/Programs/betascan/")

#open betascan results
betascan <- read.csv("combined.betascores", header = TRUE, stringsAsFactors = FALSE, colClasses = c('numeric', 'numeric', 'character'), sep = '\t')
#make index column
betascan$index <- 1:nrow(betascan)
# filter betascan outliers
betascan_outliers <- betascan[betascan$Beta >= 30, ]
betascan_outliers <- na.omit(betascan_outliers)

# example scaffolds here, expand to more
scaff_184 <- betascan[betascan$scaffold1.betascores.txt == "scaffold184.betascores.txt",]
scaff_period <- betascan[betascan$scaffold1.betascores.txt == "scaffold532.betascores.txt",]
scaff_pdfr <- betascan[betascan$scaffold1.betascores.txt == "scaffold87.betascores.txt",]
scaff_275 <- betascan[betascan$scaffold1.betascores.txt == "scaffold275.betascores.txt",]
#plots and calculation of quantiles
plot(scaff_184$Position, scaff_184$Beta, pch = 19, cex = 0.3)
plot(scaff_period$Position, scaff_period$Beta, pch = 19, cex = 0.3, ylim = c(0, 18))
quantile(scaff_period$Beta, na.rm = TRUE, probs = seq(0, 1, 0.01))
top_pdfr <- scaff_pdfr[scaff_pdfr$Beta > 14.5, ]
plot(scaff_pdfr$Position, scaff_pdfr$Beta, pch = 19, cex = 0.3, ylim = c(0, 23))
quantile(scaff_pdfr$Beta, na.rm = TRUE, probs = seq(0, 1, 0.01))
plot(scaff_275$Position, scaff_275$Beta, pch = 19, cex = 0.3, ylim = c(0, 30))
plot(scaff_275$Position, scaff_275$Beta, pch = 19, cex = 0.3, ylim = c(0, 30), xlim = c(750000, 800000))
quantile(scaff_period$Beta, na.rm = TRUE, probs = seq(0, 1, 0.01))
top_period <- scaff_period[scaff_period$Beta > 14.58, ]
plot(scaff_period$Position, scaff_period$Beta, pch = 19, cex = 0.3, ylim = c(0, 18), xlim = c(70000, 150000))

betascan_outliers_275 <- scaff_275 [scaff_275$Beta >= 25, ]
betascan_outliers <- na.omit(betascan_outliers)

# changing 0s to NAs for plotting 
betascan_sub <- betascan[betascan$Beta >= 0, ]
betascan_sub$change <- NA
count = 0
for (i in 1:(nrow(betascan_sub)))
for (i in 1:nrow(betascan_sub))
{
  betascan_sub[i, 5] <- count %% 2
  if(is.na(betascan_sub[i, 1]) == TRUE) {count <- count+1}
  
  
}

match_df <- matrix(nrow = 2, ncol = 2)
match_df[, 1] <- c(1, 0)
match_df[, 2] <- c("black", "white")
colnames(match_df) <- c("V5", "color")
# example plot
pdf("~/Desktop/scaff_184.pdf", height = 4, width = 6)
plot(scaff_184$Position, scaff_184$Beta, pch = 19, cex = 0.3, ylim = c(0,41), main = "Betascan of scaffold 184", ylab = "Betascore", xlab = "Position on scaffold")
dev.off()
