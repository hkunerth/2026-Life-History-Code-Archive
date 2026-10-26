# script to look for sinusoidal patterns in allele frequency data
#written to run on the cornell HPC

setwd("/workdir/hdk/")
library(pbapply)
source("JTK_CYCLEv3.1.R")

project <- "poolseq_cbsu_all_scaffolds"

options(stringsAsFactors=FALSE)
annot <- read.delim("ecb_9_1.annot")
data <- read.delim("ecb_9_1.data")

rownames(data) <- 1:nrow(data)
data <- data[,-1]
data <- cbind(data, annot)
data[data == 0] <- NA
data <- data[complete.cases(data), ]
annot <- as.data.frame(data[, 11])
data <- data[, -11]


jtkdist(ncol(data))

periods <- 4:6
jtk.init(periods,1);

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- pbapply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res,data)
  results <- results[order(res$ADJ.P,-res$AMP),]

})
print(st)

save(results,file=paste("JTK",project,"rda",sep="."))
write.table(results,file=paste("JTK",project,"txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")

