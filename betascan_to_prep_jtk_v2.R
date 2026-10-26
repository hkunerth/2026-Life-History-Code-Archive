# visualize betascan output
library(tidyr)

setwd("~/Desktop/Life_history/poolseq_analysis_v2/")



#read the primary data file
rc_file <- read.csv("AllFiles.14.200.m10.idf_rc", sep = "\t", header = TRUE)
#split up scaffold names
rc_file <- separate(data = rc_file, col = X..chr, into = c("Scaffold", "length"), sep = "\\|size")
#handle the bric-a-brac scaffold which is named differently
rc_file$Scaffold<-gsub("v02","",as.character(rc_file$Scaffold))
rc_file$length<-gsub("v02","",as.character(rc_file$length))
rc_file$length<-replace_na(rc_file$length,"414000")
#length filter
rc_file <- rc_file[which(as.numeric(rc_file$length) >= 5000), ]

# add an index column
rc_file$observation <- 1:nrow(rc_file)
#filter for  basepairs
selected = c("A", "T", "C", "G")
rc_file <- rc_file[rc_file$rc %in% selected, ]
#bind outputs, could be cleaned up but it works
rc_file <- rc_file[, c(1:4, 9:20, 22:31, 33)]
rc_file <- separate(data = rc_file, col = maa_1, into = c("maa1", "maa1x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_2, into = c("maa2", "maa2x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_3, into = c("maa3", "maa3x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_4, into = c("maa4", "maa4x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_5, into = c("maa5", "maa5x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_6, into = c("maa6", "maa6x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_7, into = c("maa7", "maa7x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_8, into = c("maa8", "maa8x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_9, into = c("maa9", "maa9x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_10, into = c("maa10", "maa10x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_1, into = c("mia1", "mia1x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_2, into = c("mia2", "mia2x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_3, into = c("mia3", "mia3x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_4, into = c("mia4", "mia4x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_5, into = c("mia5", "mia5x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_6, into = c("mia6", "mia6x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_7, into = c("mia7", "mia7x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_8, into = c("mia8", "mia8x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_9, into = c("mia9", "mia9x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_10, into = c("mia10", "mia10x"), sep = "/")

#rc_file <- head(rc_file, -1)


rc_file = subset(rc_file, select = -c(maa1x, maa2x, maa3x, maa4x, maa5x, maa6x, maa7x, maa8x, maa9x, maa10x, mia1x, mia2x, mia3x, mia4x, mia5x, mia6x, mia7x, mia8x, mia9x, mia10x) )
for(y in 7:27){
rc_file[, y] <- as.numeric(rc_file[, y])
}

nrows = nrow(rc_file)
results <- matrix(0, nrows, 10)
#convert file to results that JTKcycle can use
#res_alt <- numeric(nrow(rc_file))
for(j in 1:nrows){
for(i in c(1:10)){
       if(rc_file$rc[j] == substring(rc_file$major_alleles.maa.[j], i, i)){ 
           results[j, i] <- (rc_file[j, i+6] / (rc_file[j, i+6] + rc_file[j, i+16]))
       }
      if(rc_file$rc[j] == substring(rc_file$minor_alleles.mia.[j], i, i)){ 
           results[j, i] <- (rc_file[j, i+16] / (rc_file[j, i+6] + rc_file[j, i+16]))
       }
  }

  if(j %% 10000 == 0){print(j)
    print(Sys.time())
    }
}
# label and order results
results_labeled <- cbind(rc_file$observation, results)
colnames(results_labeled) <- c("Index", "N25", "N75", "N125", "N175", "N225", "S25", "S75", "S125", "S175", "S225")
results_ordered <- cbind(results_labeled[, 1], results_labeled[, 6], results_labeled[, 5], results_labeled[, 4], results_labeled[, 3], results_labeled[, 2], results_labeled[, 7:11])
colnames(results_ordered) <- c("Index", "N225", "N175", "N125", "N75", "N25", "S25", "S75", "S125", "S175", "S225")
#make annotation file
annot <- rc_file[, 1:4]

all_res <- cbind(annot, results_ordered)
all_res[all_res == 0] <- NA
all_res <- all_res[complete.cases(all_res), ]

annot_out <- all_res[, 1:4]
data_out <- all_res[, 5:15]
# write results
write.table(annot_out, "ecb_9_1.annot", quote = FALSE, sep = '\t')
write.table(data_out, "ecb_9_1.data", quote = FALSE, sep = '\t')

