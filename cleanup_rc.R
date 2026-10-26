# visualize betascan output
library(tidyr)

setwd("~/Desktop/Programs/betascan/")




# grab largest 600 scaffolds
for(i in 5:600) {
  infile <- paste("~/Desktop/for_betascan/scaffold",i,sep="")
  outfile <- paste("scaffold",i,".betascan",sep="")

# open main rc file and rename columns
rc_file <- read.csv(infile, sep = "\t", header = TRUE)
colnames(rc_file) <- c('chr',	'pos',	'rc',	'allele_count',	'allele_states',	'deletion_sum',	'snp_type',	'major_alleles_maa',	'minor_alleles_mia',	'maa_1',	'maa_2',	'maa_3',	'maa_4',	'maa_5',	'maa_6',	'maa_7',	'maa_8',	'maa_9',	'maa_10',	'maa_11',	'mia_1',	'mia_2',	'mia_3',	'mia_4',	'mia_5',	'mia_6',	'mia_7',	'mia_8',	'mia_9',	'mia_10',	'mia_11')
# clean up scaffolds and rename bric-a-brac scaffold
rc_file <- separate(data = rc_file, col = chr, into = c("Scaffold", "length"), sep = "\\|size")
rc_file$Scaffold<-gsub("v02","",as.character(rc_file$Scaffold))
rc_file$length<-gsub("v02","",as.character(rc_file$length))
rc_file$length<-replace_na(rc_file$length,"414000")


#filter basepairs
rc_file$observation <- 1:nrow(rc_file)
selected = c("A", "T", "C", "G")
rc_file <- rc_file[rc_file$rc %in% selected, ]
# separate major and minor alleles based on column
rc_file <- rc_file[, c(1:4, 9:10, 12:20, 23:31, 33)]
rc_file <- separate(data = rc_file, col = maa_2, into = c("maa2", "maa2x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_3, into = c("maa3", "maa3x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_4, into = c("maa4", "maa4x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_5, into = c("maa5", "maa5x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_6, into = c("maa6", "maa6x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_7, into = c("maa7", "maa7x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_8, into = c("maa8", "maa8x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_9, into = c("maa9", "maa9x"), sep = "/")
rc_file <- separate(data = rc_file, col = maa_10, into = c("maa10", "maa10x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_2, into = c("mia2", "mia2x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_3, into = c("mia3", "mia3x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_4, into = c("mia4", "mia4x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_5, into = c("mia5", "mia5x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_6, into = c("mia6", "mia6x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_7, into = c("mia7", "mia7x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_8, into = c("mia8", "mia8x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_9, into = c("mia9", "mia9x"), sep = "/")
rc_file <- separate(data = rc_file, col = mia_10, into = c("mia10", "mia10x"), sep = "/")

#covert to numeric
for(y in 7:43){
rc_file[, y] <- as.numeric(rc_file[, y])
}
#subset file
rc_file = subset(rc_file, select = -c(maa2x, maa3x, maa4x, maa5x, maa6x, maa7x, maa8x, maa9x, maa10x, mia2x, mia3x, mia4x, mia5x, mia6x, mia7x, mia8x, mia9x, mia10x) )

res_ref <- numeric(nrow(rc_file))

res_alt <- numeric(nrow(rc_file))
#loop through all SNPs and calculate allele frequencies
for(j in 1:nrow(rc_file)){
for(i in c(1:9)){
       if(substring(rc_file$major_alleles_maa[j], i+1, i+1) != 'N'){
          (ifelse((rc_file$rc[j] == substring(rc_file$major_alleles_maa[j], i+1, i+1)), res_ref[j] <- rc_file[j, i+6] + res_ref[j] , res_alt[j] <- rc_file[j, i+6] + res_alt[j]))
       }
         if(substring(rc_file$minor_alleles_mia[j], i+1, i+1) != 'N'){ 
           (ifelse((rc_file$rc[j] == substring(rc_file$minor_alleles_mia[j], i+1, i+1)), res_ref[j] <- rc_file[j, i+15] + res_ref[j] , res_alt[j] <- rc_file[j, i+15] + res_alt[j]))
         }
}
  if(j %% 10000 == 0){print(j)
    print(Sys.time())
    }
}
# write outputs
output <- cbind(rc_file$pos, res_alt, (res_alt+res_ref))
write.table(output, outfile, quote = FALSE, col.names = FALSE, row.names = FALSE)
}
