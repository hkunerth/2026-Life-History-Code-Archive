#example of using popoolation to calculate PI for a single population
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N175.idf.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 120 --output N175_pi1kb.txt --snp-output N175_onlypi1kb.snps
