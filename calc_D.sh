# run D calc from popoolation
#written to run on the cornell HPC

perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input N225.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N225_D1filtered.kb.txt --snp-output N225_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input N175.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N175_D1filtered.kb.txt --snp-output N175_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input N125.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N125_D1filtered.kb.txt --snp-output N125_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input N75.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N75_D1filtered.kb.txt --snp-output N75_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input N25.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N25_D1filtered.kb.txt --snp-output N25_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input S25.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S25_D1filtered.kb.txt --snp-output S25_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input S75.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S75_D1filtered.kb.txt --snp-output S75_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input S125.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S125_D1filtered.kb.txt --snp-output S125_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input S175.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S175_D1filtered.kb.txt --snp-output S175_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input S225.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S225_D1filtered.kb.txt --snp-output S225_onlyD1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure D --input UZ1.idf.filtered.mpileup --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output UZ1_D1filtered.kb.txt --snp-output UZ1_onlyD1kb.filtered.snps

command; echo "popoolation done" | mail -s "D calc filtering done" hdk26@cornell.edu


