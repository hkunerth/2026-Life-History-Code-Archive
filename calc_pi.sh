# run pi calc from popoolation
#written to run on the cornell HPC

perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N225.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N225_pi1filtered.kb.txt --snp-output N225_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N175.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N175_pi1filtered.kb.txt --snp-output N175_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N125.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N125_pi1filtered.kb.txt --snp-output N125_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N75.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N75_pi1filtered.kb.txt --snp-output N75_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N25.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output N25_pi1filtered.kb.txt --snp-output N25_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S25.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S25_pi1filtered.kb.txt --snp-output S25_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S75.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S75_pi1filtered.kb.txt --snp-output S75_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S125.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S125_pi1filtered.kb.txt --snp-output S125_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S175.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S175_pi1filtered.kb.txt --snp-output S175_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S225.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output S225_pi1filtered.kb.txt --snp-output S225_onlypi1kb.filtered.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input UZ1.idf.filtered.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 50 --output UZ1_pi1filtered.kb.txt --snp-output UZ1_onlypi1kb.filtered.snps

command; echo "popoolation done" | mail -s "pi calc filtering done" hdk26@cornell.edu



