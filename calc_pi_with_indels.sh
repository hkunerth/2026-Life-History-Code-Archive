# run pi calc from popoolation

perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N225.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output N225_pi1kb_ind.txt --snp-output N225_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N175.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output N175_pi1kb_ind.txt --snp-output N175_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N125.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output N125_pi1kb_ind.txt --snp-output N125_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N75.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output N75_pi1kb_ind.txt --snp-output N75_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input N25.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output N25_pi1kb_ind.txt --snp-output N25_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S25.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output S25_pi1kb_ind.txt --snp-output S25_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S75.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output S75_pi1kb_ind.txt --snp-output S75_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S125.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output S125_pi1kb_ind.txt --snp-output S125_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S175.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output S175_pi1kb_ind.txt --snp-output S175_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input S225.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output S225_pi1kb_ind.txt --snp-output S225_onlypi1kb_ind.snps
perl popoolation_1.2.2/Variance-sliding.pl --fastq-type sanger --measure pi --input UZ1.mpileup --min-count 2 --min-coverage 20 --max-coverage 200 --window-size 1000 --step-size 1000 --min-covered-fraction 0.2 --pool-size 80 --output UZ1_pi1kb_ind.txt --snp-output UZ1_onlypi1kb_ind.snps

command; echo "popoolation done" | mail -s "pi calc filtering done" hdk26@cornell.edu



