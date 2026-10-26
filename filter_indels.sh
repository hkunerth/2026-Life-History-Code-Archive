# run indel filters from popoolation
#written to run on the cornell HPC

perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input N225.mpileup --gtf AllFiles.indels.gtf --output N225.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input N175.mpileup --gtf AllFiles.indels.gtf --output N175.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input N125.mpileup --gtf AllFiles.indels.gtf --output N125.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input N75.mpileup --gtf AllFiles.indels.gtf --output N75.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input N25.mpileup --gtf AllFiles.indels.gtf --output N25.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input S25.mpileup --gtf AllFiles.indels.gtf --output S25.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input S75.mpileup --gtf AllFiles.indels.gtf --output S75.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input S125.mpileup --gtf AllFiles.indels.gtf --output S125.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input S175.mpileup --gtf AllFiles.indels.gtf --output S175.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input S225.mpileup --gtf AllFiles.indels.gtf --output S225.idf.mpileup
perl popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input UZ1.mpileup --gtf AllFiles.indels.gtf --output UZ1.idf.mpileup

command; echo "popoolation done" | mail -s "indel filtering done" hdk26@cornell.edu
