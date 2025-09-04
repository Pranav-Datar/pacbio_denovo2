conda activate sratoolkit
prefetch SRR8797220 --output-directory /home/adminiisc/komodo_sra
#the sra file is downloaded in /home/adminiisc/komodo_sra

fasterq-dump SRR8797220.sra
#creates fastq file from sra file

gzip SRR8797220.sra.fastq
#converts fastq to fq.gz
