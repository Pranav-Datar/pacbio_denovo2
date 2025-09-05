conda activate sratoolkit
prefetch SRR8797220 --output-directory /home/adminiisc/komodo_sra
#the sra file is downloaded in /home/adminiisc/komodo_sra

fasterq-dump SRR8797220.sra
#creates fastq file from sra file

#Now to gunzip it, do the following
conda install -c conda-forge pigz
pv SRR8797220.sra.fastq | pigz -c > SRR8797220.sra.fastq.gz
#pv: pipe viewer
#pigz = parallel implementation of gzip
#-c: write the compressed data to standard output (stdout), instead of replacing the original file

#create read length table
echo "platform,length" > length.csv
#creates platform and length columns in length.csv file

conda install bioconda::bioawk

pv SRR8797220.sra.fastq.gz | bioawk -c fastx '{print "PacBio_HiFi," length($seq)}' >> length.csv
#pv: pipe viewer
#pv SRR8797220.sra.fastq.gz: pv reads the zipped fastq file
#bioawk: version of awk that understands fasta/fastq files
#-c fastx: parses the input as fasta/fastq
#length($seq) = number of bases in that read.
#"PacBio_HiFi," is just a label (hardcoded string)
