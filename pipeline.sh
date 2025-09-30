conda activate sratoolkit
prefetch SRR8797220 --output-directory /home/pranav/komodo_sra
#the sra file is downloaded in /home/adminiisc/komodo_sra

sam-dump SRR8797220.sra > SRR8797220.sam
#creates sam file from sra file

conda activate samtools
samtools view -bS SRR8797220.sam > SRR8797220.bam
#creates bam file from sam file
#view: The primary command for viewing and converting alignment files
#-b: Specifies that the output format should be BAM
#-S: input is sam format

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
wc -l length.csv
#prints number of lines in the file (20,59,680)

# Visualising read-length distribution

#download the length.csv on the local computer, and do the mentioned steps in R
setwd("//wsl.localhost/Ubuntu/home/pranavdatar")
library(ggplot2)
library(cowplot) #provides addition to ggplot2
library(scales) #map data to aesthetics

read_length_df <- read.csv("length.csv") #reads the data into the dataframe
mean_length <- mean(read_length_df$length) #calculates mean length

total_length_plot <- ggplot(read_length_df, aes(x = length)) +
  geom_histogram(binwidth = 100, fill = "#377eb8", color = "black", alpha = 0.5) +
  geom_vline(xintercept = mean_length, linetype = "dashed", color = "red", size = 0.3) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(x = "Read length (bp)", y = "Count",
       title = "PacBio HiFi Read Length Distribution") +
  theme_bw()
print(total_length_plot)

short_length_plot <- ggplot(read_length_df, aes(x = length)) +
  geom_histogram(binwidth = 50, fill = "#377eb8", color = "black", alpha = 0.5) +
  geom_vline(xintercept = mean_length, linetype = "dashed", color = "red", size = 0.3) +
  scale_x_continuous(labels = comma, limits = c(0, 20000)) +
  scale_y_continuous(labels = comma) +
  labs(x = "Read length (bp)", y = "Count",
       title = "PacBio HiFi Reads ≤ 20 kb") +
  theme_bw()
print(short_length_plot)

plot <- plot_grid(total_length_plot, short_length_plot, ncol = 1)
print(plot) #print both the plots in a single image

#assembly statistics
conda install -c conda-forge -c bioconda assembly-stats
conda update assembly-stats
assembly-stats SRR8797220.sra.fastq >> N50_stat
cat N50_stat 

#FastQC
conda install bioconda::fastqc
fastqc SRR8797220.sra.fastq.gz

#nanoplot
conda create -n nanoplot_env -c bioconda -c conda-forge nanoplot
conda activate nanoplot_env
NanoPlot --fastq SRR8797220.sra.fastq -o nanoplot_result

#adapter trimming
conda install -n base -c conda-forge mamba
mamba create -n hififilt_env -c bioconda -c conda-forge hifiadapterfilt
conda activate hififilt_env
hifiadapterfilt.sh -p SRR8797220 -l 44 -m 97 -o hifi_filtered SRR8797220.sra.fastq.gz ##make sure that the only input file with the prefix is in current directory. move all other files with the same prefix in some other directory, as this algorithm uses some loop with the prefix)
#hifiadapterfilt.sh: main script from the GitHub webpage
#-p prefix for output files
#-l 44: minimum adapter length match (default: 44. Reads must have ≥44 bp of adapter sequence to be flagged as contaminated.)
#-m 97: Minimum percent identity between the read segment and the adapter sequence (97%). This ensures only real adapter remnants are removed, not reads with coincidental matches.
#-o hifi_filtered output directory prefix
#SRR8797220.sra.fastq.gz: input file

#lenth trimming
conda activate chopper_env
chopper -l 1000 -i SRR8797220.sra.filt.fastq.gz | gzip > SRR8797220.sra.filt.lenfilt.fastq.gz
#filters reads less than 1000 bp

#QC post-trimming
conda activate nanoplot_env
NanoPlot --fastq SRR8797220.sra.filt.lenfilt.fastq.gz -o nanoplot_result_lengthfilt

#subsampling upto certain depth for few pilot runs
seqtk sample -s100 SRR16080541.sra.lenfilt10k.fastq 0.33 > SRR16080541.sra.lenfilt10k.subsampled0.33.fastq
#sample: subcommand that randomly selects a subset of reads
#-s: sets a random seed for reproducibility (here, 100).
#0.33: fraction of reads to be kept (here, 33%)

#assembly using hifiasm
hifiasm -o assembly_q7 -t 24 SRR8797220.sra.filt.lenfilt2qualfilt7.fastq 2> hifiasm_q7.log
#-o assembly_q7: states that "assembly_q7" would be the prefix of every output file
#-t 24: threads
#SRR8797220.sra.filt.lenfilt2qualfilt7.fastq: input file
#2> hifiasm_q7.log: stores all the running script in 2> "hifiasm_q7.log" file
