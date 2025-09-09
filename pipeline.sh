conda activate sratoolkit
prefetch SRR8797220 --output-directory /home/pranav/komodo_sra
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
