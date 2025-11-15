conda create -n sratoolkit-c bioconda -c conda-forge sra-tools
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

filtlong --target_bases 60000000000 SRR16080541.fastq.gz > SRR16080541_40x.fastq

#assembly using hifiasm
hifiasm -o assembly_q7 -t 24 SRR8797220.sra.filt.lenfilt2qualfilt7.fastq 2> hifiasm_q7.log
#-o assembly_q7: states that "assembly_q7" would be the prefix of every output file
#-t 24: threads
#SRR8797220.sra.filt.lenfilt2qualfilt7.fastq: input file
#2> hifiasm_q7.log: stores all the running script in 2> "hifiasm_q7.log" file

#converting the assembled genome file to fasta format
awk '/^S/{print ">"$2"\n"$3}' assembly.bp.p_ctg.gfa > primary.fasta


#quality check post-assembly
conda create -n quast_env python=3.8 -y
conda activate quast_env
conda install -c bioconda -c conda-forge quast -y

quast primary.fasta -o quast_primary



