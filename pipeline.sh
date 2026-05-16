#downloading files from google drive
gdown --folder "google_drive_folder_link"

#converting bam to fastq
conda activate samtools
samtools fastq input.bam > output.fastq

conda create -n sratoolkit -c bioconda -c conda-forge sra-tools
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

#adapter trimming: not compulsorily required since pacbio sequencer itself cuts off the adapter sequences. HiFiadapterfilt not working, giving database error. try using cutadapt
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

#Count total barcode hits in the sequences
seqkit locate -p TGCATACTGCGAGTAT 1_A01_output.fastq | wc -l
#divide the output by the total number of reads from nanoplot output, which will give the % of reads with the respective sequence present

#selectively trim out the barcode sequences found at the start and the end of the reads. This is optional, just for eyeballing
seqkit locate -p "barcode_sequence" brevicauda_combined.fastq  > barcode_hits.tsv #get the positions of adapter sequence existence
seqkit fx2tab -n -l brevicauda_combined.fastq > read_lengths.tsv #extract out the read lengths
#classify the positions of adapters; the start ones and the end ones. this step is optional and not compulsory
awk 'NR==FNR{len[$1]=$2; next}
{
  start=$5; end=$6; L=len[$1];

  if (start <= 30) print $1 > "start_hits.txt";
  else if (end >= (L-30)) print $1 > "end_hits.txt";
}' read_lengths.tsv barcode_hits.tsv

#trim out the start and end adapter sequences
conda activate cutadapt_env

cutadapt \
-j 16 \
-g barcode sequence \
-a barcode sequence \
-O 14 \
-e 0.05 \
-o brevicauda_combined_adaptertrimmed.fastq \
brevicauda_combined.fastq
#-j 16: 16 threads
#-g: read start
#-a: read end
#-O 14: only trim in atleast 14 bases match
#-e 0.05: allow upto 1% mismatches

#correct command for barcode trimming only from the start and end of the read
cutadapt \
  -j 16 \
  -g ^TGCATACTGCGAGTAT \
  -a TGCATACTGCGAGTAT$ \
  -O 14 \
  -e 0.05 \
  -o brevicauda_combined_adaptertrimmed.fastq \
  brevicauda_combined.fastq


#universal blunt adapter trimming

#check if it is present
seqkit locate -p ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT input.fastq | head

cutadapt \
  -j 16 \
  -g ^ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT \
  -g ^ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT \
  -a ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT$ \
  -a ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT$ \
  -O 40 \
  -e 0.05 \
  --no-indels \
  -o reads_adaptertrimmed.fastq \
  input.fastq

#length trimming
conda activate chopper_env
chopper -l 1000 -i SRR8797220.sra.filt.fastq.gz | gzip > SRR8797220.sra.filt.lenfilt.fastq.gz
#filters reads less than 1000 bp

#contamination check

#run the fastq file through kraken2
conda activate kraken2
 kraken2 \
--db ~/k2_pluspfp_db \
--threads 48 \
--report 1_B01.kraken.report.txt \
--output 1_B01.kraken.output.txt \
DRR818294.fastq

#interpreting the output
awk '
$4=="U"{unclass=$2}
$6=="root"{root=$2}
$6=="Archaea"{arch=$2}
$6=="Bacteria"{bac=$2}
$6=="Viruses"{vir=$2}
$6=="Fungi"{fun=$2}
$6=="Viridiplantae"{pla=$2}
$6=="Homo" && $7=="sapiens"{hum=$2}
END{
total=unclass+root
printf "Category\tPercent\n"
printf "Unclassified\t%.4f%%\n",(unclass/total)*100
printf "Archaea\t%.4f%%\n",(arch/total)*100
printf "Bacteria\t%.4f%%\n",(bac/total)*100
printf "Viruses\t%.4f%%\n",(vir/total)*100
printf "Fungi\t%.4f%%\n",(fun/total)*100
printf "Plants\t%.4f%%\n",(pla/total)*100
printf "Homo_sapiens\t%.4f%%\n",(hum/total)*100
}' 1_B01.kraken.report.txt

#QC post-trimming
conda activate nanoplot_env
NanoPlot --fastq SRR8797220.sra.filt.lenfilt.fastq.gz -o nanoplot_result_lengthfilt

#subsampling upto certain depth for few pilot runs
seqtk sample -s100 SRR16080541.sra.lenfilt10k.fastq 0.33 > SRR16080541.sra.lenfilt10k.subsampled0.33.fastq
#sample: subcommand that randomly selects a subset of reads
#-s: sets a random seed for reproducibility (here, 100).
#0.33: fraction of reads to be kept (here, 33%)

filtlong --target_bases 60000000000 SRR16080541.fastq.gz > SRR16080541_40x.fastq
#selectively keeps top reads based on length and quality until it is 60 Gb

#assembly using hifiasm
hifiasm -o assembly_q7 -t 24 SRR8797220.sra.filt.lenfilt2qualfilt7.fastq 2> hifiasm_q7.log
#-o assembly_q7: states that "assembly_q7" would be the prefix of every output file
#-t 24: threads
#SRR8797220.sra.filt.lenfilt2qualfilt7.fastq: input file
#2> hifiasm_q7.log: stores all the running script in 2> "hifiasm_q7.log" file

#converting the assembled genome file to fasta format
awk '/^S/{print ">"$2"\n"$3}' assembly.bp.p_ctg.gfa > primary.fasta
#it reads the gfa file and writes out the contig sequences in fasta format

#quality check post-assembly
conda create -n quast_env python=3.8 -y
conda activate quast_env
conda install -c bioconda -c conda-forge quast -y

quast primary.fasta -o quast_primary

# BUSCO
conda create -n busco_env -c conda-forge -c bioconda busco=6.0.0
conda activate busco_env
busco -i primary.fasta -l sauropsida_odb12 -m genome -o busco_primary

#compleasm (it is faster than BUSCO, but not recommended for distant genome assemblies (non-model organisms)
conda create -n compleasm_env -c conda-forge -c bioconda compleasm
conda activate compleasm_env
compleasm run -a primary.fasta -o compleasm_output_homo -l primates

#snail plot
conda activate blobtoolkit_env
nano meta.json
#write the following in it:
#{
#  "id": "komodo_genome",
#  "name": "Varanus komodoensis"
#}

blobtools create \
  --fasta /home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta \
  --meta meta.json \
  my_dataset

blobtools view \
  --plot \
  --view snail \
  my_dataset

