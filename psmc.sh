##Indexing the genome

conda activate pbmm2_env
pbmm2 index V_komo_genome.fasta reference.mmi --preset HIFI
#V_komo_genome.fasta: input genome
#reference.mmi: output index

#aligning the raw pacbio reads to the genome
pbmm2 align reference.mmi SRR8797220.fastq V_komo_aligned.bam --preset HIFI --sort
#SRR8797220.fastq: input reads
#V_komo_aligned_bam: output
#--sort: for sorting
conda deactivate
#Now, variant calling. But before that, ensure that you have Sorted BAM, Indexed BAM (.bai), Reference FASTA and Indexed FASTA (.fai)

#remove sex chromosomes
#first, calculate the average genome depth
conda activate qualimap_env
qualimap bamqc -bam V_salvator_aligned.bam -outdir qualimap_out -nt 40 --java-mem-size=50G

#second, calculate per contig depth
conda activate samtools
samtools coverage V_salvator_aligned.bam > coverage.tsv

awk '$3 > 100000' coverage.tsv > large_scaffolds.tsv #keep only large contigs as smaller ones can bring in noise.

awk '{print $1, $7}' large_scaffolds.tsv > scaffold_depths.tsv #create a new file with scaffold name and its depth

conda activate r_env
R
x <- read.table("scaffold_depths.tsv", header=FALSE)
hist(x$V2, breaks=40)
q()
#save the image
#if 2 peaks then the smaller peak could correspond to W chromosome
#the sex chromosomes (especially the W chromosome in lizards) will have ~half the total coverage

#now, specifically for Z chromosome, since it is diploid in males, but still with different effective population size,

#download the komodo dragon genome assembly report, which contains information about sex chromosome contigs
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/798/865/GCF_004798865.1_ASM479886v1/GCF_004798865.1_ASM479886v1_assembly_report.txt

#filter out those contig names and save them in a .txt file
grep -E "^scaffold(79|80|89|99|113|164)[[:space:]]" komodo_assembly_report.txt | awk '{print $7}' > komodo_Z_ids.txt

#filter out the sequences for those sex chromosome contigs and save them in a .fasta file
conda activate samtools
samtools faidx V_komo_genome.fasta $(cat komodo_Z_ids.txt) > komodo_Z.fa

#note: komodo_Z.fa is in /home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome

#figure out the z-linked scaffolds in your assembly 
conda activate minimap2
minimap2 -x asm5 -t 40 reference.mmi /home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/komodo_Z.fa > z_alignments.paf
#-x asm5: assembly-to-assembly alignment. Note: minimap2 is good for general sequence alignment while pbmm2 is specialized for pacbio read alignment

#Extract candidate salvator Z scaffolds, filter strong alignments. This gives candidate Z scaffolds in our genome assembly
awk '$11 > 100000 {print $6}' z_alignments.paf | sort | uniq > salvator_Z_scaffolds.txt

#now create a file which contains all the sex linked scaffolds, which are to be removed from the input genome.
#for W scaffolds
awk '$3 > 100000 && $7 > 27 && $7 < 37 {print $1}' coverage.tsv > salvator_W_scaffolds.txt
#This keeps scaffold length > 100000, and depth between 27X-37X (since i saw the smaller peak corresponding to W scaffolds at 27-37X the the graph plotted above.

#For indexing the genome,
conda activate samtools
samtools faidx #genome file
conda deactivate

#variant calling for finding out where are the variants across the diploid genome
conda activate longshot_env
longshot --bam V_komo_aligned.bam --ref V_komo_genome.fasta --out variants.raw.vcf
conda deactivate

conda activate bcftools_env

bgzip variants.raw.vcf > variants.raw.vcf.gz
#for compressing the vcf file

#index the vcf file
tabix -p vcf variants.raw.vcf.gz

bcftools filter \
>   -i 'QUAL>=30 && FORMAT/DP>=40 && FORMAT/GQ>=10 && MQ>=20' \
>   variants.raw.vcf.gz \
>   -Oz -o variants.filtered.vcf.gz

#-i: keeps variants only if the below conditions are true
#QUAL>=30 site level confidence 
#FORMAT/DP>=4  depth per sample (number of reads covering per position)
#FORMAT/GQ>=10 genotype quality
#Oz: output in compressed format

#index the filtered vcf file
tabix -p vcf variants.filtered.vcf.gz

#now using the filtered vcf file, make a diploid consensus genome for psmc which will have the info about homozygotes and heterozygotes based on IUPAC codes
bcftools consensus \
>   -f V_komo_genome.fasta \
>   -I \
>   variants.filtered.vcf.gz \
>   > consensus.fa

conda deactivate

#psmc requires fastq file and not fasta file. so add dummy score "I" and convert it to fastq
conda activate seqtk_env
seqtk seq -F 'I' consensus.fa > consensus.fq
conda deactivate

#convert to psmc input file
conda activate psmc_env
fq2psmcfa -q20 consensus.fq > consensus.psmcfa
#-q0 for preventing unecessary masking as qualities are anyways fake

#split the file for bootstrapping
splitfa consensus.psmcfa > split.psmcfa

psmc -N25 -t10 -r5 -p "4+10*2+4+6" \
  -o result.psmc \
  consensus.psmcfa

  #N25: 25 iterations
  #t10: the upper timing bound
  #-r5: recombination/mutation ratio
  #-p "4+10*2+4+6": division of the time bins 

  seq 1 100 | xargs -I{} -P 10 \
  psmc -N25 -t10 -r5 -b -p "4+10*2+4+6" \
  -o bootstrap_{}.psmc \
  split.psmcfa

  #rereuns psmc 100 times
  #P 10: runs 10 jobs at a time
  #-b: bootstrap run

cat result.psmc bootstrap_*.psmc > combined.psmc

psmc_plot.pl  -u 7.9e-9 -g 12 komodo combined.psmc

epstopdf komodo.eps
