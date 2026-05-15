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

#the sex chromosomes (especially the W chromosome in lizards) will have ~half the total coverage

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
