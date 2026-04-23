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

#Now, variant calling. But before that, ensure that you have Sorted BAM, Indexed BAM (.bai), Reference FASTA and Indexed FASTA (.fai)
conda activate longshot_env
longshot --bam V_komo_aligned.bam --ref V_komo_genome.fasta --out variants.vcf

conda activate bcftools_env

bgzip variants.vcf

bcftools consensus \
  -f V_komo_genome.fasta \
  -I \
  variants.filtered.vcf.gz \
  > consensus.fa
tabix -p vcf variants.vcf.gz

bcftools filter \
>   -i 'QUAL>=30 && FORMAT/DP>=4 && FORMAT/GQ>=10' \
>   variants.vcf.gz \
>   -Oz -o variants.filtered.vcf.gz

tabix -p vcf variants.filtered.vcf.gz

bcftools consensus \
>   -f V_komo_genome.fasta \
>   -I \
>   variants.filtered.vcf.gz \
>   > consensus.fa

seqtk seq -F 'I' consensus.fa > consensus.fq

fq2psmcfa -q20 consensus.fq > consensus.psmcfa

splitfa consensus.psmcfa > split.psmcfa

psmc -N25 -t10 -r5 -p "4+10*2+4+6" \
  -o result.psmc \
  consensus.psmcfa

  seq 1 100 | xargs -I{} -P 10 \
  psmc -N25 -t10 -r5 -b -p "4+10*2+4+6" \
  -o bootstrap_{}.psmc \
  split.psmcfa

cat result.psmc bootstrap_*.psmc > combined.psmc

psmc_plot.pl  -u 7.9e-9 -g 12 komodo combined.psmc

epstopdf komodo.eps
