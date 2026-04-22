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
tabix -p vcf variants.vcf.gz

bcftools filter \
>   -i 'QUAL>=30 && FORMAT/DP>=4 && FORMAT/GQ>=10' \
>   variants.vcf.gz \
>   -Oz -o variants.filtered.vcf.gz

tabix -p vcf variants.filtered.vcf.gz
