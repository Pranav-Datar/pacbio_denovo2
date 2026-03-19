conda create -n repeatmodeler -c bioconda -c conda-forge repeatmodeler
conda activate repeatmodeler
BuildDatabase -name Vkomo_db GCF_004798865.1_ASM479886v1_genomic.fasta 2>&1 | tee repeatmodeler.log
RepeatModeler -threads 40 -database Vkomo_db/Vkomo_db 2>&1 | tee repeatmodeler.log

#download transcriptome data for komodo dragon
#download genome annotation file (gff) of komodo dragon from ncbi
gffread GCF_004798865.1_ASM479886v1_genomic.gff -T -o annotation.gtf

#transcriptome alignment and assembly
conda create -n star_env -c bioconda -c conda-forge star
conda activate star_env
STAR \
--runThreadN 50 \
--runMode genomeGenerate \
--genomeDir star_index \
--genomeFastaFiles genome/V_komo_genome.fasta \
--sjdbGTFfile genome/annotation.gtf \
--sjdbOverhang 149

STAR \
--runThreadN 50 \
--genomeDir star_index \
--readFilesIn transcriptome_data/SRR8735152_1.fastq transcriptome_data/SRR8735152_2.fastq \
--outFileNamePrefix sample1_ \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--quantMode TranscriptomeSAM

grep "Uniquely mapped reads %" sample1_Log.final.out
grep "% of reads mapped to multiple loci" sample1_Log.final.out
grep "% of reads unmapped" sample1_Log.final.out

mkdir star_output
#move all the files here


Trinity \
  --genome_guided_bam sample1_Aligned.sortedByCoord.out.bam \
  --genome_guided_max_intron 10000 \
  --max_memory 50G \
  --CPU 40

  #downloaded the protein sequences for reptiles mentioned in komodo dragon genome assembly paper. i went to each genome assembly in NCBI, and wget downloaded the protein.faa.gz file from ftp. the protein sequences were decompressed.

#make sure that there is no selenocysteine (U) amino acid in it. if its there, replace it with (C), which is the closest and workable substitute

#for aligning protein sequences to komodo genome,
exonerate --model protein2genome --query combined_reptile_proteins.faa --target /home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta --maxintron 10000 --
showtargetgff yes > exonerate.gff
