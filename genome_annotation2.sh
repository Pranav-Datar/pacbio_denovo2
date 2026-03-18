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
