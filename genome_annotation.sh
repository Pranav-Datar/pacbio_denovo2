conda create -n repeatmodeler -c bioconda -c conda-forge repeatmodeler

BuildDatabase -name Vkomo_db GCF_004798865.1_ASM479886v1_genomic.fasta
#converts genome.fasta into an indexed database, and prepares it for fast discovery
#BuildDatabase: program included within RepeatModeler
#-name Vkomo_db: constructs database "Vkomo_db"
#GCF_004798865.1_ASM479886v1_genomic.fasta: input genome sequence
