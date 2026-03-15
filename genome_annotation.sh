conda create -n repeatmodeler -c bioconda -c conda-forge repeatmodeler

BuildDatabase -name Vkomo_db GCF_004798865.1_ASM479886v1_genomic.fasta
#converts genome.fasta into an indexed database, and prepares it for fast discovery
#BuildDatabase: program included within RepeatModeler
#-name Vkomo_db: constructs database "Vkomo_db"
#GCF_004798865.1_ASM479886v1_genomic.fasta: input genome sequence

#output:
#Vkomo_db.nsq: contains the actual genome sequence data stored in database format
#Vkomo_db.nhr: contains sequence metadata
#Vkomo_db.nin: stores the index of sequences 
#Vkomo_db.nnd: additional indexing information
#Vkomo_db.nni: nucleotide index index
#Vkomo_db.nog: sequence grouping information
#Vkomo_db.njs: job scheduling metadata
#Vkomo_db.translation: contains translation information
