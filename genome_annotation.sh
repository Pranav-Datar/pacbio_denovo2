conda create -n repeatmodeler -c bioconda -c conda-forge repeatmodeler
conda activate repeatmodeler
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

mkdir Vkomo_db
mv Vkomo_db* ~/Vkomo_db

RepeatModeler -database Vkomo_db/Vkomo_db -threads 40 -LTRStruct > out.log
#RepeatModeler: discovers transposable elements and repeats de novo
#database Vkomo_db: tells database to analyse the mentioned database
#LTRStruct: activates an additional LTR retrotransposon detection pipeline

#output:
#You will se many RM_* files. The first RM file generated, will have consensi.fa and consensi.fa.classified files, which are the main files.

grep -c ">" consensi.fa.classified
# This is to count the number of repeat families discovered.

RepeatMasker -pa 40 -lib RM_2500375.SunMar151620062026/consensi.fa.classified -gff -dir repeatmasker_out GCF_004798865.1_ASM479886v1_genomic.fasta
#RepeatMasker: Program that identifies and annotates repetitive elements in the genome.
#-pa 40: 40 threads
#-lib RM_2500375.SunMar151620062026/consensi.fa.classified: path the the consensi.fa.classified file which has the repeat library to be used in the analysis
#-gff: produce output in GFF format (general feature format). it consists of columns like chr1  RepeatMasker  LINE/L1   120034  120850
#-dir repeatmasker_out: output file directory
#GCF_004798865.1_ASM479886v1_genomic.fasta: input genome

conda create -n maker -c bioconda -c conda-forge maker

maker -CTL
#open maker_opts.ctl by 
nano maker_opts.ctl

#make sure the following:
#genome=put the fasta.masked from repeatmasker output
#model_org=remove "all"

#the remaining section should look like:
#model_org=
#rmlib=
#repeat_protein=
#softmask=1
#snaphmm=
#gmhmm=
#augustus_species=

#to save the file, Ctrl + O, enter and Ctrl + X.

#run maker
maker

#output: GCF_004798865.1_ASM479886v1_genomic.fasta_master_datastore_index.log. this file tracks all the annotation jobs

gff3_merge -d GCF_004798865.1_ASM479886v1_genomic.fasta_master_datastore_index.log
fasta_merge -d GCF_004798865.1_ASM479886v1_genomic.fasta_master_datastore_index.log
maker2zff -n GCF_004798865.1_ASM479886v1_genomic.fasta.all.gff
#it scans the datastore index, collects all gene annotations and combines them into one gff3 file
#look out for GCF_004798865.1_ASM479886v1_genomic.fasta.all.gff file. 

fasta_merge -d GCF_004798865.1_ASM479886v1_genomic.fasta_master_datastore_index.log
