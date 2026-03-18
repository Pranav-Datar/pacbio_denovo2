conda create -n repeatmodeler -c bioconda -c conda-forge repeatmodeler
conda activate repeatmodeler
BuildDatabase -name Vkomo_db GCF_004798865.1_ASM479886v1_genomic.fasta 2>&1 | tee repeatmodeler.log
RepeatModeler -threads 40 -database Vkomo_db/Vkomo_db 2>&1 | tee repeatmodeler.log

wget https://www.dfam.org/releases/Dfam_3.8/families/Dfam.embl.gz
