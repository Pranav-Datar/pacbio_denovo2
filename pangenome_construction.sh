#PanSV naming
seqkit replace \
-p '^(\S+)' \
-r 'acanthurus#0#${1}' \
/home/pranav/genome_assemblies/NCBI_data/V_acanthurus_ncbi/GCA_050042745.1_ZJU_Vac_2.0_genomic.fasta \
> acanthurus.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'komodoensis#0#${1}' \
/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta \
> komodoensis.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'salvator#0#${1}' \
/home/pranav/genome_assemblies/NCBI_data/V_salvator_2_ncbi/genome/GWHHKDC00000000.1.genome.fasta \
> salvator.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'brevicauda#1#${1}' \
/home/pranav/genome_assemblies/primary_data/V_brevicauda/brevicauda_combined_reanalysis2/assembly_files/brevicauda.assembly.bp.hap1.p_ctg.fasta \
> brevicauda.hap1.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'brevicauda#2#${1}' \
/home/pranav/genome_assemblies/primary_data/V_brevicauda/brevicauda_combined_reanalysis2/assembly_files/brevicauda.assembly.bp.hap2.p_ctg.fasta \
> brevicauda.hap2.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'panoptes#1#${1}' \
/home/pranav/genome_assemblies/primary_data/V_panoptes/panoptes_reanalysis2/assembly_files/panoptes.assembly.bp.hap1.p_ctg.fasta \
> panoptes.hap1.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'panoptes#2#${1}' \
/home/pranav/genome_assemblies/primary_data/V_panoptes/panoptes_reanalysis2/assembly_files/panoptes.assembly.bp.hap2.p_ctg.fasta \
> panoptes.hap2.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'scalaris#1#${1}' \
/home/pranav/genome_assemblies/primary_data/V_scalaris/scalaris_reanalysis2/assembly_files/scalaris.assembly.bp.hap1.p_ctg.fasta \
> scalaris.hap1.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'scalaris#2#${1}' \
/home/pranav/genome_assemblies/primary_data/V_scalaris/scalaris_reanalysis2/assembly_files/scalaris.assembly.bp.hap2.p_ctg.fasta \
> scalaris.hap2.pansn.fa




conda create -n pggb_env -c conda-forge -c bioconda pggb
conda activate pggb_env


