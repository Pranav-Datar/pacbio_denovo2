#PanSV naming
seqkit replace \
-p '^(\S+)' \
-r 'acanthurus#0#${1}' \
/home/pranav/genome_assemblies/NCBI_data/V_acanthurus_ncbi/GCA_050042745.1_ZJU_Vac_2.0_genomic.fasta \
> /home/pranav/genome_assemblies/pangenome_input/acanthurus.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'komodoensis#0#${1}' \
/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta \
> /home/pranav/genome_assemblies/pangenome_input/komodoensis.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'salvator#0#${1}' \
/home/pranav/genome_assemblies/NCBI_data/V_salvator_2_ncbi/genome/GWHHKDC00000000.1.genome.fasta \
> /home/pranav/genome_assemblies/pangenome_input/salvator.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'brevicauda#1#${1}' \
/home/pranav/genome_assemblies/primary_data/V_brevicauda/brevicauda_combined_reanalysis2/assembly_files/brevicauda.assembly.bp.hap1.p_ctg.fasta \
> /home/pranav/genome_assemblies/pangenome_input/brevicauda.hap1.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'brevicauda#2#${1}' \
/home/pranav/genome_assemblies/primary_data/V_brevicauda/brevicauda_combined_reanalysis2/assembly_files/brevicauda.assembly.bp.hap2.p_ctg.fasta \
> /home/pranav/genome_assemblies/pangenome_input/brevicauda.hap2.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'panoptes#1#${1}' \
/home/pranav/genome_assemblies/primary_data/V_panoptes/panoptes_reanalysis2/assembly_files/panoptes.assembly.bp.hap1.p_ctg.fasta \
> /home/pranav/genome_assemblies/pangenome_input/panoptes.hap1.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'panoptes#2#${1}' \
/home/pranav/genome_assemblies/primary_data/V_panoptes/panoptes_reanalysis2/assembly_files/panoptes.assembly.bp.hap2.p_ctg.fasta \
> /home/pranav/genome_assemblies/pangenome_input/panoptes.hap2.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'scalaris#1#${1}' \
/home/pranav/genome_assemblies/primary_data/V_scalaris/scalaris_reanalysis2/assembly_files/scalaris.assembly.bp.hap1.p_ctg.fasta \
> /home/pranav/genome_assemblies/pangenome_input/scalaris.hap1.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'scalaris#2#${1}' \
/home/pranav/genome_assemblies/primary_data/V_scalaris/scalaris_reanalysis2/assembly_files/scalaris.assembly.bp.hap2.p_ctg.fasta \
> /home/pranav/genome_assemblies/pangenome_input/scalaris.hap2.pansn.fa

seqkit replace \
-p '^(\S+)' \
-r 'scalaris#2#${1}' \
/home/pranav/genome_assemblies/primary_data/V_scalaris/scalaris_reanalysis2/assembly_files/scalaris.assembly.bp.hap2.p_ctg.fasta \
> scalaris.hap2.pansn.fa

#concatenation
cat /home/pranav/genome_assemblies/pangenome_input/*.pansn.fa \
> /home/pranav/genome_assemblies/pangenome_input/varanus.fa

#compression
bgzip -@ 16 /home/pranav/genome_assemblies/pangenome_input/varanus.fa

#indexing
conda activate samtools
samtools faidx varanus.fa.gz

#check if fragmented contigs are present,  and if present trim them off
seqkit stats /home/pranav/genome_assemblies/pangenome_input/*.pansn.fa

conda create -n pggb_env -c conda-forge -c bioconda pggb
conda activate pggb_env
mkdir -p /home/pranav/genome_assemblies/pggb_partition

partition-before-pggb \
    -i /home/pranav/genome_assemblies/pangenome_input/varanus.fa.gz \
    -o /home/pranav/genome_assemblies/pggb_partition \
    -t 32 \
    -p 90 \
    -s 10k
#-p 90: minimum average nucleotide identity for segments
#-s 10k: segmenting length for scaffolding the graph (starting seed)

mkdir -p commands scripts logs results

awk '/^pggb -i/{flag=1} flag{print}' \
varanus.fa.gz.c325321.11fba48.8625f7d.smooth.07-29-2026_23:48:35.params.yml \
> commands/all_commands.txt

perl -0pe 's/\\\n\s*/ /g' commands/all_commands.txt > commands/all_commands_oneline.txt

nano scripts/benchmark_large.slurm

#!/bin/bash
#SBATCH --job-name=pggb_large
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=300G
#SBATCH --output=logs/pggb_large_%j.out
#SBATCH --error=logs/pggb_large_%j.err

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pggb_env

cd /home/pranav/genome_assemblies/pggb_partition

pggb \
    -i varanus.fa.gz.c325321.community.1.fa \
    -o varanus.fa.gz.c325321.community.1.fa.out \
    -s 10000 \
    -l 50000 \
    -p 90 \
    -c 1 \
    -K 19 \
    -F 0.001 \
    -g 30 \
    -k 23 \
    -f 0 \
    -B 10M \
    -n 9 \
    -j 0 \
    -e 0 \
    -G 700,1100 \
    -O 0.001 \
    -d 100 \
    -Q Consensus_ \
    -Y "#" \
    --threads 64 \
    --poa-threads 64

#ctrl + o, enter, ctrl + x

sbatch scripts/benchmark_large.slurm

squeue -j #JOBID
scontrol show job #JOBID

odgi squeeze \
    -f gfa_list.txt \
    -o merged.og \
    -O \
    -t 32 \
    -P

    odgi sort \
-i merged.og \
-o results/merged.sorted.og \
-P \
-t 64

odgi layout \
-i results/merged.sorted.og \
-o results/merged.layout.og \
-t 64

odgi draw \
    -i results/merged.sorted.og \
    -c results/merged.layout.lay \
    -p results/merged.draw.png
