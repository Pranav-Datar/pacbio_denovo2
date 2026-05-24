conda activate minimap2_env

minimap2 -x asm10 -t 30 /home/pranav/NCBI_FCS/contaminants_combined.fa primary.fasta > assembly_vs_contaminants.paf
awk '($10 > 5000) && (($10/$2) > 0.3) && (($10/$11) > 0.9)' assembly_vs_contaminants.paf > suspicious_contaminants.paf

#$10 > 5000: atleast 5kb of matching sequence
#(($10/$2) > 0.3): fraction_of_scaffold_matching > 30%
#(($10/$11) > 0.9): %identity (matching bases divided by alignment length)
