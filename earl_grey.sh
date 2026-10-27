conda create -n earl_grey_env -c conda-forge -c bioconda earlgrey
conda activate earl_grey_env
earlGrey -g /home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta -s Varanus_komodoensis -o earlgrey -t 16 -d yes -v yes -e yes

#-s: species name, consistent species identifier
#-d yes: creates softmasked genome at the end
#-v yes: validates resume by checksun
#-e yes: Run HELIANO as an optional step to detect Helitrons
