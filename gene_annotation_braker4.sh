###############################################################################
#                    BRAKER4 GENE ANNOTATION PIPELINE
###############################################################################
#
# SPECIES:              Varanus komodoensis
# HPC OS:               Rocky Linux 8.10
# JOB SCHEDULER:        SLURM 22.05.9
# CONTAINER RUNTIME:    Apptainer 1.5.5
# WORKFLOW MANAGER:     Snakemake 8.18.2
# SLURM PLUGIN:         snakemake-executor-plugin-slurm 2.6.0
#
# ANNOTATION APPROACH:
#   BRAKER4 Evidence-based GeneMark-ETP (ETP)
#
# EVIDENCE:
#   1. Varanus komodoensis genome assembly
#   2. Repeat-masked genome
#   3. OrthoDB v12 Vertebrata protein database
#   4. Pogona vitticeps protein FASTA
#   5. Eublepharis macularius protein FASTA
#   6. Alligator mississippiensis protein FASTA
#   7. Paired-end RNA-seq reads (SRR8735152)
#
# RESOURCE SETTINGS:
#   CPUs:             26
#   SLURM partition:  compute
#   Default memory:   120000 MB
#   Maximum runtime:  4320 minutes
#
# HOW TO USE THIS DOCUMENT:
#
#   This is a MANUAL, COPY-PASTE-READY LABORATORY PROTOCOL.
#
#   Do not execute the entire document as one script.
#   Read each section and copy-paste the commands into the terminal.
#
#   [ONCE]          Run only during initial software/database setup.
#   [PER GENOME]    Repeat when annotating another genome.
#   [PER DATASET]   Repeat when introducing new evidence data.
#   [EACH SESSION]  Repeat when opening a new terminal session.
#   [PER RUN]       Execute for each new annotation workflow.
#   [MONITOR]       Execute whenever required during annotation.
#
###############################################################################



###############################################################################
# STEP 1. CHECK HPC ENVIRONMENT
#
# FREQUENCY: Once per HPC account; repeat when moving to another cluster.
#
# Identify the operating system, scheduler and available compute resources.
###############################################################################

# Check operating system
cat /etc/os-release

# Check SLURM version
sinfo --version

# Display available partitions
sinfo

# Display partition, availability, time limit, nodes, CPUs and memory
sinfo -o "%P %a %l %D %c %m"

# Check Python version
python3 --version

# Check Git installation
git --version



###############################################################################
# STEP 2. INSTALL APPTAINER
#
# FREQUENCY: Once per HPC account.
#
# Apptainer is the container runtime used by BRAKER4.
# It is the successor to SingularityCE.
#
# You do not need to install both Singularity and Apptainer.
#
# Your existing installation:
#   ~/software/apptainer/bin
#
# IMPORTANT:
# The following installation commands are for a NEW HPC account.
# Your current Apptainer installation is already completed.
#
# Source compilation requires Go, GCC, make and the required system
# development libraries. Consult the official documentation and HPC
# administrator if dependencies or privileges are missing.
###############################################################################

# Create software directory
mkdir -p ~/software

# Enter software directory
cd ~/software

# Download Apptainer source version 1.5.5
wget https://github.com/apptainer/apptainer/releases/download/v1.5.5/apptainer-1.5.5.tar.gz

# Extract source archive
tar -xzf apptainer-1.5.5.tar.gz

# Enter source directory
cd apptainer-1.5.5

# Configure user-space installation
./mconfig --prefix="$HOME/software/apptainer"

# Compile Apptainer
make -C builddir

# Install Apptainer
make -C builddir install

# Add Apptainer to PATH permanently
echo 'export PATH="$HOME/software/apptainer/bin:$PATH"' >> ~/.bashrc

# Reload shell configuration
source ~/.bashrc

# Verify installation
apptainer --version

# Expected version from this workflow:
# apptainer version 1.5.5-1.el8



###############################################################################
# STEP 3. CONFIGURE APPTAINER
#
# FREQUENCY:
#   Installation: Once.
#   Environment variables: Every new terminal session.
###############################################################################

# Add Apptainer to PATH
export PATH="$HOME/software/apptainer/bin:$PATH"

# Check Apptainer version
apptainer --version

# Create dedicated Apptainer cache directory
mkdir -p /home/pranav/apptainer_cache

# Set cache directory
export APPTAINER_CACHEDIR=/home/pranav/apptainer_cache

# Test cache write permissions
touch "$APPTAINER_CACHEDIR/test_file"

# Remove test file
rm "$APPTAINER_CACHEDIR/test_file"

# Test container execution
apptainer exec docker://alpine:latest cat /etc/os-release

# Add cache variable permanently to bash configuration
echo 'export APPTAINER_CACHEDIR="$HOME/apptainer_cache"' >> ~/.bashrc

# Reload shell configuration
source ~/.bashrc



###############################################################################
# STEP 4. INSTALL SNAKEMAKE
#
# FREQUENCY: Once per HPC account.
#
# Snakemake manages the BRAKER4 workflow and submits individual jobs to SLURM.
#
# Versions used:
#   Snakemake: 8.18.2
#   SLURM executor plugin: 2.6.0
###############################################################################

# Create software directory
mkdir -p ~/software

# Enter software directory
cd ~/software

# Create Python virtual environment
python3 -m venv snakemake_env

# Activate virtual environment
source ~/software/snakemake_env/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install Snakemake
pip install snakemake==8.18.2

# Install pandas
pip install pandas

# Install SLURM executor plugin
pip install snakemake-executor-plugin-slurm==2.6.0

# Verify Snakemake
snakemake --version

# Verify SLURM plugin
pip show snakemake-executor-plugin-slurm

# Verify pandas
pip show pandas

# Expected Snakemake version:
# 8.18.2


# FOR FUTURE SESSIONS:
# Only activate the existing environment.
#
# source ~/software/snakemake_env/bin/activate



###############################################################################
# STEP 5. DOWNLOAD BRAKER4
#
# FREQUENCY: Once per HPC account.
#
# BRAKER4 repository contains:
#   - Main Snakefile
#   - Workflow rules
#   - Configuration and supporting files
#
# The bioinformatics tools used by BRAKER4 are executed inside containers.
###############################################################################

# Enter software directory
cd ~/software

# Clone official BRAKER4 repository
git clone https://github.com/Gaius-Augustus/BRAKER4.git

# Enter BRAKER4 directory
cd ~/software/BRAKER4

# Inspect repository
ls

# Verify main Snakefile
ls -lh Snakefile

# Record repository commit
git log -1 --oneline

# Record repository URL
git remote -v

# IMPORTANT:
# Do not clone BRAKER4 again when starting another genome annotation.
# Reuse the existing repository.
#
# To update BRAKER4, first check compatibility with your existing workflow
# and then deliberately update the repository.



###############################################################################
# STEP 6. DOWNLOAD ORTHODB v12 DATABASE
#
# FREQUENCY: Once per shared database.
#
# Database:
#   OrthoDB v12 Vertebrata protein FASTA
#
# Purpose:
#   Provides orthologous protein evidence for gene prediction.
#
# IMPORTANT:
# Download the Vertebrata partition from the official OrthoDB v12
# distribution. Verify the release and file integrity before use.
#
# The exact download endpoint and checksum should be taken from the
# official database distribution page; do not assume an old URL or
# checksum is valid for every mirror.
###############################################################################

# Create OrthoDB directory
mkdir -p ~/genome_assemblies/protein_evidence/OrthoDB_v12

# Enter directory
cd ~/genome_assemblies/protein_evidence/OrthoDB_v12

# Download Vertebrata.fa.gz from the official OrthoDB v12 distribution.
# If using a direct download URL supplied by OrthoDB, use:
#
# wget -c -O Vertebrata.fa.gz "OFFICIAL_DOWNLOAD_URL"
#
# Verify the downloaded archive using the checksum published by OrthoDB.
#
# Example checksum verification command:
#
# echo "PUBLISHED_MD5_CHECKSUM  Vertebrata.fa.gz" | md5sum -c -

# Extract the database while retaining the compressed archive
gunzip -c Vertebrata.fa.gz > Vertebrata.fa

# Check extracted file size
ls -lh Vertebrata.fa

# Count protein records
grep -c "^>" Vertebrata.fa

# Inspect FASTA headers
head -n 10 Vertebrata.fa

# IMPORTANT:
# Do not download or extract this database again for another genome
# if the same OrthoDB version is being used.



###############################################################################
# STEP 7. PREPARE ADDITIONAL REPTILE PROTEIN EVIDENCE
#
# FREQUENCY: Once per evidence dataset.
#
# Additional protein evidence:
#   1. Pogona vitticeps
#   2. Eublepharis macularius
#   3. Alligator mississippiensis
#
# These protein FASTAs are additional evidence beyond OrthoDB.
#
# Their original accessions and download URLs must be recorded by the user.
###############################################################################

# Create reptile protein directory
mkdir -p ~/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins

# Enter directory
cd ~/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins

# Place the following files in this directory:
#
# P_vitticeps_protein.faa
# E_macularius_protein.faa
# A_mississippiensis_protein.faa

# Check file sizes
ls -lh *.faa

# Count Pogona protein records
grep -c "^>" P_vitticeps_protein.faa

# Count Eublepharis protein records
grep -c "^>" E_macularius_protein.faa

# Count Alligator protein records
grep -c "^>" A_mississippiensis_protein.faa

# IMPORTANT:
# Record species, accession, source database and release version
# for each protein FASTA.



###############################################################################
# STEP 8. PREPARE GENOME ASSEMBLY
#
# FREQUENCY: Repeat for each new genome.
#
# Input:
#   Varanus komodoensis genome assembly
#
# Genome:
#   V_komo_genome.fasta
###############################################################################

# Define genome path
GENOME=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta

# Check file size
ls -lh "$GENOME"

# Count assembled sequences
grep -c "^>" "$GENOME"

# Inspect sequence headers
grep "^>" "$GENOME" | head

# IMPORTANT:
# Use the same assembly version throughout the annotation.
# Record assembly accession and version in the project README.



###############################################################################
# STEP 9. REPEAT MASKING
#
# FREQUENCY: Once per genome assembly.
#
# Your current workflow uses a soft-masked genome prepared using RepeatMasker.
#
# The repeat-masking analysis is an upstream workflow and should be completed
# before BRAKER4.
#
# The exact RepeatMasker installation and execution commands are not included
# here because the original commands and database versions were not recorded.
###############################################################################

# Define soft-masked genome
MASKED_GENOME=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/repeat_analysis/RM_files/repeatmaskerout_softmaskedgenome/V_komo_genome.fasta.masked

# Check masked genome size
ls -lh "$MASKED_GENOME"

# Count sequences
grep -c "^>" "$MASKED_GENOME"

# IMPORTANT:
# The genome and masked genome must correspond to the same assembly.
# Do not use a masked genome generated from a different assembly version.



###############################################################################
# STEP 10. PREPARE RNA-SEQ EVIDENCE
#
# FREQUENCY: Repeat for each new RNA-seq dataset.
#
# RNA-seq accession:
#   SRR8735152
#
# Evidence type:
#   Paired-end FASTQ
#
# BRAKER4 will perform HISAT2 alignment internally.
###############################################################################

# Define forward reads
RNA_R1=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/transcriptome_data/SRR8735152_1.fastq

# Define reverse reads
RNA_R2=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/transcriptome_data/SRR8735152_2.fastq

# Check file sizes
ls -lh "$RNA_R1" "$RNA_R2"

# Count forward reads
zcat -f "$RNA_R1" | awk 'END {print "R1 reads:", NR/4}'

# Count reverse reads
zcat -f "$RNA_R2" | awk 'END {print "R2 reads:", NR/4}'

# IMPORTANT:
# Paired-end files should have equal read counts.
# Verify read pairing before annotation.
#
# Do not supply the previous STAR BAM for this particular configuration.
# The sample sheet uses FASTQ evidence and BRAKER4 performs HISAT2 alignment.



###############################################################################
# STEP 11. CREATE BRAKER4 WORKING DIRECTORY
#
# FREQUENCY: Repeat for each genome.
#
# This directory stores:
#   - samples.csv
#   - config.ini
#   - AUGUSTUS configuration
#   - Snakemake metadata
#   - Container cache
#   - Annotation results
#   - SLURM logs
###############################################################################

# Create working directory
mkdir -p /home/pranav/genome_assemblies/braker_work

# Enter working directory
cd /home/pranav/genome_assemblies/braker_work

# Create additional logs directory
mkdir -p logs

# IMPORTANT:
# Use a separate working directory for each species/genome.
# Do not mix Snakemake metadata between independent annotation projects.



###############################################################################
# STEP 12. CREATE samples.csv
#
# FREQUENCY:
#   Create once per genome.
#   Do not recreate when resuming the same workflow.
#
# The sample sheet contains 13 columns.
#
# protein_fasta:
#   Multiple protein FASTAs are separated using colons.
#
# bam_files:
#   Empty because FASTQ input is used.
#
# fastq_r1 / fastq_r2:
#   Paired-end RNA-seq reads.
#
# busco_lineage:
#   sauropsida_odb12
###############################################################################

# Enter working directory
cd /home/pranav/genome_assemblies/braker_work

# Create sample sheet
cat > samples.csv << 'EOF'
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage
V_komodoensis,/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta,/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/repeat_analysis/RM_files/repeatmaskerout_softmaskedgenome/V_komo_genome.fasta.masked,/home/pranav/genome_assemblies/protein_evidence/OrthoDB_v12/Vertebrata.fa:/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins/P_vitticeps_protein.faa:/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins/E_macularius_protein.faa:/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins/A_mississippiensis_protein.faa,,/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/transcriptome_data/SRR8735152_1.fastq,/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/transcriptome_data/SRR8735152_2.fastq,,,,,,sauropsida_odb12
EOF

# Display sample sheet
column -s, -t samples.csv

# Check column count
awk -F',' '{print NF}' samples.csv

# Expected:
# 13 columns in header
# 13 columns in sample row

# IMPORTANT:
# Do not overwrite samples.csv when resuming an annotation.
# The file should only be created once for the genome.



###############################################################################
# STEP 13. CREATE config.ini
#
# FREQUENCY:
#   Create once per genome.
#   Do not recreate when resuming the same workflow.
#
# PARAMETER EXPLANATIONS:
#
# min_contig:
#   Minimum contig length considered by the workflow.
#
# use_compleasm_hints:
#   Enables compleasm-derived hints.
#
# skip_busco:
#   0 means BUSCO is not skipped.
#
# skip_optimize_augustus:
#   0 means AUGUSTUS optimization is enabled.
#
# cpus_per_task:
#   CPU count configured for BRAKER4 rules.
#
# mem_of_node:
#   Memory request in MB.
#
# max_runtime:
#   Maximum runtime in minutes.
###############################################################################

# Enter working directory
cd /home/pranav/genome_assemblies/braker_work

# Create configuration file
cat > config.ini << 'EOF'
[paths]
samples_file = samples.csv
augustus_config_path = augustus_config

[PARAMS]
fungus = 0
min_contig = 10000
use_varus = 0
skip_optimize_augustus = 0
skip_single_exon_downsampling = 0
use_dev_shm = 0
use_compleasm_hints = 1
skip_busco = 0
run_omark = 0
translation_table = 1
run_ncrna = 0
masking_tool = repeatmasker
no_cleanup = 0

[SLURM_ARGS]
cpus_per_task = 26
mem_of_node = 120000
max_runtime = 4320
EOF

# Verify configuration
grep -E '^(cpus_per_task|mem_of_node|max_runtime)' config.ini

# Expected:
# cpus_per_task = 26
# mem_of_node = 120000
# max_runtime = 4320

# IMPORTANT:
# Do not recreate config.ini when resuming a workflow.



###############################################################################
# STEP 14. VALIDATE ALL INPUT FILES
#
# FREQUENCY: Before every new annotation.
#
# Verify that all genome and evidence files exist and are non-empty.
###############################################################################

# Define genome and evidence paths
GENOME=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta

MASKED_GENOME=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/repeat_analysis/RM_files/repeatmaskerout_softmaskedgenome/V_komo_genome.fasta.masked

ORTHODB=/home/pranav/genome_assemblies/protein_evidence/OrthoDB_v12/Vertebrata.fa

POGONA=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins/P_vitticeps_protein.faa

EUBLEPHARIS=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins/E_macularius_protein.faa

ALLIGATOR=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins/A_mississippiensis_protein.faa

RNA_R1=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/transcriptome_data/SRR8735152_1.fastq

RNA_R2=/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/transcriptome_data/SRR8735152_2.fastq

# Validate all files
for file in "$GENOME" "$MASKED_GENOME" "$ORTHODB" \
            "$POGONA" "$EUBLEPHARIS" "$ALLIGATOR" \
            "$RNA_R1" "$RNA_R2"
do
    if [[ -s "$file" ]]; then
        echo "OK: $file"
    else
        echo "MISSING OR EMPTY: $file"
    fi
done

# IMPORTANT:
# Do not proceed if a required file is missing or empty.



###############################################################################
# STEP 15. TEST APPTAINER ON COMPUTE NODE
#
# FREQUENCY:
#   First use on HPC.
#   Repeat after changes to Apptainer or the cluster environment.
#
# This creates a short SLURM allocation and tests container execution
# on a compute node.
###############################################################################

srun -p compute -N 1 -n 1 --mem=1G -t 00:05:00 \
    apptainer exec docker://alpine:latest cat /etc/os-release



###############################################################################
# STEP 16. ACTIVATE SOFTWARE ENVIRONMENT
#
# FREQUENCY: Every new terminal session.
#
# Activate Snakemake and configure Apptainer before running BRAKER4.
###############################################################################

# Activate Snakemake environment
source ~/software/snakemake_env/bin/activate

# Add Apptainer to PATH
export PATH="$HOME/software/apptainer/bin:$PATH"

# Set Apptainer cache
export APPTAINER_CACHEDIR=/home/pranav/apptainer_cache

# Enter BRAKER4 working directory
cd /home/pranav/genome_assemblies/braker_work

# Check Snakemake
snakemake --version

# Check Apptainer
apptainer --version



###############################################################################
# STEP 17. RUN BRAKER4 DRY RUN
#
# FREQUENCY:
#   Before submitting each new workflow.
#   Also repeat after significant configuration changes.
#
# PURPOSE:
#   Construct the Snakemake dependency graph without submitting jobs.
#
# Your previous successful dry run identified 38 jobs.
###############################################################################

snakemake \
    --snakefile ~/software/BRAKER4/Snakefile \
    --directory /home/pranav/genome_assemblies/braker_work \
    --dry-run \
    --executor slurm \
    --default-resources slurm_partition=compute mem_mb=120000 \
    --cores 26 \
    --jobs 26 \
    --use-apptainer \
    --apptainer-prefix .apptainer_cache \
    --apptainer-args "-B /home -B /scratch" \
    --latency-wait 120 \
    --restart-times 3

# Review:
#   - Total number of jobs
#   - Rule dependencies
#   - CPU requirements
#   - Memory requirements
#   - Input and output paths
#
# A dry run does not submit actual annotation jobs.



###############################################################################
# STEP 18. SUBMIT BRAKER4 TO SLURM
#
# FREQUENCY: Once per annotation workflow.
#
# IMPORTANT:
# Run this section inside a Screen session.
# Do not submit another copy while the existing workflow is active.
###############################################################################

# Start a named Screen session
screen -S braker4_komodo

# Inside Screen, enter working directory
cd /home/pranav/genome_assemblies/braker_work

# Activate Snakemake
source ~/software/snakemake_env/bin/activate

# Configure Apptainer
export PATH="$HOME/software/apptainer/bin:$PATH"

export APPTAINER_CACHEDIR=/home/pranav/apptainer_cache

# Submit BRAKER4 workflow
snakemake \
    --snakefile ~/software/BRAKER4/Snakefile \
    --directory /home/pranav/genome_assemblies/braker_work \
    --executor slurm \
    --default-resources slurm_partition=compute mem_mb=120000 \
    --cores 26 \
    --jobs 26 \
    --use-apptainer \
    --apptainer-prefix .apptainer_cache \
    --apptainer-args "-B /home -B /scratch" \
    --latency-wait 120 \
    --restart-times 3

# IMPORTANT:
# The command above starts the Snakemake master process.
# Individual BRAKER4 rules are submitted as SLURM jobs.
#
# To detach from Screen:
#
# Press Ctrl+A
# Then press D
#
# This does not terminate Snakemake or the submitted SLURM jobs.



###############################################################################
# STEP 19. MONITOR SLURM JOBS
#
# FREQUENCY: Whenever required during annotation.
###############################################################################

# Display active and queued jobs
squeue -u pranav

# Display job ID, name, state, elapsed time and node
squeue -u pranav -o "%.18i %.20j %.8T %.10M %.10l %.6D %R"

# Monitor memory usage of a running job
# Replace JOBID with the actual SLURM job ID.
sstat -j JOBID.batch --format=JobID,AveRSS,MaxRSS,AveVMSize,MaxVMSize

# Inspect completed job statistics
sacct -j JOBID --format=JobID,JobName,State,Elapsed,ReqMem,MaxRSS

# Display node-level RAM
# This reports the node where the command is executed.
free -h



###############################################################################
# STEP 20. MONITOR SCREEN AND SNAKEMAKE
#
# FREQUENCY: Whenever required.
###############################################################################

# List Screen sessions
screen -ls

# Reattach to BRAKER4 Screen
screen -r braker4_komodo

# Check whether Snakemake master process is running
pgrep -af snakemake

# Detach from Screen without terminating the process:
#
# Ctrl+A, then D



###############################################################################
# STEP 21. INSPECT BRAKER4 LOGS
#
# FREQUENCY: During annotation and troubleshooting.
#
# SLURM logs are stored inside the Snakemake working directory.
###############################################################################

# Enter working directory
cd /home/pranav/genome_assemblies/braker_work

# List all SLURM logs
find .snakemake/slurm_logs -type f

# Inspect a particular job log
# Replace JOBID with the actual job ID.
#
# Example: prepare_genome job 7338
tail -n 30 .snakemake/slurm_logs/rule_prepare_genome/V_komodoensis/7338.log

# Follow a log in real time
tail -f .snakemake/slurm_logs/rule_prepare_genome/V_komodoensis/7338.log



###############################################################################
# STEP 22. RESUME AN INTERRUPTED WORKFLOW
#
# FREQUENCY: Only if Snakemake stops unexpectedly.
#
# Snakemake tracks completed outputs and workflow metadata.
# Re-running the same command in the same working directory allows it
# to determine which jobs remain unfinished.
#
# FIRST:
#   Check that the previous Snakemake master process has stopped.
#   Check SLURM jobs.
#   Investigate and resolve the original failure.
###############################################################################

# Check Snakemake process
pgrep -af snakemake

# Check active SLURM jobs
squeue -u pranav

# Return to the same working directory
cd /home/pranav/genome_assemblies/braker_work

# Activate environment
source ~/software/snakemake_env/bin/activate

# Configure Apptainer
export PATH="$HOME/software/apptainer/bin:$PATH"

export APPTAINER_CACHEDIR=/home/pranav/apptainer_cache

# Re-run the same Snakemake submission command from STEP 18.
#
# IMPORTANT:
# Do not delete:
#   - .snakemake/
#   - samples.csv
#   - config.ini
#   - Completed intermediate outputs
#
# Do not start a second Snakemake master process while the original
# process is still running.



###############################################################################
# STEP 23. INSPECT FINAL ANNOTATION OUTPUTS
#
# FREQUENCY: After workflow completion.
###############################################################################

# Enter working directory
cd /home/pranav/genome_assemblies/braker_work

# Locate GFF3, GTF and FASTA outputs
find . -type f \
    \( -name "*.gff3" -o -name "*.gtf" \
       -o -name "*.faa" -o -name "*.fna" \)

# Locate BUSCO and compleasm results
find . \( -iname "*busco*" -o -iname "*compleasm*" \)

# Check directory sizes
du -sh *

# IMPORTANT OUTPUTS TO RETAIN:
#
# 1. Gene annotation GFF3/GTF
# 2. Predicted protein FASTA
# 3. Predicted CDS FASTA
# 4. Predicted transcript FASTA
# 5. BUSCO completeness results
# 6. Compleasm results
# 7. Snakemake and SLURM logs



###############################################################################
# STEP 24. REPRODUCIBILITY RECORD
#
# FREQUENCY: For every genome annotation.
#
# Record software versions, input accessions and database versions.
###############################################################################

# Record operating system
cat /etc/os-release

# Record Apptainer version
apptainer --version

# Record Snakemake version
snakemake --version

# Record SLURM version
sinfo --version

# Record BRAKER4 Git commit
git -C ~/software/BRAKER4 log -1 --oneline

# Record OrthoDB file size
ls -lh ~/genome_assemblies/protein_evidence/OrthoDB_v12/Vertebrata.fa

# Record genome file size
ls -lh /home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta

# Record RNA-seq accessions
echo "RNA-seq accession: SRR8735152"

# Record BUSCO lineage
echo "BUSCO lineage: sauropsida_odb12"


###############################################################################
# END OF BRAKER4 GENE ANNOTATION PROTOCOL
###############################################################################


###############################################################################
# FREQUENCY SUMMARY
###############################################################################
#
# ONCE PER HPC ACCOUNT:
#   STEP 1  - HPC environment checks
#   STEP 2  - Apptainer installation
#   STEP 3  - Apptainer configuration
#   STEP 4  - Snakemake installation
#   STEP 5  - BRAKER4 download
#
# ONCE PER SHARED DATABASE:
#   STEP 6  - OrthoDB v12 download
#
# ONCE PER EVIDENCE DATASET:
#   STEP 7  - Additional reptile protein datasets
#
# REPEAT FOR EACH GENOME:
#   STEP 8  - Genome assembly
#   STEP 9  - Repeat masking
#   STEP 10 - RNA-seq evidence
#   STEP 11 - Working directory
#   STEP 12 - samples.csv
#   STEP 13 - config.ini
#
# BEFORE EACH ANNOTATION:
#   STEP 14 - Input validation
#   STEP 15 - Apptainer compute-node test when needed
#   STEP 16 - Environment activation
#   STEP 17 - Snakemake dry run
#   STEP 18 - BRAKER4 submission
#
# DURING ANNOTATION:
#   STEP 19 - SLURM monitoring
#   STEP 20 - Screen monitoring
#   STEP 21 - Log inspection
#
# WHEN INTERRUPTED:
#   STEP 22 - Resume workflow
#
# AFTER COMPLETION:
#   STEP 23 - Inspect annotation outputs
#   STEP 24 - Record reproducibility information
#
###############################################################################
