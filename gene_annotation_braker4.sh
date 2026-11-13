###############################################################################
#                    BRAKER4 GENE ANNOTATION PIPELINE
###############################################################################
#
# SPECIES:              Varanus komodoensis
# COMPUTE NODE:         node01
# HPC OS:               Rocky Linux 8.10
# EXECUTION MODE:       Local execution (no SLURM submission)
# CONTAINER RUNTIME:    Apptainer 1.5.5
# WORKFLOW MANAGER:     Snakemake 8.18.2
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
#   Memory limit:     120000 MB
#   Concurrent jobs:  1
#   Execution:        Local
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

# Check available CPU threads
nproc

# Check available memory
free -h

# Check disk space
df -h /home

# Check operating system
cat /etc/os-release

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
mkdir -p "$HOME/apptainer_cache"

# Set cache directory
export APPTAINER_CACHEDIR="$HOME/apptainer_cache"

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
# Snakemake manages the BRAKER4 workflow locally.
# Individual rules are executed directly on the compute node.
#
# Versions used:
#   Snakemake: 8.18.2
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

# Verify Snakemake
snakemake --version

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

# Create shared reptile protein evidence directory
mkdir -p ~/genome_assemblies/protein_evidence/reptile_proteins

cd ~/genome_assemblies/protein_evidence/reptile_proteins

# Place the protein FASTAs selected for this annotation here.
# Examples:
# V_komodoensis_protein.faa
# A_carolinensis_protein.faa
# N_naja_protein_decompressed.faa
# G_gallus_protein_decompressed.faa
# M_musculus_proteins_decompressed.faa
# H_sapiens_proteins_decompressed.faa

ls -lh

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
#   - Local Snakemake and workflow logs
###############################################################################

# Define species-specific working directory
SPECIES=V_komodoensis

WORKDIR=/home/pranav/genome_assemblies/NCBI_data/${SPECIES}_ncbi/braker_work

# Create and enter working directory
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# Create logs directory
mkdir -p logs

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
cd "$WORKDIR"

# Create sample sheet
cat > samples.csv << 'EOF'
sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage
V_komodoensis,/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/genome/V_komo_genome.fasta,/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/repeat_analysis/RM_files/repeatmaskerout_softmaskedgenome/V_komo_genome.fasta.masked,/home/pranav/genome_assemblies/protein_evidence/OrthoDB_v12/Vertebrata.fa:/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins/P_vitticeps_protein.faa:/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins/E_macularius_protein.faa:/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/reptile_proteins/A_mississippiensis_protein.faa,,/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/transcriptome_data/SRR8735152_1.fastq,/home/pranav/genome_assemblies/NCBI_data/V_komodoensis_ncbi/transcriptome_data/SRR8735152_2.fastq,,,,,,sauropsida_odb12
EOF


###############################################################################
# VALIDATE PROTEIN EVIDENCE
#
# Check protein FASTA formatting before combining databases.
###############################################################################

# Define protein evidence files
PROTEINS=(
    "$HOME/genome_assemblies/protein_evidence/OrthoDB_v12/Vertebrata.fa"
    "$HOME/genome_assemblies/protein_evidence/reptile_proteins/P_vitticeps_protein.faa"
    "$HOME/genome_assemblies/protein_evidence/reptile_proteins/E_macularius_protein.faa"
    "$HOME/genome_assemblies/protein_evidence/reptile_proteins/A_mississippiensis_protein.faa"
)

# Check file existence and non-zero size
for file in "${PROTEINS[@]}"; do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Missing or empty protein file: $file"
        exit 1
    fi

    echo "OK: $file"
done

# Check for non-standard protein sequence characters
python3 - "${PROTEINS[@]}" <<'PY'
import sys
import re

allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*-")

for filename in sys.argv[1:]:
    print(f"\nChecking: {filename}")

    records = 0
    errors = 0
    current_id = None
    has_sequence = False

    with open(filename) as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip():
                continue

            if line.startswith(">"):
                if current_id is not None and not has_sequence:
                    print(f"Empty sequence: {current_id}")
                    errors += 1

                current_id = line[1:].split()[0] if line[1:].split() else None
                has_sequence = False
                records += 1

            else:
                seq = re.sub(r"\s+", "", line).upper()

                if current_id is None:
                    print(f"Sequence before header at line {line_no}")
                    errors += 1
                    continue

                invalid = set(seq) - allowed

                if invalid:
                    print(f"Invalid characters at line {line_no}: {invalid}")
                    errors += 1

                has_sequence = True

    if current_id is not None and not has_sequence:
        print(f"Empty sequence: {current_id}")
        errors += 1

    print(f"Records: {records}")
    print(f"Errors: {errors}")

    if errors:
        sys.exit(1)
PY


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
cd "$WORKDIR"

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

# IMPORTANT:
# BRAKER4 internally reads these configuration keys.
# They are retained even when Snakemake runs locally.
#
# cpus_per_task:
#   CPU count used by BRAKER4 rules.
#
# mem_of_node:
#   Memory limit used in BRAKER4 rule resource definitions.
#
# max_runtime:
#   Runtime parameter retained for BRAKER4 configuration.
#   SLURM is not used by the local execution command.

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
# CHECK GENOME AND SOFTMASKED GENOME CONSISTENCY
###############################################################################

echo "===== Genome sequence count ====="
grep -c "^>" "$GENOME"

echo "===== Softmasked genome sequence count ====="
grep -c "^>" "$MASKED_GENOME"

echo "===== Genome and masked genome sizes ====="
du -h "$GENOME" "$MASKED_GENOME"

###############################################################################
# STEP 15. TEST APPTAINER
#
# FREQUENCY:
#   First use on a compute node.
#   Repeat after changes to Apptainer or the system environment.
#
# This tests container execution directly on the current node.
###############################################################################

# Test Apptainer directly on node01
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
cd "$WORKDIR"

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
# The successful V. komodoensis workflow completed 30 steps.
# The number of jobs may vary with input data and workflow configuration.
###############################################################################

snakemake \
    --snakefile /home/pranav/software/BRAKER4/Snakefile \
    --directory "$WORKDIR" \
    --dry-run \
    --cores 26 \
    --jobs 1 \
    --resources mem_mb=120000 \
    --default-resources "tmpdir='${WORKDIR}/tmp'" \
    --use-apptainer \
    --apptainer-prefix .apptainer_cache \
    --apptainer-args "-B /home" \
    --latency-wait 120 \
    --printshellcmds
# Review:
#   - Total number of jobs
#   - Rule dependencies
#   - CPU requirements
#   - Memory requirements
#   - Input and output paths
#
# A dry run does not submit actual annotation jobs.


###############################################################################
# CHECK GENEMARK-ETP PROTEIN DATABASE CONSISTENCY
#
# FREQUENCY:
#   Before resuming a workflow after changing protein evidence.
#   Before rerunning a failed GeneMark-ETP stage.
###############################################################################

# Define sample and GeneMark directories
SAMPLE=V_komodoensis

GM_DIR="output/${SAMPLE}/GeneMark-ETP"

# Check whether GeneMark configuration exists
if [[ -d "$GM_DIR" ]]; then

    GM_CONFIG=$(find "$GM_DIR" -maxdepth 1 -name "*.yaml" -print -quit)

    if [[ -n "$GM_CONFIG" ]]; then

        echo "===== GeneMark protein configuration ====="

        PROTEIN_DB=$(awk -F': ' '$1=="protdb_path" {
            print $2
        }' "$GM_CONFIG")

        echo "Configured database: $PROTEIN_DB"

        # Check configured database
        if [[ ! -s "$PROTEIN_DB" ]]; then
            echo "ERROR: Configured protein database missing or empty."
            exit 1
        fi

        # Check existing GeneMark intermediates
        for file in \
            "$GM_DIR/arx/proteins.fa" \
            "$GM_DIR/data/proteins.fa"
        do
            if [[ -e "$file" ]]; then
                echo "Existing intermediate: $file"
                ls -lh "$file"
            else
                echo "Not present: $file"
            fi
        done

        # Compare original input with GeneMark's copied input
        if [[ -s "$GM_DIR/arx/proteins.fa" ]]; then

            if cmp -s "$PROTEIN_DB" "$GM_DIR/arx/proteins.fa"; then
                echo "OK: GeneMark input matches configured database."
            else
                echo "WARNING: GeneMark input differs from configured database."
                echo "Investigate stale intermediate before resuming."
            fi

        fi

    else
        echo "No existing GeneMark YAML configuration."
        echo "This appears to be a first GeneMark run."
    fi

else
    echo "No existing GeneMark-ETP directory."
    echo "This appears to be a first annotation run."
fi.


###############################################################################
# IF THE PROTEIN DATABASE HAS CHANGED:
#
# Preserve and move stale GeneMark intermediates before rerunning.
# Execute only after confirming that the configured database is correct.
###############################################################################

# Example for V. komodoensis:
#
# mkdir -p "$GM_DIR/stale_backup"
#
# mv "$GM_DIR/arx/proteins.fa" \
#    "$GM_DIR/stale_backup/proteins.fa.previous"
#
# mv "$GM_DIR/data/proteins.fa" \
#    "$GM_DIR/stale_backup/proteins.fa.reformatted_previous"
#
# Do not execute these commands automatically during ordinary resumption.



# Submit BRAKER4 workflow
###############################################################################
# STEP 18. RUN BRAKER4 LOCALLY
#
# FREQUENCY: Once per annotation workflow.
#
# IMPORTANT:
# Run this section inside a Screen session.
# The workflow executes directly on node01.
# No SLURM jobs are submitted.
###############################################################################

###############################################################################
# STEP 18. RUN BRAKER4 LOCALLY
#
# FREQUENCY: Once per annotation workflow.
#
# EXECUTION MODE:
#   Local execution on node01 using Snakemake and Apptainer.
#   No SLURM jobs are submitted.
#
# IMPORTANT:
# Run this section inside a Screen session.
# Do not start another Snakemake instance while the workflow is active.
###############################################################################

# Start a named Screen session
screen -S braker_new

# Enter working directory
cd "$WORKDIR"

# Activate Snakemake environment
source ~/software/snakemake_env/bin/activate

# Configure Apptainer
export PATH="$HOME/software/apptainer/bin:$PATH"
export APPTAINER_CACHEDIR="$HOME/apptainer_cache"

# Run BRAKER4 locally
/home/pranav/software/snakemake_env/bin/snakemake \
    --snakefile /home/pranav/software/BRAKER4/Snakefile \
    --directory "$WORKDIR" \
    --cores 26 \
    --jobs 1 \
    --resources mem_mb=120000 \
    --default-resources "tmpdir='${WORKDIR}/tmp'" \
    --use-apptainer \
    --apptainer-prefix .apptainer_cache \
    --apptainer-args "-B /home" \
    --latency-wait 120 \
    --printshellcmds

# Detach from Screen without terminating the workflow:
# Press Ctrl+A, then D.
#
# IMPORTANT:
#   --cores 26: Maximum CPU threads available to Snakemake.
#   --jobs 1: Run one Snakemake job at a time.
#   --resources mem_mb=120000: Snakemake memory scheduling limit (120 GB).
#   --use-apptainer: Execute tools inside Apptainer containers.
#
# The workflow runs directly on node01.
# Individual rules are executed locally, not submitted to SLURM.
###############################################################################



###############################################################################
# STEP 19. MONITOR LOCAL COMPUTATION
#
# FREQUENCY: Whenever required during annotation.
###############################################################################

# Monitor CPU and memory consumption by your processes
ps -u "$USER" -o pid,ppid,%cpu,%mem,etime,stat,cmd --sort=-%cpu | head -20

# Check active BRAKER4 processes
ps -u "$USER" -o pid,%cpu,%mem,etime,cmd | \
    grep -E 'snakemake|hisat2|samtools|stringtie|gmetp|GeneMark|augustus|optimize_augustus' | \
    grep -v grep

# Check Snakemake master process
pgrep -afu "$USER" snakemake

# Check available memory
free -h

# Check disk space
df -h /home

# Monitor BRAKER4 output directory
watch -n 30 'du -sh output/V_komodoensis/'

# Monitor system load
uptime



###############################################################################
# STEP 20. MONITOR SCREEN AND SNAKEMAKE
#
# FREQUENCY: Whenever required.
###############################################################################

# List Screen sessions
screen -ls

# Reattach to BRAKER4 Screen
screen -r braker_new

# Check whether Snakemake master process is running
pgrep -afu "$USER" snakemake

# Detach without terminating the workflow:
# Press Ctrl+A, then D


###############################################################################
# STEP 21. INSPECT BRAKER4 LOGS
#
# FREQUENCY: During annotation and troubleshooting.
###############################################################################

# Enter working directory
cd /home/pranav/genome_assemblies/braker_work

# List Snakemake execution logs
ls -ltr .snakemake/log/

# Inspect the latest Snakemake log
LATEST_LOG=$(ls -t .snakemake/log/*.snakemake.log | head -1)
tail -n 50 "$LATEST_LOG"

# List BRAKER4 workflow logs
find logs/V_komodoensis -type f -printf '%TY-%Tm-%Td %TH:%TM %p\n' | sort

# Inspect GeneMark-ETP log, if available
tail -n 30 logs/V_komodoensis/genemark_etp/genemark_etp.log


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

###############################################################################
# STEP 22. RESUME AN INTERRUPTED WORKFLOW
#
# FREQUENCY: Only if Snakemake stops unexpectedly.
###############################################################################

# Check whether the previous Snakemake process has stopped
pgrep -afu "$USER" snakemake

# Check that no previous workflow instance is still running.
# Do not launch another instance while the original is active.

# Enter the same working directory
cd /home/pranav/genome_assemblies/braker_work

# Activate environment
source ~/software/snakemake_env/bin/activate

# Configure Apptainer
export PATH="$HOME/software/apptainer/bin:$PATH"
export APPTAINER_CACHEDIR=/home/pranav/apptainer_cache

# Resume the workflow locally
/home/pranav/software/snakemake_env/bin/snakemake \
    --snakefile /home/pranav/software/BRAKER4/Snakefile \
    --directory /home/pranav/genome_assemblies/braker_work \
    --cores 26 \
    --jobs 1 \
    --resources mem_mb=120000 \
    --default-resources "tmpdir='/home/pranav/genome_assemblies/braker_work/tmp'" \
    --use-apptainer \
    --apptainer-prefix .apptainer_cache \
    --apptainer-args "-B /home" \
    --latency-wait 120 \
    --printshellcmds

# IMPORTANT:
# Do not delete:
#   - .snakemake/
#   - samples.csv
#   - config.ini
#   - Completed intermediate outputs
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

# Inspect final annotation files
ls -lh output/V_komodoensis/results/

# Inspect quality-control outputs
ls -lh output/V_komodoensis/results/quality_control/

# Inspect completion marker
ls -lh output/V_komodoensis/results/.done

# Check directory size
du -sh output/V_komodoensis/results/

# IMPORTANT OUTPUTS TO RETAIN:
#
# 1. Gene annotation GFF3/GTF
# 2. Predicted protein FASTA
# 3. Predicted CDS FASTA
# 4. Predicted transcript FASTA
# 5. BUSCO completeness results
# 6. Compleasm results
# 7. Snakemake and BRAKER4 execution logs



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

# Record local execution configuration
echo "Execution mode: Local Snakemake"
echo "CPU limit: 26"
echo "Concurrent jobs: 1"
echo "Memory scheduling limit: 120000 MB"
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
# RECORD INPUT FILE CHECKSUMS
###############################################################################

# Genome and masked genome
sha256sum "$GENOME" "$MASKED_GENOME" \
    > input_genome_checksums.sha256

# Protein evidence
sha256sum "${PROTEINS[@]}" \
    > input_protein_checksums.sha256

# RNA-seq evidence
sha256sum "$RNA_R1" "$RNA_R2" \
    > input_rnaseq_checksums.sha256

# Record BRAKER4 repository commit
git -C ~/software/BRAKER4 rev-parse HEAD \
    > braker4_commit.txt
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
#   STEP 15 - Apptainer test when needed
#   STEP 16 - Environment activation
#   STEP 17 - Snakemake dry run
#   STEP 18 - Local BRAKER4 execution
#
# DURING ANNOTATION:
#   STEP 19 - Local computation monitoring
#   STEP 20 - Screen monitoring
#   STEP 21 - Log inspection
#
# WHEN INTERRUPTED:
#   STEP 22 - Resume local workflow
#
# AFTER COMPLETION:
#   STEP 23 - Inspect annotation outputs
#   STEP 24 - Record reproducibility information
#
###############################################################################
