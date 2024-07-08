#!/bin/bash
# Author: Diana Proctor: code was cleaned up by ChatGPT
# Revised July 5, 2024

# Error handling
error_exit() {
    echo "$1" 1>&2
    return 1
}

############ Define paths to input
DATA_DIR="/data/proctordm/ESKAPE/data"
BASE_DIR="$DATA_DIR"
READS_DIR="/data/proctordm/ESKAPE/00_reads"

BINNING_DIR="$DATA_DIR/02_binning"
BINNING_SCRIPT="$BINNING_DIR/binning.sh"
LOG_DIR="$BINNING_DIR/logs"
FAILED_SAMPLES_LIST="$BINNING_DIR/failed_samples_list.txt"

# Define container path
CONTAINER_DIR="/data/proctordm/ESKAPE/container"
SIF_FILE="$CONTAINER_DIR/MAG_wf_containers_metawrap.sif"

# Define paths for megahit
MEGAHIT_FINAL_CONTIGS_LIST="$BASE_DIR/01_assembly_megahit/final_contigs_list.txt"
MEGAHIT_BINNING_DIR="$BINNING_DIR/02_megahit_bins"
MEGAHIT_SWARM_FILE="$BINNING_DIR/megahit_swarm_file"
MEGAHIT_FAILED_SAMPLES_LIST="$BINNING_DIR/megahit_failed_samples_list.txt"

# Define paths for spades output
SPADES_FINAL_CONTIGS_LIST="$BASE_DIR/01_assembly_spades/final_contigs_list.txt"
SPADES_BINNING_DIR="$BINNING_DIR/02_spades_bins"
SPADES_SWARM_FILE="$BINNING_DIR/spades_swarm_file"
SPADES_FAILED_SAMPLES_LIST="$BINNING_DIR/spades_failed_samples_list.txt"

############ Define job parameters
THREADS=32
MEMORY=128
TIME="24:00:00"

# Define metawrap parameters
MIN_CONTIG_LENGTH=5000
METAWRAP_THREADS=32
METAWRAP_MEMORY=128

# Ensure directories exist
mkdir -p "$BASE_DIR"
mkdir -p "$BINNING_DIR"
mkdir -p "$SPADES_BINNING_DIR"
mkdir -p "$MEGAHIT_BINNING_DIR"
mkdir -p "$LOG_DIR"

# Load Singularity module
echo "Loading Singularity module..."
module load singularity || error_exit "Error: Failed to load Singularity module"
if [ $? -ne 0 ]; then
    error_exit "Error: Singularity module did not load properly."
fi

# Create container directory if it doesn't exist
echo "Creating container directory at $CONTAINER_DIR..."
mkdir -p "$CONTAINER_DIR" || error_exit "Error: Failed to create container directory"

# Check if the Singularity image file exists
if [ -f "$SIF_FILE" ]; then
    echo "Singularity image file $SIF_FILE already exists. Skipping pull."
else
    echo "Pulling the metawrap container from Shub..."
    singularity pull "$SIF_FILE" shub://sskashaf/MAG_wf_containers:metawrap || error_exit "Error: Failed to pull metawrap container"
fi

# Change to base directory
echo "Changing to base directory at $BASE_DIR..."
cd "$BASE_DIR" || error_exit "Error: Failed to change to base directory"

# Ensure final contigs lists exist
if [ ! -f "$MEGAHIT_FINAL_CONTIGS_LIST" ]; then
    error_exit "Error: final contigs list for MegaHit does not exist at $MEGAHIT_FINAL_CONTIGS_LIST"
fi

if [ ! -f "$SPADES_FINAL_CONTIGS_LIST" ]; then
    error_exit "Error: final contigs list for SPAdes does not exist at $SPADES_FINAL_CONTIGS_LIST"
fi

# Read assemblies into arrays
mapfile -t MEGAHIT_ASSEMBLIES < "$MEGAHIT_FINAL_CONTIGS_LIST"
mapfile -t SPADES_ASSEMBLIES < "$SPADES_FINAL_CONTIGS_LIST"

# Debugging: Print the contents of the assemblies arrays
echo "Contents of MEGAHIT_ASSEMBLIES array:"
printf '%s\n' "${MEGAHIT_ASSEMBLIES[@]}"

echo "Contents of SPADES_ASSEMBLIES array:"
printf '%s\n' "${SPADES_ASSEMBLIES[@]}"

############ Create lists of forward and reverse reads
echo "Creating lists of forward and reverse reads..."
cd $BASE_DIR
ls -d $READS_DIR/*_1.fastq > READ1.list
ls -d $READS_DIR/*_2.fastq > READ2.list

# Check if READ1.list and READ2.list were created successfully
if [[ -f "READ1.list" && -f "READ2.list" ]]; then
    echo "READ1.list and READ2.list created successfully."
else
    error_exit "Failed to create READ1.list and/or READ2.list."
fi

# Read lists into arrays
echo "Reading lists into arrays..."
mapfile -t READ1 < READ1.list
mapfile -t READ2 < READ2.list

# Extract sample names from reads
read_samples=()
for read in "${READ1[@]}"; do
    read_samples+=($(basename "$read" | sed 's/_1.fastq//'))
done

# Extract sample names from assemblies
megahit_samples=()
for assembly in "${MEGAHIT_ASSEMBLIES[@]}"; do
    megahit_samples+=($(basename "$(dirname "$assembly")" | sed 's/_megahit_out//'))
done

spades_samples=()
for assembly in "${SPADES_ASSEMBLIES[@]}"; do
    spades_samples+=($(basename "$(dirname "$assembly")" | sed 's/_spades_out//'))
done

# Check for mismatches
mismatch_found=false
mismatch_details=()

for sample in "${read_samples[@]}"; do
    if [[ ! " ${megahit_samples[@]} " =~ " ${sample} " ]]; then
        mismatch_details+=("Sample $sample found in reads but not in MegaHit assemblies.")
        mismatch_found=true
    fi
    if [[ ! " ${spades_samples[@]} " =~ " ${sample} " ]]; then
        mismatch_details+=("Sample $sample found in reads but not in SPAdes assemblies.")
        mismatch_found=true
    fi
done

for sample in "${megahit_samples[@]}"; do
    if [[ ! " ${read_samples[@]} " =~ " ${sample} " ]]; then
        mismatch_details+=("Sample $sample found in MegaHit assemblies but not in reads.")
        mismatch_found=true
    fi
done

for sample in "${spades_samples[@]}"; do
    if [[ ! " ${read_samples[@]} " =~ " ${sample} " ]]; then
        mismatch_details+=("Sample $sample found in SPAdes assemblies but not in reads.")
        mismatch_found=true
    fi
done

if [ "$mismatch_found" = true ]; then
    echo "Warning: There are mismatches between the sample lists:"
    printf '%s\n' "${mismatch_details[@]}"
fi

# Create a single binning script
echo "Creating binning script at $BINNING_SCRIPT..."
cat << EOF > "$BINNING_SCRIPT"
#!/usr/bin/bash
ASSEMBLY=\$1
READ1=\$2
READ2=\$3
OUTPUT=\$4
module load singularity
source /usr/local/current/singularity/app_conf/sing_binds
singularity run $SIF_FILE \\
metawrap binning --metabat2 --maxbin2 --concoct -l $MIN_CONTIG_LENGTH -t $METAWRAP_THREADS -m $METAWRAP_MEMORY -a \${ASSEMBLY} -o \${OUTPUT} \${READ1} \${READ2}
EOF

# Check if the binning script was created successfully
if [[ -f "$BINNING_SCRIPT" ]]; then
    echo "Binning script created successfully at $BINNING_SCRIPT."
    chmod +x "$BINNING_SCRIPT" || error_exit "Error: Failed to make binning script executable"
else
    error_exit "Error: Failed to create binning script at $BINNING_SCRIPT"
fi

############ Initialize the swarm files
> "$MEGAHIT_SWARM_FILE"
> "$SPADES_SWARM_FILE"

############ Generate output directory names
echo "Generating output directory names..."
> "$BINNING_DIR/BINNING_OUT.list"
for read1 in "${READ1[@]}"; do
    out_dir=$(basename "$read1" | sed 's/_1.fastq/_binning_out/')
    echo "$out_dir" >> "$BINNING_DIR/BINNING_OUT.list"
done
mapfile -t BINNING_OUT < "$BINNING_DIR/BINNING_OUT.list"

# Debugging: Print the contents of BINNING_OUT array
echo "Contents of BINNING_OUT array:"
printf '%s\n' "${BINNING_OUT[@]}"

# Clean up temporary files
echo "Cleaning up temporary files..."
rm READ1.list READ2.list "$BINNING_DIR/BINNING_OUT.list"

# Check if arrays are populated
if [[ ${#READ1[@]} -eq 0 || ${#READ2[@]} -eq 0 || ${#BINNING_OUT[@]} -eq 0 ]]; then
    error_exit "Error: Arrays are not populated correctly."
else
    echo "Arrays properly populated."
fi

############ Create the Megahit swarm file
echo "Creating the Megahit swarm file..."
for (( i=0; i<${#READ1[@]}; i++ )); do
    read_base=$(basename "${READ1[i]}" "_1.fastq" | sed 's/_.*//')
    echo "Processing read base for Megahit: $read_base" >&2

    matching_assembly=""
    output_dir=""
    found_match=false

    for assembly in "${MEGAHIT_ASSEMBLIES[@]}"; do
        assembly_base=$(basename "$(dirname "$assembly")" | sed 's/_megahit_out//')
        echo "Comparing with Megahit assembly base: $assembly_base" >&2
        if [[ $assembly_base == "$read_base" ]]; then
            matching_assembly=$assembly
            output_dir="$MEGAHIT_BINNING_DIR/${read_base}_binning_out"
            echo "bash $BINNING_SCRIPT $matching_assembly ${READ1[i]} ${READ2[i]} $output_dir" >> "$MEGAHIT_SWARM_FILE"
            echo "Added to Megahit swarm: bash $BINNING_SCRIPT $matching_assembly ${READ1[i]} ${READ2[i]} $output_dir" >&2
            found_match=true
            break
        fi
    done

    if [[ "$found_match" = false ]]; then
        echo "No matching Megahit assembly found for read base: $read_base" >&2
    fi
done

# Check if Megahit swarm file was created successfully
if [[ -f "$MEGAHIT_SWARM_FILE" && -s "$MEGAHIT_SWARM_FILE" ]]; then
    echo "Megahit swarm file created at $MEGAHIT_SWARM_FILE"
else
    error_exit "Failed to create Megahit swarm file."
fi

############ Create the Spades swarm file
echo "Creating the Spades swarm file..."
for (( i=0; i<${#READ1[@]}; i++ )); do
    read_base=$(basename "${READ1[i]}" "_1.fastq" | sed 's/_.*//')
    echo "Processing read base for Spades: $read_base" >&2

    matching_assembly=""
    output_dir=""
    found_match=false

    for assembly in "${SPADES_ASSEMBLIES[@]}"; do
        assembly_base=$(basename "$(dirname "$assembly")" | sed 's/_spades_out//')
        echo "Comparing with Spades assembly base: $assembly_base" >&2
        if [[ $assembly_base == "$read_base" ]]; then
            matching_assembly=$assembly
            output_dir="$SPADES_BINNING_DIR/${read_base}_binning_out"
            echo "bash $BINNING_SCRIPT $matching_assembly ${READ1[i]} ${READ2[i]} $output_dir" >> "$SPADES_SWARM_FILE"
            echo "Added to Spades swarm: bash $BINNING_SCRIPT $matching_assembly ${READ1[i]} ${READ2[i]} $output_dir" >&2
            found_match=true
            break
        fi
    done

    if [[ "$found_match" = false ]]; then
        echo "No matching Spades assembly found for read base: $read_base" >&2
    fi
done

# Check if Spades swarm file was created successfully
if [[ -f "$SPADES_SWARM_FILE" && -s "$SPADES_SWARM_FILE" ]]; then
    echo "Spades swarm file created at $SPADES_SWARM_FILE"
else
    error_exit "Failed to create Spades swarm file."
fi

############ Submit the swarm jobs
echo "Submitting the Megahit swarm job..."
megahit_swarmid=$(swarm -f $MEGAHIT_SWARM_FILE --job-name megahit_binning -t $THREADS -g $MEMORY --time $TIME --logdir $LOG_DIR)

echo "Submitting the Spades swarm job..."
spades_swarmid=$(swarm -f $SPADES_SWARM_FILE --job-name spades_binning -t $THREADS -g $MEMORY --time $TIME --logdir $LOG_DIR)

# Check if swarm jobs were submitted successfully
if [[ -z "$megahit_swarmid" ]]; then
    error_exit "Failed to submit Megahit swarm jobs."
else
    echo "Megahit swarm jobs submitted with ID $megahit_swarmid"
fi

if [[ -z "$spades_swarmid" ]]; then
    error_exit "Failed to submit Spades swarm jobs."
else
    echo "Spades swarm jobs submitted with ID $spades_swarmid"
fi

########################### you may or may not be able to run all this on your sinteractive node
# Wait for swarm jobs to finish
echo "Waiting for Megahit swarm jobs to finish..."
while squeue -u $USER -j $megahit_swarmid > /dev/null 2>&1; do
    sleep 60
done

echo "Waiting for Spades swarm jobs to finish..."
while squeue -u $USER -j $spades_swarmid > /dev/null 2>&1; do
    sleep 60
done

# Write message indicating that swarm jobs finished successfully
echo "Swarm jobs completed successfully" >> "$LOG_DIR/swarm_completion.log"

# Analyze Megahit log files for failures and create a list of failed samples
echo "Analyzing Megahit log files for failures..."
> $MEGAHIT_FAILED_SAMPLES_LIST
for logfile in $LOG_DIR/*.o*; do
    if ! grep -q "PIPELINE SUCCESSFULLY FINISHED" "$logfile"; then
        sample=$(grep -oP '(?<=final.contigs.fa ).*(?= )' "$logfile" | head -n 1)
        echo "$sample" >> $MEGAHIT_FAILED_SAMPLES_LIST
    fi
done

# Check if Megahit failed samples list was created successfully and write "none" if empty
if [[ -s "$MEGAHIT_FAILED_SAMPLES_LIST" ]]; then
    echo "Megahit failed samples list created at $MEGAHIT_FAILED_SAMPLES_LIST"
else
    echo "none" > $MEGAHIT_FAILED_SAMPLES_LIST
fi

# Analyze Spades log files for failures and create a list of failed samples
echo "Analyzing Spades log files for failures..."
> $SPADES_FAILED_SAMPLES_LIST
for logfile in $LOG_DIR/*.o*; do
    if ! grep -q "PIPELINE SUCCESSFULLY FINISHED" "$logfile"; then
        sample=$(grep -oP '(?<=final.contigs.fa ).*(?= )' "$logfile" | head -n 1)
        echo "$sample" >> $SPADES_FAILED_SAMPLES_LIST
    fi
done

# Check if Spades failed samples list was created successfully and write "none" if empty
if [[ -s "$SPADES_FAILED_SAMPLES_LIST" ]]; then
    echo "Spades failed samples list created at $SPADES_FAILED_SAMPLES_LIST"
else
    echo "none" > $SPADES_FAILED_SAMPLES_LIST
fi

echo "Script completed successfully."
