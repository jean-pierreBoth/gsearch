#!/bin/bash

# Function to resolve full path from relative path
resolve_full_path() {
    echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
}

# Function to check if GNU Parallel is installed
check_parallel_installed() {
    if ! command -v parallel > /dev/null 2>&1; then
        echo "GNU Parallel is not installed. Please install it to use this script." >&2
        exit 1
    fi
}

# Check if GNU Parallel is installed
check_parallel_installed

# Check if sufficient arguments are provided or if help is requested
if [[ "$1" == "" || "$2" == "" || "$3" == "" ]]; then
    echo "
    Usage: $0 <source_folder> <destination_folder> <number_of_subfolders>

    source_folder             Input folder contains files. Supports files with extensions:
                              *.fna, *.fna.gz, *.fasta, *.fasta.gz, *.fa, *.fa.gz, *.faa, *.faa.gz.
    destination_folder        Output folder where files will be distributed into subfolders.
    number_of_subfolders      Number of subfolders to split files into.
    " >&2
    exit 1
fi

SOURCE_FOLDER=$(resolve_full_path "$1")
DESTINATION_FOLDER=$(resolve_full_path "$2")
NUM_FOLDERS="$3"

# Initialize the random number generator with a fixed seed for reproducibility
RANDOM_SEED=42 # You can change this number to any other integer value
RANDOM=$RANDOM_SEED

# Check if the source folder exists and is a directory
if [ ! -d "$SOURCE_FOLDER" ]; then
    echo "Error: Source folder '$SOURCE_FOLDER' does not exist or is not a directory." >&2
    exit 1
fi

# Create destination folder if it doesn't exist
mkdir -p "$DESTINATION_FOLDER"

# Create the specified number of sub-folders
for (( i=1; i<=NUM_FOLDERS; i++ )); do
    SUBFOLDER="$DESTINATION_FOLDER/folder_$i"
    mkdir -p "$SUBFOLDER"
done

# Export NUM_FOLDERS and DESTINATION_FOLDER so it's available in the environment for GNU Parallel
export NUM_FOLDERS
export DESTINATION_FOLDER

# Function to perform the copy operation
copy_file() {
    local file=$1
    local folder_num=$((RANDOM % NUM_FOLDERS + 1))
    local destination_subfolder="$DESTINATION_FOLDER/folder_$folder_num"
    cp "$file" "$destination_subfolder/"
}

export -f copy_file

# Find files with specific extensions and distribute them using GNU Parallel
find "$SOURCE_FOLDER" -type f \( \
    -name "*.fna" -o \
    -name "*.fna.gz" -o \
    -name "*.fasta" -o \
    -name "*.fasta.gz" -o \
    -name "*.fa" -o \
    -name "*.fa.gz" -o \
    -name "*.faa" -o \
    -name "*.faa.gz" \
\) | parallel copy_file

echo "Files have been distributed into $NUM_FOLDERS sub-folders in $DESTINATION_FOLDER."
