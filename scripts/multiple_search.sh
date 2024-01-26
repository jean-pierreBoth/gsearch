#!/bin/bash

# Function to resolve full path from a given path
resolve_full_path() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS specifics
        echo "$(cd "$(dirname "$1")"; pwd -P)/$(gbasename "$1")"
    else
        # Other Unix-like systems
        echo "$(cd "$(dirname "$1")"; pwd -P)/$(basename "$1")"
    fi
}

# Function to get the number of processors
get_num_processors() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # For macOS
        command -v gnproc > /dev/null 2>&1 && echo $(gnproc) || echo 1
    else
        # For Linux and other Unix-like systems
        echo $(nproc)
    fi
}

# Function to check if gsearch is installed
check_gsearch_installed() {
    if ! command -v gsearch > /dev/null 2>&1; then
        echo "gsearch is not installed. Please install it to use this script."
        exit 1
    fi
}

# Check if gsearch is installed
check_gsearch_installed

# Check for the existence of reformat.sh
if ! command -v reformat.sh >/dev/null 2>&1; then
    echo "reformat.sh is not found. Please ensure it is installed and in your PATH."
    exit 1
fi

# Check if input folder, new_folder, and output file name are provided
if [[ "$1" == "" || "$2" == "" || "$3" == "" ]]; then
    echo "
    Usage: $0 <input_folder> <new_folder> <output_file_name>
    
    input_folder         Input folder contains database files from multiple_build.sh.
    new_folder           Folder contains files to search against input folder.
    output_file_name     Output file name to store search results. 

    
    " >&2
    exit 1
fi

# Resolve full paths
INPUT_FOLDER=$(resolve_full_path "$1")
NEW_FOLDER=$(resolve_full_path "$2")
COMBINED_FILE=$(resolve_full_path "$3")

# Clear or create the combined file
> "$COMBINED_FILE"

# Get the number of processors
NUM_PROCESSORS=$(get_num_processors)

# Store the header (assuming the first line of the first file is the header)
HEADER_LINE=""

# Iterate over each subfolder in the input folder
for folder in "$INPUT_FOLDER"/*/; do
    folder_name=$(basename "$folder")

    # Run the gsearch command
    gsearch --pio 1000 --nbthreads $NUM_PROCESSORS request -b "$folder" -r "$NEW_FOLDER" -n 50

    # Check if gsearch.neighbors.txt was created
    if [[ ! -f "./gsearch.neighbors.txt" ]]; then
        echo "gsearch did not produce gsearch.neighbors.txt for $folder_name"
        continue
    fi

    # Run the reformat.sh script
    reformat.sh 16 1 ./gsearch.neighbors.txt ./clean.txt

    # Store the header if it's not already stored
    if [[ -z "$HEADER_LINE" ]]; then
        HEADER_LINE=$(head -n 1 ./clean.txt)
    fi

    # Append the cleaned data (excluding the header) to the combined file
    tail -n +2 ./clean.txt >> "$COMBINED_FILE"

    # Remove the intermediate files
    rm ./clean.txt ./gsearch.neighbors.txt ./gsearch.matches
done

# Add the header back to the top of the sorted file
{
    echo "$HEADER_LINE";
    if [[ "$OSTYPE" == "darwin"* ]]; then
        gsort -k1,1n -k2,2n -g <(cat "$COMBINED_FILE")
    else
        sort -k1,1n -k2,2n -g <(cat "$COMBINED_FILE") 
    fi
} > "${COMBINED_FILE}_sorted.txt"

rm "$COMBINED_FILE"

echo "Processing complete. Sorted results are in ${COMBINED_FILE}_sorted.txt."
