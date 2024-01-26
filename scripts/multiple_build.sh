#!/bin/bash

# Function to get the number of processors
get_num_processors() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # For macOS
        if command -v gnproc > /dev/null 2>&1; then
            echo $(gnproc)
        else
            echo "gnproc not installed. Please install it using Homebrew."
            exit 1
        fi
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

#!/bin/bash

# Check if input and destination folders are provided
if [[ "$1" == "" || "$2" == "" ]]; then
    echo "
    Usage: $0 <input_folder> <destination_folder>
    input_folder           Folder contains multiple folders with files in each folder
    destination_folder     Output folder prefix for output folder from gsearch for each subfolder in input_folder
    " >&2
    exit 1
fi

INPUT_FOLDER="$1"
DESTINATION_FOLDER="$2"
OUTPUT_FOLDER="${DESTINATION_FOLDER}_output"

# Create the output folder
mkdir -p "$OUTPUT_FOLDER"

# Get the absolute path of the input folder
INPUT_FOLDER_ABS=$(realpath "$INPUT_FOLDER")

# Get the number of processors
NUM_PROCESSORS=$(nproc)

# Iterate over each subfolder in the input folder
for subfolder in "$INPUT_FOLDER"/*/; do
    # Extract the subfolder name
    subfolder_name=$(basename "$subfolder")

    # Create a corresponding output subfolder
    output_subfolder="$OUTPUT_FOLDER/$subfolder_name"
    mkdir -p "$output_subfolder"

    # Change to the output subfolder
    pushd "$output_subfolder" > /dev/null

    # Construct the relative path from the output subfolder back to the input subfolder
    rel_path_to_input_subfolder="$(realpath --relative-to="$PWD" "$INPUT_FOLDER_ABS/$subfolder_name")"

    # Run the gsearch command
    gsearch --pio 4000 --nbthreads $NUM_PROCESSORS tohnsw -d "$rel_path_to_input_subfolder" -k 16 -s 15000 -n 128 --ef 1600 --algo optdens

    # Return to the original directory
    popd > /dev/null
done

echo "Processing complete. Output available in $OUTPUT_FOLDER."
