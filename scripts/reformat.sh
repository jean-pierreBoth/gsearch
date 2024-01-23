#!/bin/bash

# Function to check if a command exists
command_exists() {
    type "$1" &> /dev/null
}

BASENAME="basename"
if [[ "$(uname)" == "Darwin" ]]; then
    # macOS system detected, check for gawk and gbasename
    AWK="gawk"
    BASENAME="gbasename"
    if ! command_exists $AWK; then
        echo "gawk is not installed. Please install gawk using 'brew install gawk'."
        exit 1
    fi
    if ! command_exists $BASENAME; then
        echo "gbasename is not installed. Please install coreutils using 'brew install coreutils'."
        exit 1
    fi
elif [[ "$(uname)" == "Linux" ]]; then
    # Linux system detected, check for awk
    AWK="awk"
    if ! command_exists $AWK; then
        echo "awk is not installed. Please install awk."
        exit 1
    fi
else
    echo "Unsupported operating system. This script supports macOS and Linux."
    exit 1
fi


# Check if correct number of arguments is provided
if [[ "$1" == "" || "$2" == "" || "$3" == "" || "$4" == "" || "$4" == "-h" ]]; then
    echo "
    Usage: ./script.sh kmer model input_file output_file

    kmer          The kmer value used for ANI calculation.
    model         The model to be used for ANI calculation (1 or 2), corresponding to Poisson model or Binomial model.
    input_file    File containing the data to be transformed into tabular format.
    output_file   File where the tabular output will be saved.
    " >&2
    exit 1
fi

# Assign command-line arguments to variables
kmer="$1"
model="$2"
input_file="$3"
output_file="$4"

# Process the file with awk/gawk
# You can use other evolutionary models such as Binomial model. Change the ANI line as you want: ANI = ((1-$4)*2/(1-$4+1))^(1/16)
$AWK -v basename_cmd="$BASENAME" -v kmer="$kmer" -v model="$model" -F '\t' '
function basename(path, cmd, _result) {
    cmd = cmd " \"" path "\""
    cmd | getline _result
    close(cmd)
    return _result
}
function calculate_ani(distance, kmer, model) {
    if (model == 1) {
        return (1 + log((1-distance)*2/(1-distance+1))/kmer)*100
    } else if (model == 2) {
        return ((1-distance)*2/(1-distance+1))^(1/kmer)*100
    } else {
        return "Invalid Model"
    }
}

BEGIN {
    print "Query_Name\tDistance\tNeighbor_Fasta_name\tNeighbor_Seq_Len\tANI"
}
/^query_id:/ {
    query_id = basename($2, basename_cmd)
    distance = $4
    answer_fasta_path = basename($6, basename_cmd)
    answer_seq_len = $8
    ANI = calculate_ani(distance, kmer, model)
    print query_id "\t" distance "\t" answer_fasta_path "\t" answer_seq_len "\t" ANI
}
' "$input_file" > "$output_file"