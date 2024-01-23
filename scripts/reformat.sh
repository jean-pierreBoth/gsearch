#!/bin/bash

# Function to check if a command exists
command_exists() {
    type "$1" &> /dev/null
}

# Determine the awk/gawk command based on the operating system
if [[ "$(uname)" == "Darwin" ]]; then
    # macOS system detected, check for gawk
    AWK="gawk"
    if ! command_exists $AWK; then
        echo "gawk is not installed. Please install gawk using 'brew install gawk'."
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
if [[ "$1" == "" || "$1" == "-h" || "$2" == "" ]]; then
    echo "
    Usage: ./reformat.sh gsearch.neighbor.txt output_tabular.txt

    input_file    File containing the data to be transformed into tabular format.
    output_file   File where the tabular output will be saved.
    " >&2
    exit 1
fi

# Assign command-line arguments to variables
input_file="$1"
output_file="$2"

# Process the file with awk/gawk
# You can use other evolutionary models such as Binomial model. Change the ANI line as you want: ANI = ((1-$4)*2/(1-$4+1))^(1/16)
$AWK -F '\t' '
BEGIN {
    print "Query ID\tDistance\tAnswer FASTA Path\tAnswer Seq Len\tANI"
}
/^query_id:/ {
    query_id = $2
    distance = $4
    answer_fasta_path = $6
    answer_seq_len = $8
    ANI = (1 + log((1-$4)*2/(1-$4+1))/16*100
    print query_id "\t" distance "\t" answer_fasta_path "\t" answer_seq_len "\t" ANI
}
' "$input_file" > "$output_file"
