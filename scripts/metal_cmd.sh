#!/usr/bin/env bash

# Usage:
#   metal_cmd.sh <command_file> <outfile_prefix> <input1> [input2 ...]
#
# Example:
#   metal_cmd.sh \
#       results/meta_analysis/BMI_metal_commands.txt \
#       results/BMI_ewas_meta_analysis_results_ \
#       results/F/F_BMI_ewas_bacon_results.csv.gz \
#       results/M/M_BMI_ewas_bacon_results.csv.gz

if [[ "$#" -lt 3 ]]; then
    echo "Usage: $0 <command_file> <outfile_prefix> <input1> [input2 ...]" >&2
    exit 2
fi

# First argument: file to create containing METAL commands.
output_script="$1"
shift

# Second argument: METAL output prefix.
# METAL appends its numbered output suffix, e.g. <prefix>1.txt.
outfile_prefix="$1"
shift

# All remaining arguments are bacon-adjusted EWAS result files.
input_files=( "$@" )

# Ensure the generated-command-file directory exists.
mkdir -p "$(dirname "$output_script")"

{
    echo "SCHEME STDERR"
    echo "COLUMNCOUNTING LENIENT"
    echo "SEPARATOR COMMA"
    echo "PVALUELABEL bacon.pval"
    echo "EFFECTLABEL bacon.es"
    echo "STDERRLABEL bacon.se"
    echo "MARKER cpgid"
    echo "WEIGHT n"
    echo

    for input_file in "${input_files[@]}"; do
        echo "PROCESSFILE $input_file"
    done

    echo
    echo "OUTFILE $outfile_prefix .txt"
    echo "ANALYZE HETEROGENEITY"
    echo "CLEAR"
} > "$output_script"

echo "METAL command file created: $output_script"