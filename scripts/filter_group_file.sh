#!/bin/bash

# Script to filter and format group files for a specific gene
# Usage: filter_group_file.sh <gene> <output_file> <group_file1> [group_file2] ...

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <gene> <output_file> <group_file1> [group_file2] ..." >&2
    exit 1
fi

gene="$1"
output="$2"
shift 2
group_files=("$@")

# Initialize output file
> "$output"

# Extract lines from input groups
for group in "${group_files[@]}"; do
    if [[ "$group" == *.gz ]]; then
        zcat "$group" | grep -m1 -A1 "$gene" >> "$output" || true
    else
        grep -m1 -A1 "$gene" "$group" >> "$output" || true
    fi
done

touch "$output"

# Validate number of lines
nlines=$(wc -l < "$output")
if [[ "$nlines" -ne 2 ]]; then
    echo "ERROR: Expected exactly 2 lines in $output, found $nlines" >&2
    exit 1
fi

# Read first line and second line
read -r first < "$output"
read -r second < <(tail -n +2 "$output")

# Split first line: first two columns, and the rest
set -- $first
gene_col=$1
col_var=$2
shift 2
first_var=$1
shift
rest="$*"

# Fix first variant with logging
fixed_first_var=$(echo "$first_var" | awk '{ 
    n = split($0,p,/[^[:alnum:]]+/); 
    if(n==4){ 
        chr=p[1]; pos=p[2]; ref=p[3]; alt=p[4]; 
        if(chr!~ /^chr/) chr="chr"chr;
        if(chr=="chr23") chr="chrX";
        if(chr=="chr24") chr="chrY";  
        printf("%s:%s:%s:%s\n",chr,pos,ref,alt) 
    } else { 
        print $0 
    } 
}')

if [[ "$first_var" != "$fixed_first_var" ]]; then
    echo "First variant reformatted: $first_var -> $fixed_first_var" >&2
    echo "Warning: please check that this reformatting makes sense. Note that if it is not, the result group-based result will be different" >&2
else
    echo "First variant format OK: $first_var" >&2
fi

# Fix rest of variants (silent)
fixed_rest=$(echo "$rest" | awk '{ 
    for(i=1;i<=NF;i++) { 
        n = split($i, p, /[^[:alnum:]]+/); 
        if(n==4) { 
            chr=p[1]; pos=p[2]; ref=p[3]; alt=p[4]; 
            if(chr !~ /^chr/) chr="chr" chr; 
            if(chr=="chr23") chr="chrX";
            if(chr=="chr24") chr="chrY";  
            printf "%s:%s:%s:%s", chr,pos,ref,alt 
        } else { 
            printf "%s", $i 
        } 
        if(i<NF) printf " " 
    } 
}')

# Write output
echo "$gene_col $col_var $fixed_first_var $fixed_rest" > "$output"
echo "$second" >> "$output"
