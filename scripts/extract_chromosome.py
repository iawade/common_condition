import pandas as pd
import sys
import time
import gzip
import argparse

def get_gene_chr(ensembl_id, file='data/all_genes.tsv.gz', is_gzipped=True):
    # Dynamically create the output filename
    open_input = gzip.open if is_gzipped else open
    with open_input(file, 'rt') as infile:
        for line in infile:
            columns = line.rstrip('\n').split()
            if len(columns) >= 3 and columns[3] == ensembl_id:
                return columns[0]

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Retrieve chr for a gene.')
    parser.add_argument('--ensembl_id', required=True, help='Ensembl gene ID')
    args = parser.parse_args()
    # Run the function with provided arguments
    chr=get_gene_chr(args.ensembl_id)
    print(chr)