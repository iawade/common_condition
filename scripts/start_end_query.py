import argparse
import sys
import time
import logging
import gzip

def get_gene_coordinates(ensembl_id, file='data/all_genes.tsv.gz', is_gzipped=True):

    # Dynamically create the output filename
    output_file = f"run_files/{ensembl_id}.bed"
    open_input = gzip.open if is_gzipped else open

    with open_input(file, 'rt') as infile, open(output_file, 'wt') as output:
        for line in infile:
            columns = line.rstrip('\n').split()
            if len(columns) >= 3 and columns[3] == ensembl_id:
                output.write(line)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Retrieve gene coordinates all_genes.tsv.gz file.')
    parser.add_argument('--gene_file', required=False, default='data/all_genes.tsv.gz', help='lookup file of all gene (start, end) pairs')
    parser.add_argument('--not_gzipped', action='store_true', help='is the lookup file gzipped? If not, use this flag')
    parser.add_argument('--ensembl_id', required=True, help='Ensembl gene ID')
    args = parser.parse_args()

    # Run the function with provided arguments
    get_gene_coordinates(args.ensembl_id, args.gene_file, not args.not_gzipped)