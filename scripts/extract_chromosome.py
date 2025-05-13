import pandas as pd
import sys
import time
import gzip

def get_gene_chr(ensembl_id, file='data/all_genes.tsv.gz', is_gzipped=True):
    # Dynamically create the output filename
    open_input = gzip.open if is_gzipped else open
    with open_input(file, 'rt') as infile:
        for line in infile:
            columns = line.rstrip('\n').split()
            if len(columns) >= 3 and columns[3] == ensembl_id:
                gene_line=line
    return gene_line[0]

file = "data/gene_phenotype_pairs_170425.csv"
df = pd.read_csv(file, sep=",")

for i in range(len(df.Region)):
	print(f'{df.phenotype[i]}, {get_gene_coordinates(df.Region[i])}')

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Retrieve chr for a gene.')
    parser.add_argument('--ensembl_id', required=True, help='Ensembl gene ID')
    args = parser.parse_args()
    # Run the function with provided arguments
    get_gene_chr(args.ensembl_id)