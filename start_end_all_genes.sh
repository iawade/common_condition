# Download the .gtf file
wget ftp://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz

# Extract the chromsome, start, and end positions
# Note: using gawk and gzip here as I'm on a mac, but awk and zcat work on linux

# Mac
gzcat Homo_sapiens.GRCh38.114.gtf.gz | gawk '$3 == "gene"' |
gawk -F'\t' '{
    match($9, /gene_id "([^"]+)"/, a);
    match($9, /gene_name "([^"]+)"/, b);
    print "chr" $1 "\t" $4 "\t" $5 "\t" a[1]
}' | gzip > data/all_genes.tsv.gz
