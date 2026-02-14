wget ftp://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz
gunzip Homo_sapiens.GRCh38.114.gtf.gz

gawk '$3 == "gene"' Homo_sapiens.GRCh38.114.gtf |
gawk -F'\t' '{
    match($9, /gene_id "([^"]+)"/, a);
    match($9, /gene_name "([^"]+)"/, b);
    print a[1] "\t" $1 "\t" $4 "\t" $5 "\t" b[1]
}' > all_genes.tsv

gzip all_genes.tsv


gawk -F'\t' '
$3=="gene" && ($9 ~ /gene_biotype "protein_coding"/ || $9 ~ /gene_type "protein_coding"/) {
    match($9, /gene_id "([^"]+)"/, a)
    match($9, /gene_name "([^"]+)"/, b)
    print a[1] "\t" $1 "\t" $4 "\t" $5 "\t" b[1]
}' Homo_sapiens.GRCh38.114.gtf > protein_coding_genes.tsv

gzip protein_coding_genes.tsv