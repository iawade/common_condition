library(data.table)
library(dplyr)

system("gsutil cp ${WORKSPACE_BUCKET}/snakemake_amr/variance/amr_snp_wise_pca_covariates_* .")
system("rm *null*")
# find the
for (file in dir()) {
	fwrite(
		fread(file) %>% filter(V2!="sparse"),
		sep = " ", quote=FALSE, col.names=FALSE,
		file = gsub("\\.varianceRatio\\.txt$", ".null.varianceRatio.txt", file)
	)
}

system("gsutil cp *null.varianceRatio.txt ${WORKSPACE_BUCKET}/snakemake_amr/variance/")