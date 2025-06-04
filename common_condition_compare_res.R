# Load libraries
library(data.table)
library(ggplot2)

# Read files
sep23 <- fread("Stroke_SAIGE_combined.txt", sep = "\t", header = TRUE)
cond <- fread("Stroke_SAS_combined.txt", sep = "\t", header = FALSE)

# Assign proper column names to conditional file
setnames(cond, c(
  "Region", "Group", "max_MAF", "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
  "BETA_Burden", "SE_Burden", "Pvalue_cond", "Pvalue_Burden_cond",
  "Pvalue_SKAT_cond", "BETA_Burden_cond", "SE_Burden_cond", "MAC",
  "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare"
))

# Rename relevant p-value columns in cond to avoid collision
setnames(cond, old = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT"), 
              new = c("Pvalue_condfile", "Pvalue_Burden_condfile", "Pvalue_SKAT_condfile"))

# Merge on Region and Group
merged <- merge(sep23, cond, by = c("Region", "Group", "max_MAF"))

# Remove rows where Group contains "Cauchy"
merged_clean <- merged[!grepl("Cauchy", Group)]

# Write merged clean data to a file
fwrite(merged_clean, "merged_burden_pvalues_filtered.tsv", sep = "\t", quote = FALSE)

# Optional: Plot burden p-values
ggplot(merged_clean, aes(x = -log10(Pvalue_Burden_sep23), y = -log10(Pvalue_Burden_condfile))) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Pre-conditioning Burden P-value Comparison",
    x = "-log10(Pvalue_Burden) from Sep23",
    y = "-log10(Pvalue_Burden) from Conditional run"
  ) +
  theme_minimal()
