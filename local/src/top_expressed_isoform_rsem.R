#! /usr/bin/env Rscript
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

df <- read.delim(input_file, header = TRUE)
print(head(df))
df$TPM <- as.numeric(df$TPM)

top_tpm <- df %>%
  filter(TPM > 0) %>%
  group_by(gene_id) %>%
  filter(TPM == max(TPM, na.rm = TRUE)) %>%
  ungroup()

write.table(top_tpm, file = gzfile(output_file), sep = "\t", quote = FALSE, row.names = FALSE)
