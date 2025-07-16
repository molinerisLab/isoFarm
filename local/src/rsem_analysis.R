#! /usr/bin/env Rscript
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
df_file <- args[1]
meta_file <- args[2]
out1_file <- args[3]
out2_file <- args[4]

df <- read.delim(df_file)
meta <- read.delim(meta_file)

# all expressed transcripts
df_exp <- df %>% 
  mutate(TPM = as.numeric(TPM)) %>%
  left_join(meta, by = "sample") %>%
  group_by(group, transcript_id, gene_id, gene_symbol, gene_type, transcript_type) %>%
  summarise(mean_TPM = mean(TPM)) %>%
  filter(mean_TPM >= 1) %>%
  ungroup()
write.table(df_exp, gzfile(out1_file), quote = F, sep = "\t", row.names = F)

# most expressed isoform
df_top_isoform <- df_exp %>%
  group_by(gene_id, group) %>%
  slice_max(order_by = mean_TPM, n = 1, with_ties = FALSE) %>%
  ungroup()
write.table(df_top_isoform, gzfile(out2_file), quote = F, sep = "\t", row.names = F)
