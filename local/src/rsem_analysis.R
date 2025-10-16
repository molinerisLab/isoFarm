#! /usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(rlang)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    "Usage: rsem_analysis.R <rsem_isoforms_tsv(.gz)> <metadata_tsv> <group_col> <out_expressed.gz> <out_top_iso.gz>\n",
    "Example: rsem_analysis.R rsem.isoforms.results.gene_symbol.gz metadata.txt genotype_general_time out1.gz out2.gz"
  )
}

rsem_file  <- args[1]
meta_file  <- args[2]
group_col  <- args[3]
out1_file  <- args[4]
out2_file  <- args[5]

# ---- Read input ----
df <- read.delim(rsem_file, check.names = FALSE)
meta <- read.delim(meta_file, check.names = FALSE)

# Sanity check
sample_col_df <- names(df)[1] 

if (sample_col_df %in% names(meta)) {
  sample_col_meta <- sample_col_df
} else if ("sample" %in% names(meta)) {
  meta <- meta %>% rename(!!sample_col_df := "sample")
  sample_col_meta <- sample_col_df
} else if (identical(names(meta)[1], "sample")) {
  meta <- meta %>% rename(!!sample_col_df := all_of("sample"))
  sample_col_meta <- sample_col_df
} else {
  stop(
    "Could not find a sample key to join on.\n",
    "- RSEM first column is '", sample_col_df, "'.\n",
    "- Metadata must contain a column named '", sample_col_df, "' OR a column named 'sample'.\n",
    "Please rename your metadata accordingly."
  )
}

if (!(group_col %in% names(meta))) {
  stop(
    "Grouping column '", group_col, "' not found in metadata. Available columns: ",
    paste(names(meta), collapse = ", ")
  )
}

missing_samples <- setdiff(unique(df[[sample_col_df]]), unique(meta[[sample_col_meta]]))
if (length(missing_samples) > 0) {
  stop(
    "The following samples in the RSEM table are missing in metadata: ",
    paste(missing_samples, collapse = ", ")
  )
}

# all expressed transcripts
df_exp <- df %>%
  mutate(TPM = suppressWarnings(as.numeric(TPM))) %>%
  left_join(meta, by = setNames(sample_col_meta, sample_col_df)) %>%
  filter(!is.na(.data[[group_col]])) %>%
  group_by(
    .data[[group_col]],
    transcript_id, gene_id, gene_symbol, gene_type, transcript_type
  ) %>%
  summarise(mean_TPM = mean(TPM), .groups = "drop_last") %>%
  filter(mean_TPM >= 1) %>%
  ungroup()
# save
write.table(df_exp, gzfile(out1_file), quote = FALSE, sep = "\t", row.names = FALSE)


# most expressed isoform
df_top_isoform <- df_exp %>%
  group_by(gene_id, .data[[group_col]]) %>%
  slice_max(order_by = mean_TPM, n = 1, with_ties = FALSE) %>%
  ungroup()
# save
write.table(df_top_isoform, gzfile(out2_file), quote = FALSE, sep = "\t", row.names = FALSE)
