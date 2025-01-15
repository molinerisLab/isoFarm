#! /usr/bin/env Rscript
library(optparse)
library(dplyr)

option_list <- list(
  make_option(c("-m", "--multi_salmon"), action="store",
              help= "Path to multi_SalmonQuant.gz"),
  make_option(c("-a", "--annotation"), action="store",
              help= "Path to ensg 2 enst 2 gene_name [3 columns, no header]"),
  make_option(c("-a", "--annotation"), action="store",
              help= "Path to metadata.txt [first column named sample]"),
  make_option(c("-c", "--condition"), action="store",
              help= "Column name to be summarized."))
opt <- parse_args(OptionParser(option_list = option_list))


# read options ----
multi_salmon <- read.delim(opt$multi_salmon, header = T)
est_to_esg_to_genename <- read.delim(opt$annotation, header = F, col.names = c("ensg","gene_name","enst"))
# "/home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/37/primary_assembly.annotation.ensg2enst2gene_name.map"
meta <- read.delim(opt$metadata, header = T)

condition <- opt$condition

# merge data ----
multi_salmon$Name <- gsub("(ENST[^|]*)\\|.*", "\\1", multi_salmon$Name)
multi_salmon <- merge(multi_salmon, meta, by.x = "Sample", by.y = "sample")
multi_salmon$TPM <- as.numeric(multi_salmon$TPM)

multi_salmon_mean <- multi_salmon %>%
  group_by(.data[[condition]], Name) %>%
  summarise(mean_TPM = mean(TPM), Length) %>%
  ungroup()

multi_salmon_mean <- merge(multi_salmon_mean, est_to_esg_to_genename, by.x="Name", by.y="enst") 

multi_salmon_mean_gene <- multi_salmon_mean %>%
  group_by(.data[[condition]], gene_name) %>%
  filter(mean_TPM == max(mean_TPM)) %>%
  ungroup() %>%
  distinct()

outfile <- gsub(".gz$",".top_expressed_isoforms.gz",opt$multi_salmon)
write.table(multi_salmon_mean_gene, gzfile(outfile), sep = "\t", quote = F, row.names = F, col.names = T)
