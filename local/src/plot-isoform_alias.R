#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

# Read input files ----
tpm <- read.delim(snakemake@input[["tpm_matrix"]], header = T)
features <- read.delim(snakemake@input[["alias"]], header = T)

tpm_sub <- merge(tpm, features, by.x="Name", by.y="GENCODE")

# Barplot ----
pl <- ggplot(tpm_sub, aes(x=isoform, y=TPM, fill=isoform)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(.~Sample, nrow = 3) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        strip.background = element_blank(),
        strip.text = element_text(size = 12)) +
  geom_point(size = .5, show.legend = F) +
  xlab("") 

if(length(unique(tpm_sub$isoform))<=12){
  pl <- pl + scale_fill_manual(values = RColorBrewer::brewer.pal(n=length(unique(tpm_sub$isoform)), name="Set3"))
}

# Save ----
if(!dir.exists("Salmon_plot")) dir.create("Salmon_plot")
pdf(snakemake@output[["out"]], paper="a4r", width = 8, height = 7) 
print(pl)
dev.off()

