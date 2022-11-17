####################
#
# Plot Isoforms
#
#

# 0. Resources ----
library(RNAseqRtools)

expression.unit <- "TPM"
group <- "condition"
group_ref <- "Control"
shape_by <- NULL

out_dir <- "/sto1/epigen/GambarottaRonchi_GutPNS/RNAseq/dataset/v3_211209_DRG/isoforms/"
path_counts   <- "/sto1/epigen/GambarottaRonchi_GutPNS/RNAseq/dataset/v3_211209_DRG/isoforms/ALL_quant.TPM_matrix.Nrg1_isoforms.tsv"
path_metadata <- "/sto1/epigen/GambarottaRonchi_GutPNS/RNAseq/dataset/v3_211209_DRG/metadata.txt"
path_transcript_info <- "/sto1/epigen/GambarottaRonchi_GutPNS/RNAseq/dataset/v3_211209_DRG/Nrg1_isoforms.tsv"

counts <- read.delim(path_counts, row.names = 1)
metadata  <- read.delim(path_metadata)
metadata$condition <- factor(metadata$condition)
metadata$condition <- relevel(metadata$condition, ref=group_ref)

transcript_info <- read.delim(path_transcript_info)


rownames(transcript_info) <- transcript_info$transcript_ID
#transcript_info <- transcript_info[rownames(counts),]


# Create DGElist obj
y <- RNAseqRtools::processRNAseqEdgeR(m = counts
                        , experimental_info = metadata
                        , gene_info = transcript_info
                        , group = group
                        , reference = group_ref
                        , filter.expr.th = 1
                        , filter.sample.th = 1
                        )

# Set color palette
pal_1 <- c("grey","lightskyblue","deepskyblue4")

# Create list of plots
plot_list <- list()
for (transcript in row.names(y[["counts"]])) {
  p <- RNAseqRtools::plotExpression(y
                               , expression.unit = "counts"
                               , gene = transcript
                               , experimental_info = y$samples
                               , group.by = group
                               , pal = pal_1
                               , error.type = "se"
                               , point.height = 0
                               , show.names = T) + ylab("TPM")
  plot_list[[transcript]] <- p
}

# Save results
outfile <- paste0(out_dir, "/transcript_expression_Nrg1_isoforms_tpm.pdf")
pdf(file = outfile, paper = "a4", useDingbats = F)
ggpubr::ggarrange(plotlist=plot_list, nrow = 4, ncol = 4, common.legend = T)
dev.off()





