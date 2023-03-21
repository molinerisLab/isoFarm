#! /usr/bin/env Rscript

# 0. libraries ----
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(plyr))

# 1. create parser object ----
parser <- ArgumentParser(description = "Create lineplot from RNA counts data")

# add options
parser$add_argument("-c", "--count_matrix", action="store", required = T,
                    help = "feature count matrix path")
parser$add_argument("-m", "--metadata", action="store", required = T,
                    help = "metadata file path")
parser$add_argument("-f", "--feature", action="store", required = T,
                    help = "list of feature for count matrix subsetting")
parser$add_argument("-g", "--group_col", action="store", required = T,
                    help = "name of metadata column to use for averaging counts")
parser$add_argument("-d", "--condition_col", action="store",required = T,
                    help = "name of metadata column to use as line groups")
parser$add_argument("-t", "--time_point", action="store", required = T,
                    help = "name of metadata column to use as x axis in lineplot")
parser$add_argument("--time_order_col", action="store", default = NULL,
                    help = "name of metadata column to use for x axis ordering [Default: alphanumeric order]")
parser$add_argument("--pal", action="store", default=NULL,
                    help = "custom made palette [Default: RColorBrewer]")

# create args object
args <- parser$parse_args()


# 2. read input files ----
counts <- read.delim(args$count_matrix, header = T, row.names = 1)
features <- scan(arguments$enst_file, character(), quote = "")
meta <- read.delim(arguments$metadata, header = T)

# filter TPM matrix with enst of interest
counts <- counts[row.names(counts) %in% features,]

#group <- arguments$group_col
group_col <- "time_condition"

# average across group
toplot <- reshape2::melt(as.matrix(z), varnames = c("feature","sample"), value.name = "TPM")
toplot <- merge(toplot, meta, by = "sample")
toplot <- ddply(toplot, c("feature", group_col)
                , summarize
                , av = mean(TPM)
                , se = sd(TPM)/sqrt(length(TPM)), 
                .drop = F)

# add condition col
#condition_col <- "time"
#if(!is.null(opt$condition_col)){
#  opt$condition_col <- condition_col
#  toplot <- merge(toplot, meta[, c(group_col, condition_col)])
#  colnames(toplot)[condition_col] <- "cond"
#} else {
#  toplot$cond <- "none"
#}
timepoint <- "time"
condition_col <- "condition"
toplot <- merge(toplot, meta[,c(group_col,timepoint,condition_col)])

# set group order
if(!is.null(opt$group_order)){
  group_order <- order(meta[,opt$group_order])
  #group_order <- order(meta[,"day_number"])
  toplot[,group_col] <- factor(toplot[,group_col], levels = unique(meta[,group_col][group_order]))
} else {
  toplot[,group_col] <- factor(toplot[,group_col])
}

toplot$time <- factor(toplot$time, levels = c("H9","EBs_Day6","EBs_Day12"))

# palette ----
pal <- readRDS("../../../palette.rds")

##################
#
#  Line plot ----
#

p <- ggplot(toplot, aes_string(x=timepoint, y="av", color="feature", group="feature")) +
  geom_line() +
  geom_point() +
  facet_wrap(.~condition) +
  #geom_errorbar(aes(ymax = av + se, ymin = av - se, group=feature_col), width = 0.3, lwd = 0.25) +
  theme_pubr(base_size = 8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("average TPM") +
  scale_color_discrete(name="") +
  scale_y_continuous(trans='log2')

p <- ggplot(toplot, aes_string(x=timepoint, y="av", color=condition_col, group=condition_col)) +
  geom_line() +
  geom_point() +
  facet_wrap(.~feature, scales = "free_y") +
  #geom_errorbar(aes(ymax = av + se, ymin = av - se, group=feature_col), width = 0.3, lwd = 0.25) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("average TPM") +
  scale_color_discrete(type = pal[[condition_col]])


pdf("Salmon_plot/DNMT3B.enst.line.pdf", paper="a4", width = 5, height = 4.5) 
print(p)
dev.off()





                    

# Create and arrange list of plots
plot_list <- list()
for (feature_ID in row.names(rds[["CPM"]])) {
  p <- RNAseqRtools::plotExpression(rds
                               , expression.unit = "CPM"
                               , gene = feature_ID
                               , experimental_info = rds$samples
                               , group.by = group
                               , pal = pal
                               , error.type = "se"
                               , point.height = 0
                               , show.names = T) + 
    geom_col(col="black", linewidth=.1) + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  plot_list[[feature_ID]] <- p
}
plot_arrange <- ggpubr::ggarrange(plotlist=plot_list, common.legend = T)


# Save results ----
pdf(file = stdout(), paper = "a4", useDingbats = F)
print(plot_arrange)
dev.off()












