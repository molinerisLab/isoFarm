#! /usr/bin/env Rscript

# 0. libraries ----
suppressPackageStartupMessages(library(argparse))
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
                    help = "custom made palette as RDS [Default: RColorBrewer]")
parser$add_argument("--outdir", action="store", default=".",
                    help = "output directory [Default: .]")
# create args object
args <- parser$parse_args()


# 2. assign arguments to obj  ----
if(!is.null(opt$outdir)) outdir <- opt$outdir else outdir <- "."
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
#counts <- read.delim(args$count_matrix, header = T, row.names = 1)
#features <- scan(arguments$enst_file, character(), quote = "")
#meta <- read.delim(arguments$metadata, header = T)
#group <- arguments$group_col

# files
counts <- read.delim("multi_SalmonQuant.TPM_matrix.gz", header = T, row.names = 1)
features <- scan("Salmon_plot/SOX5.enst", character(), quote = "")
meta <- read.delim("metadata.txt", header = T)

# other args
group_col <- opt$group_col
condition_col <- opt$condition_col

# convert all "-" to "."
meta <- as.data.frame(apply(meta, 2, function(x) gsub(x, pattern = "-", replacement = ".")))

# filter TPM matrix with enst of interest
counts <- counts[row.names(counts) %in% features,]

# average across group
toplot <- reshape2::melt(as.matrix(counts), varnames = c("feature","sample"), value.name = "value")
toplot <- merge(toplot, meta, by = "sample")

toplot$group <- toplot[, opt$group_col]
toplot <- ddply(toplot, .(feature, group), mutate, av = mean(value), sd = sd(value))


#toplot <- ddply(toplot, .(feature, group), summarize, av = mean(TPM), se = sd(TPM)/sqrt(length(TPM)), .drop = F)



toplot$group <- NULL
if (!is.null(ord_idx)) 
  toplot[, group.by] <- factor(toplot[, group.by], levels = ord_idx)













meta_cols <- c("sample",group_col)
if(!is.null(opt$timepoint)){
  meta_cols <- c(meta_cols,opt$timepoint)
}
if(!is.null(opt$condition_col)){
  meta_cols <- c(meta_cols,opt$condition_col)
}

toplot <- merge(toplot, meta[,meta_cols])

if(!is.null(opt$timepoint)){
  colnames(toplot)[which(colnames(toplot)==opt$timepoint)] <- "x_axis"
} else {
  colnames(toplot)[which(colnames(toplot)==opt$group_col)] <- "x_axis"
}
if(!is.null(opt$condition_col)){
  colnames(toplot)[which(colnames(toplot)==opt$condition_col)] <- "condition"
}

# set group order ----
if(!is.null(opt$group_order)){
  group_order <- order(meta[,opt$group_order])
  toplot$x_axis <- factor(toplot$x_axis, levels = unique(meta$x_axis[group_order]))
} else {
  toplot[,group_col] <- factor(toplot[,group_col])
}

# palette ----
if(!is.null(opt$pal)){
  pal <- readRDS(opt$pal)
} else {
  library(RColorBrewer)
  pal <- brewer.pal(n = length(levels(toplot$x_axis)), name = "Set1")
  names(pal) <- levels(toplot$x_axis)
}



# 2. Line plot ----
p_line <- ggplot(toplot, aes(x=x_axis, y=av, color=feature, group=feature)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymax = av + se, ymin = av - se, group=feature_col), width = 0.3, lwd = 0.25) +
  theme_pubr(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("average TPM") +
  scale_color_discrete(name="") 

if(!is.null(opt$condition_col)){
  p_line <- p_line + facet_wrap(.~condition)
}
# save pdf
pdf(paste0(outdir,"/line_plot.pdf"), paper="a4", width = 5, height = 4.5) 
print(p)
dev.off()


# 3. bar plot with transcripts faceting ----
# filter on expressed features
toplot_bar <- ddply(toplot, .(feature), mutate, av_feature = mean(value))
# plot
p_bar <- ggplot(unique(toplot[which(toplot$av_feature!=0), c("group", "av", "sd", "feature")]), aes(x = group, y = av, fill = group)) + 
  geom_col(lwd = 0.2, width = 0.75, show.legend = T, color="black") + 
  geom_errorbar(aes(ymax = av + sd, ymin = av - sd), size = 0.2, width = 0.3, linetype = "dashed", lwd = 0.25, position = position_dodge(width = 0.8)) + 
  facet_wrap(~feature, scales = "free_y", ncol = 3) + 
  theme_classic(base_size = 8, ) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  theme(legend.key.size = unit(4, "mm"), axis.title.x = element_blank()) + 
  scale_fill_manual(values = pal) +
  geom_jitter(data = toplot[which(toplot$av_feature!=0),], aes(x = group, y = value), show.legend = F, size = 0.8, width = 0.1, height = 0) +
  ylab("TPM")
# save
outdir <- "Salmon_plot"
pdf(paste0(outdir,"/bar_plot.pdf"), paper="a4", width = 5, height = 4.5) 
print(p_bar)
dev.off()

