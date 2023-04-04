#! /usr/bin/env Rscript

########
# 0. libraries ----
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(plyr))

############
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
parser$add_argument("--group_order", action="store", default = NULL,
                    help = "name of metadata column to use for x axis ordering [Default: alphanumeric order]")
parser$add_argument("--pal", action="store", default=NULL,
                    help = "custom made palette as RDS [Default: RColorBrewer]")
parser$add_argument("--outdir", action="store", default=".",
                    help = "output directory [Default: .]")
parser$add_argument("--gene_id_file", action="store", default=NULL,
                    help = "file path with 2 columns: enst;geneID [Default: NULL]")
# create args object
opt <- parser$parse_args()



########
# 2. assign arguments to obj  ----

# create outdir
if(!is.null(opt$outdir)) outdir <- opt$outdir else outdir <- "."
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
# read input files
counts <- read.delim(opt$count_matrix, header = T, row.names = 1)
features <- scan(opt$feature, character(), quote = "")
meta <- read.delim(opt$metadata, header = T)



##########
# 3. create toplot obj ----

# convert all "-" to "."
meta <- as.data.frame(apply(meta, 2, function(x) gsub(x, pattern = "-", replacement = ".")))
# filter TPM matrix with enst of interest
counts <- counts[row.names(counts) %in% features,]
# average across group
toplot <- reshape2::melt(as.matrix(counts), varnames = c("feature","sample"), value.name = "value")
toplot <- merge(toplot, meta, by = "sample")
colnames(toplot)[which(colnames(toplot)==opt$group_col)] <- "group_col"

if(!is.null(opt$gene_id_file)){
  r <- read.delim(opt$gene_id_file, header = F, sep = ";")
  toplot <- merge(toplot, r, by.x="feature", by.y="V1")
  colnames(toplot)[which(colnames(toplot)=="feature")] <- "feature_enst"
  colnames(toplot)[which(colnames(toplot)=="V2")] <- "feature"
}

toplot <- ddply(toplot, .(feature, group_col), mutate, av = mean(value), sd = sd(value))

#toplot$group_col <- NULL
#if (!is.null(opt$time_order_col)){
#  toplot[, group.by] <- factor(toplot[, group.by], levels = ord_idx)
#}

#meta_cols <- c("sample",group_col)
#if(!is.null(opt$timepoint)){
#  meta_cols <- c(meta_cols,opt$timepoint)
#}
#if(!is.null(opt$condition_col)){
#  meta_cols <- c(meta_cols,opt$condition_col)
#}
#toplot <- merge(toplot, meta[,meta_cols])

# set x axis and faceting columns
if(!is.null(opt$timepoint)){
  colnames(toplot)[which(colnames(toplot)==opt$timepoint)] <- "x_axis"
} else {
  colnames(toplot)[which(colnames(toplot)==opt$group_col)] <- "x_axis"
}
if(!is.null(opt$condition_col)){
  colnames(toplot)[which(colnames(toplot)==opt$condition_col)] <- "condition_faceting"
}

# set group order
if(!is.null(opt$group_order)){
  group_order <- order(meta[,opt$group_order])
  if(!is.null(opt$timepoint)){
    toplot$x_axis <- factor(toplot$x_axis, levels = unique(meta[group_order,opt$timepoint]))
  } else {
    toplot$x_axis <- factor(toplot$x_axis, levels = unique(meta[group_order,opt$group_col]))
  }
} else {
  toplot[,group_col] <- factor(toplot[,group_col])
}



#########
# palette ----

if(!is.null(opt$pal)){
  pal <- readRDS(opt$pal)
} else {
  library(RColorBrewer)
  pal <- list()
  pal_cond <- brewer.pal(n = length(levels(toplot$condition_faceting)), name = "Set1")
  names(pal_cond) <- levels(toplot$condition_faceting)
  pal_feature <- c(brewer.pal(n=8, name = "Dark2"), brewer.pal(n=8, name = "Set1"))
  names(pal_feature) <- unique(toplot$feature)
  pal[["condition_faceting"]] <- pal_cond
  pal[["feature"]] <- pal_feature
}



# 2. Line plot ----
p_line_cond_facet <- ggplot(toplot, aes(x=x_axis, y=av, color=feature, group=feature)) +
  geom_line(linewidth=0.3) +
  geom_point(size=0.7) +
  geom_errorbar(aes(ymax = av + sd, ymin = av - sd, group=feature), width = 0.2, lwd = 0.25, size = 0.2) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("average TPM") +
  facet_wrap(.~condition_faceting, nrow = length(unique(toplot$condition_faceting))) + 
  scale_color_manual(name="", values = pal[["feature"]])

if(!is.null(opt$condition_col)){
  p_line_cond_facet_scaled <- p_line_cond_facet + 
    facet_wrap(.~condition_faceting, scales = "free_y", nrow = length(unique(toplot$condition_faceting)))
  p_line_cond_facet <- p_line_cond_facet +
    facet_wrap(.~condition_faceting, nrow = length(unique(toplot$condition_faceting)))
}

p_line_feature_facet <- ggplot(toplot, aes(x=x_axis, y=av, color=condition_faceting, group=condition_faceting)) +
  geom_line(linewidth=0.3) +
  geom_point(size=0.7) +
  facet_wrap(.~feature, scales = "free_y") +
  geom_errorbar(aes(ymax = av + sd, ymin = av - sd, group=feature), width = 0.3, lwd = 0.25, size = 0.2) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("average TPM") +
  scale_color_manual(values = pal[["condition_faceting"]], name="") 

# save pdf
pdf(paste0(outdir,"/line_plot.pdf"), paper="a4", width = 4, height = 4.5) 
print(p_line_cond_facet)
dev.off()

pdf(paste0(outdir,"/line_plot_feature.pdf"), paper="a4", width = 5, height = 4.5) 
print(p_line_feature_facet)
dev.off()



# 3. bar plot with transcripts faceting ----
# filter on expressed features
barplot=F
if(barplot){
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
}


