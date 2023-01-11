#! /usr/bin/env Rscript

'RNA-sequencing data visualisation
 
 plot-isoform.r
 
 Create barplot and lineplot from counts data 

Usage:
   plot-isoform.r [-o <group_order>, -c <condition_col>] <input> <metadata> <group_col> <output>

Arguments:
    input  feature count matrix with header
    metadata  metadata file path
    group_col  name of metadata column to use for averaging
    output  output figure prefix
              
Options:
  -c, --condition_col  name of metadata column to use for diverse groups [Default: unique condition for all groups]
  -o, --group_order_col  name of metadata column to use for ordering groups [Default: alphanumeric order]
  -p, --palette  custom made palette [Default: RColorBrewer]

Author:
   Olilab' -> doc

opts <- docopt(doc)

suppressWarnings(suppressMessages(library(RNAseqRtools)))
suppressWarnings(suppressMessages(library(docopt)))

# File Read ----
# https://stackoverflow.com/questions/26152998/how-to-make-r-script-takes-input-from-pipe-and-user-given-parameter
# if the input is stdin one can do 
# cat fragment.txt | ./plot_atac_frag_distribution.R --poly --pdf stdin  out.pdf
# cat fragment.txt | ./plot_atac_frag_distribution.R --poly --pdf - out.pdf
# ./plot_atac_frag_distribution.R --poly --pdf <(cat fragment.txt)  out.pdf
OpenRead <- function(arg) {
  if (arg %in% c("-", "/dev/stdin")) {
    file("stdin", open = "r")
  } else if (grepl("^/dev/fd/", arg)) {
    fifo(arg, open = "r")
  } else {
    file(arg, open = "r")
  }
}

#z <- read.table(OpenRead(arguments$input), header = FALSE)
z <- read.delim("multi_SalmonQuant.TPM_matrix.DNMT3B.gz", header = T)
#meta <- read.delim(arguments$metadata, header=T)
meta <- read.delim("metadata.txt", header = T)

feature_col <- names(z)[1]
#group <- arguments$group_col
group_col <- "condition"


toplot <- melt(z, id.vars = feature_col, variable.name = "sample", value.name = "count")
toplot <- merge(toplot, meta)
toplot <- ddply(toplot, c(feature_col, group_col), summarize, 
                av = mean(count), se = sd(count)/sqrt(length(count)), 
                .drop = F)

if(!is.null(opt$condition_col)){
  condition_col <- opt$condition_col
  toplot <- merge(toplot, meta[, c(group_col,condition_col)])
  colnames(toplot)[5] <- "cond"
} else {
  toplot$cond <- "none"
}


if(!is.null(opt$group_order)){
  group_order <- order(meta[,opt$group_order])
  #group_order <- order(meta[,"day_number"])
  toplot[,group_col] <- factor(toplot[,group_col], levels = unique(meta[,group_col][group_order]))
} else {
  toplot[,group_col] <- factor(toplot[,group_col])
}





##################
#
#  Line plot ----
#

p <- ggplot(toplot, aes_string(x=condition, y="av", color=feature_col, group=feature_col)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymax = av + se, ymin = av - se, group=feature_col), width = 0.3, lwd = 0.25) +
  theme_pubr(base_size = 8)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") + ylab("average TPM") +
  scale_color_discrete(name="") +
  scale_y_continuous(trans='log2')

pdf("multi_SalmonQuant.TPM_matrix.DNMT3B_line.pdf", paper="a4", width = 4, height = 4) 
print(p)
dev.off()

n_group <- length(levels(toplot[,group_col]))
n_cond <- length(unique(toplot$cond))

# Set color palette
if(!is.null(opt$pal)){
  pal <- opt$pal
} else {
  available_pals <- c("Greys","Blues","Greens","Oranges","Purples","Reds")
  pal_all <-list()
  pals <- list()
  for(i in 1:n_cond){
    group_cond <- unique(toplot[toplot$cond==unique(toplot$cond)[i],group_col])
    pals <- c(RColorBrewer::brewer.pal(length(group_cond), "Blues")
  }
  pal[[group_col]] <- pals
    
  for (col in available_pals[1:length(cond)]) {
    pal[[group]][[col]] <- RColorBrewer::brewer.pal()
  }
}

pal <- RColorBrewer::brewer.pal(n_group, "Blues")
names(pal) <- levels(rds$samples[[group]])

pal_default <- list()
pal_1 <- RColorBrewer::brewer.pal(9, "Blues")[c(3,5,7,9)]
names(pal_1) <- paste0("shDUBR_",c("NPC","TD_Day10","TD_Day20","TD_Day30"))
pal_2 <- RColorBrewer::brewer.pal(9, "Greens")[c(3,5,7,9)]
names(pal_2) <- paste0("shLINC01578_",c("NPC","TD_Day10","TD_Day20","TD_Day30"))
pal_3 <- RColorBrewer::brewer.pal(9, "Oranges")[c(3,5,7,9)]
names(pal_3) <- paste0("shMIR99AHG_",c("NPC","TD_Day10","TD_Day20","TD_Day30"))
pal_4 <- RColorBrewer::brewer.pal(9, "Greys")[c(3,5,7,9)]
names(pal_4) <- paste0("shPLKO_",c("NPC","TD_Day10","TD_Day20","TD_Day30"))

pal_default[["cond_time"]] <- c(pal_1,pal_2,pal_3,pal_4)

                    

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




# Stacked barplot









