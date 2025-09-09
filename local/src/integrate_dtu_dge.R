library(rnaseqDTU)
library(tximport)
library(GenomicFeatures)

# SETUP ----
samps <- read.delim("metadata.txt", header = T)
samps <- samps[grepl("DNMT3B_KO|H9", samps$sample),]
samps <- samps[grepl("D0", samps$sample),]
colnames(samps)[1] <- "sample_id"

files <- file.path("Salmon", samps$sample, "quant.sf")
names(files) <- samps$sample

files <- files[grepl("DNMT3B_KO|H9", files)]
files <- files[grepl("D0", files)]

# count-scale data: column sums are equal to the number of mapped paired-end reads per experiment.
txi <- tximport(files, type="salmon", txOut = TRUE, countsFromAbundance="scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]

# transcript-to-gene mappings
gtf <- "/home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/37/primary_assembly.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf)
txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
colnames(txdf) <- c("gene_id","feature_id")

counts <- merge(txdf, cts, by.x="feature_id", by.y="row.names", all.y=T)
counts<- na.omit(counts)



# DRIMseq ----
library(DRIMSeq)
d <- dmDSdata(counts=counts, samples=samps)
# n sample
n <- nrow(d@samples)
# n of smallest group
n.small <- 2
# Filter
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
# Only genes with more than one isoforms are kept
table(table(counts(d)$gene_id))

design_full <- model.matrix(~0+genotype_general, data=DRIMSeq::samples(d))
my.contrasts <- makeContrasts(DNMT3B_vs_Ctrl=genotype_generalDNMT3BKO - genotype_generalH9, 
                              levels=design_full)

set.seed(1)
BPPARAM <- MulticoreParam(workers = 50)
d <- dmPrecision(d, design=design_full, BPPARAM=BPPARAM)
d <- dmFit(d, design=design_full,BPPARAM=BPPARAM)
d <- dmTest(d, contrast=my.contrasts, BPPARAM=BPPARAM)
res.txp <- DRIMSeq::results(d, level="feature")
saveRDS(d, "drimseq.rds")
write.table(res.txp, gzfile("drimseq_res.3BKO_D0.tab.gz"), quote = F, col.names = T, sep = "\t", row.names = F)

res <- res[which(res$adj_pvalue < 0.05),]

pScreen <- res.txp$pvalue
strp <- function(x) substr(x,1,15)
names(pScreen) <- strp(res.txp$gene_id)
pConfirmation <- matrix(res.txp$pvalue, ncol=1)
rownames(pConfirmation) <- strp(res.txp$feature_id)
tx2gene <- res.txp[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])




# STAGER ----
# which genes contain any evidence of DTU?
# which transcripts in the genes that contain some evidence may be participating in the DTU?
# no more than 5% of the genes that pass screening will either (1) not contain any DTU, so be falsely screened genes, or (2) contain a falsely confirmed transcript.
library(stageR)
stageRObj <- stageRTx(pScreen=pScreen, 
                      pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, 
                      tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05, allowNA=T) 
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=TRUE)
})

sig_genes <- getSignificantGenes(object = stageRObj)
sig_transcripts <- getSignificantTx(object = stageRObj)
write.table(sig_genes, gzfile("drimseq_res.3BKO_D0.sig_genes.gz"), quote = F, col.names = T, sep = "\t", row.names = F)
write.table(sig_transcripts, gzfile("drimseq_res.3BKO_D0.sig_transcripts.gz"), quote = F, col.names = T, sep = "\t", row.names = F)


# EDGER ----
# Integrate diff. gene expression
txi.g <- tximport(files, type="salmon", tx2gene=txdf[,2:1])
library(edgeR)
cts.g <- txi.g$counts
normMat <- txi.g$length
normMat <- normMat / exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts.g/normMat)) + log(colSums(cts.g/normMat))
y <- DGEList(cts.g)
y <- scaleOffset(y, t(t(log(normMat)) + o))
keep <- filterByExpr(y)
y <- y[keep,]

y <- estimateDisp(y, design_full)
fit <- glmFit(y, design_full)
lrt <- glmLRT(fit,contrast = my.contrasts)
tt <- topTags(lrt, n=nrow(y), sort="none")[[1]]

common <- intersect(res$gene_id, rownames(tt))
tt <- tt[common,]
res.sub <- res[match(common, res$gene_id),]

ggplot()

bigpar()
plot(-log10(tt$FDR), -log10(res.sub$adj_pvalue), col=col,
     xlab="Gene expression",
     ylab="Transcript usage")
legend("topright",
       c("DGE","DTE","DTU","null"),
       col=c(1:3,8), pch=20, bty="n")


