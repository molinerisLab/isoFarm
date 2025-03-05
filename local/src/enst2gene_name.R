library("biomaRt")
ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")

values<- scan("/home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/32/pippo", what = "character")

getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol"), filters = "ensembl_transcript_id", values = values, mart= ensembl)
