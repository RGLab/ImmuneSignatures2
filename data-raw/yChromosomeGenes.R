library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ychrom <- getBM(attributes = c("chromosome_name", "entrezgene_id", "hgnc_symbol"),
                filters = "chromosome_name",
                values = "Y",
                mart = mart)
yChromGenes <- as.character(ychrom$hgnc_symbol)
yChromGenes <- yChromGenes[ !is.na(yChromGenes) & yChromGenes != "" ]
usethis::use_data(yChromGenes, overwrite = TRUE)
