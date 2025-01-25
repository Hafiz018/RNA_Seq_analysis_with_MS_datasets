#install R packages
install.packages("tidyverse")
install.packages("openxlsx")
install.packages("ggpubr")
install.packages("umap")
install.packages("BiocManager")
#install bioconductor packages
#install R packages
install.packages("tidyverse")
install.packages("openxlsx")
install.packages("ggpubr")
install.packages("umap")
#install bioconductor packages
BiocManager::install(c("GEOquery", "TCGAbiolinks", "DESeq2"))
BiocManager::install("Biobase")

#load packages
library(tidyverse)
library(GEOquery)
library(ggpubr)
library(openxlsx)
library(TCGAbiolinks)
library(DESeq2)
library(Biobase)
library(limma)

#import dataset
gset <- getGEO("GSE21942", GSEMatrix = TRUE, AnnotGPL=TRUE)

#extract expression matrix
expression_matrix_data<- exprs(gset[[1]])

#extract meta data
metadata<- pData(gset[[1]])

# group membership for all samples
gsms <- "00000000000000011111111111111"
sml <- strsplit(gsms, split="")[[1]]


# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }


#Differential expression analysis 
#import gene table
gene_table <- read_table("data/GSE21942.top.table.tsv")

#convert data structure
gene_table$ID <- as.factor(gene_table$ID)
gene_table$Genesymbol <- as.factor(gene_table$Genesymbol)
gene_table$Genetitle <- as.factor(gene_table$Genetitle)
head(gene_table)

#identify 10 over expressed significant genes
over_expressed_genes <- gene_table |> 
  filter(adjPVal <= 0.05 & logFC > 1 ) |> 
  arrange(desc(logFC)) |> 
  head(10)
write.xlsx(over_expressed_genes, "outputs/over_expressed_top10_genes.xlsx")
# identify 10 under expressed significant genes
under_expressed_genes <- gene_table |> 
  filter(adjPVal <= 0.05 & logFC < -1  ) |> 
  arrange(desc(logFC)) |> 
  head(10)
write.xlsx(under_expressed_genes, "outputs/under_expressed_top10_genes.xlsx")
#join significant gene
significant_gene <- bind_rows(over_expressed_genes, under_expressed_genes)

#export into csv
write.csv(significant_gene, "outputs/top20_significant_genes.csv")


