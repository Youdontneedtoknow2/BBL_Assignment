# Install packages 
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "DESeq2", "edgeR", "clusterProfiler", "org.Hs.eg.db", "pheatmap"))

# R packages
install.packages(c("tidyverse", "ggthemes", "ggsci", "ggplot2", "umap"))

# Load packages
library(GEOquery)
library(DESeq2)
library(tidyverse)
library(ggthemes)
library(ggsci)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(umap)

# Download RNA-SEQ dataset from GEO
gse_data <- getGEO("GSE21942", GSEMatrix = TRUE)

# 1. Count Data/ expression data
count_data <- exprs(gse_data[[1]])  # Extract expression matrix

# 2. Meta data
metadata <- pData(gse_data[[1]])         # Extract metadata

# View the dataset structure
head(count_data)
head(metadata)

#count data processing
# data normalization
# Log transformation
normalized_count_data <- log2(count_data + 1)

# Filter low-expressed genes
low_expressed_genes <- normalized_data[rowMeans(normalized_count_data > 10) > 0.5, ]

# GENE SYMBOL
## 1.1.Finding gene symbol of count data
platform <- annotation(gse_data[[1]]) # checking platform of dataset
platform

BiocManager::install("hgu133plus2.db", version = "3.20")
library(hgu133plus2.db)

# 1.2.Map probe IDs to gene symbols
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = rownames(count_data),
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

library(dplyr)

# 1.3.Add gene symbol column to count_data
N.count_data <- as.data.frame(normalized_count_data) |> 
  mutate(GeneSymbol = gene_symbols)

# 1.4.Reshape the data into long format
data_long <- N.count_data |> 
  pivot_longer(cols = -GeneSymbol, names_to = "samples", values_to = "value")

# View the long-format data
head(data_long)


# 2.1.Subsetting  meta data

meta.data <- metadata |> 
  dplyr::select(1, 2, 33)


# join counts and meta data
final_data <- data_long |> 
  left_join(meta.data, by = c('samples' = 'geo_accession')) |> 
  rename(Sample_type = title)

# export final data
write.csv(final_data, "data/GSE21942_normalized_processed.csv, row.names = F")

# export count data
write.csv(count_data, "data/count_data.csv, row.names = F")




























































