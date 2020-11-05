install.packages(
  "e1071",
  "ggplot2",
  "gplots"
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")
BiocManager::install("sva")
BiocManager::install("SummarizedExperiment")
BiocManager::install("pcaMethods")
BiocManager::install("affy")
BiocManager::install("GEOquery")
