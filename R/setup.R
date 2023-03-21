if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SCENIC", "GENIE3", "GRNBoost2", "AUCell", "JASPAR2022", "Seurat", "igraph"))

devtools::install_github("aertslab/SCENIC")

library(SCENIC)
library(GENIE3)
library(GRNBoost2)
library(AUCell)
library(JASPAR2022)
library(Seurat)
library(igraph)
