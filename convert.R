library(Seurat)
library(SeuratObject)
library(SeuratDisk)

myargs <- commandArgs(trailingOnly=TRUE)

inputFile <- myargs[1] # ends by .rds
outputFile <- myargs[2] # ends with .h5Seurat

cat("\n-----------------------------------------------------------------------------\n")
print(paste0("Conversion of the file ", inputFile))
cat("\n-----------------------------------------------------------------------------\n")

#' Load the Seurat Object
cat("Load Seurat Object (SO)\n\n")
SO <- readRDS(inputFile) 

SO$phase <- Idents(SO)

#' Save as a h5Seurat file format
cat("Save SO as a h5Seurat file\n\n")
SeuratDisk::SaveH5Seurat(SO, outputFile, overwrite = TRUE)

#' Convert to a h5ad file format, readable by scanpy
cat("Do the conversion from h5Seurat to h5ad file format\n\n")
SeuratDisk::Convert(outputFile, dest = "h5ad", assay = "integrated", overwrite = TRUE)

file.remove(outputFile)