library(Seurat)
library(ggplot2)

source(file = "../../gastruloids_sc_Lescroart/analysis/gastruloid_timeserie_scRNA-seq/scripts/utilities/00_generalDeps.R")

#Load data
cardioWT <- readRDS("../50_rdsObjects/02_normANDscaled_cardioWT_complete.rds")
umap <- read.table(
    file = "../04_process/DR_umap_cardioWT_complete.csv",
    sep = ",",
    header = TRUE)
cardioWT[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_")


cardioKO <- readRDS("../50_rdsObjects/02_normANDscaled_cardioKO_complete.rds")
umap <- read.table(
    file = "../04_process/DR_umap_cardioKO_complete.csv",
    sep = ",",
    header = TRUE)
cardioKO[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_")


fibroWT <- readRDS("../50_rdsObjects/02_normANDscaled_fibroWT_complete.rds")
umap <- read.table(
    file = "../04_process/DR_umap_fibroWT_complete.csv",
    sep = ",",
    header = TRUE)
fibroWT[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_")


fibroKO <- readRDS("../50_rdsObjects/02_normANDscaled_fibroKO_complete.rds")
umap <- read.table(
    file = "../04_process/DR_umap_fibroKO_complete.csv",
    sep = ",",
    header = TRUE)
fibroKO[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_")

datasetsList <- as.list(c(cardioWT, cardioKO, fibroWT, fibroKO))
datasetsNames <- c("cardioWT", "cardioKO", "fibroWT", "fibroKO")
names(datasetsList) <- datasetsNames
# datasetsList <- as.list(c(cardioWT))
# datasetsNames <- c("cardioWT")
# names(datasetsList) <- datasetsNames

goiList <- read.table("../01_requirements/genesOfInterest_2023-04-13.tsv")
goiList <- sort(goiList$V1)


invisible(sapply(datasetsNames, function(DATASET) {
    print(DATASET)
    genesToRemove <- goiList[!(goiList %in% rownames(datasetsList[[DATASET]]))]
    print(genesToRemove)
    
    goiList <- goiList[goiList %in% rownames(datasetsList[[DATASET]])]
    
    pdf(file = paste0("../52_figures/", DATASET, "/50_featurePlots_", DATASET, "_", Sys.Date(), ".pdf"))

    for (goi in goiList){
        print(FeaturePlot(datasetsList[[DATASET]], goi, reduction = "umap", order = TRUE) +
                  ggtitle(paste0(goi, " expression in ", DATASET)))
    }

    dev.off()
}))

