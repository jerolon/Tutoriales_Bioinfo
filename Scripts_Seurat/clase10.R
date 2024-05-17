source("functions.R")
###############Loading libs#######################
#Cargar versiones especificas de librerias       #
#Seurat solo funciona con ggplot2 v3.4 no con 3.5#
##################Loading libs####################
while(packageVersion("ggplot2") >= 3.5){
  detachWithParents("ggplot2")
  try(library(ggplot2, lib.loc = "/cm/shared/apps/r/4.3.1-studio/lib64/R/library"), TRUE)
}
packageVersion("ggplot2")
library(Seurat)
packageVersion("Seurat")

haircells <- LoadSeuratRds("./9_qcfilters/haircells10X_clase9.RDS")
haircells <- NormalizeData(haircells)
haircells <- FindVariableFeatures(haircells, selction.method = "vst", nFeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(haircells), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(haircells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(haircells)
haircells <- ScaleData(haircells, features = all.genes)
var <- VariableFeatures(haircells)
haircells <- RunPCA(haircells, features = var)
VizDimLoadings(haircells, dims = 1:2, reduction = "pca")
ElbowPlot(haircells)
haircells <- FindNeighbors(haircells, dims = 1:50, k.param = 30)
haircells <- FindClusters(haircells, resolution = 0.3)
haircells <- RunUMAP(haircells, dims = 1:19)
haircells <- RunTSNE(haircells, reduction = "pca", dims = 1:19, perplexity = 5, step = 1000)
DimPlot(haircells, reduction = "umap")
var_genes_haircells <- subset(haircells, features = var)
var_genes_haircells[["RNA"]]$counts <- NULL
var_genes_haircells[["RNA"]]$data <- NULL
neighbour_ambient <- FindNeighbors(t(var_genes_haircells[["RNA"]]$scale.data), l2.norm = TRUE, return.neighbor = T, k.param = 30)
neighbour_pca <- FindNeighbors(haircells, dims = 1:50, reduction = "pca", l2.norm = TRUE, k.param = 30,return.neighbor = TRUE)
neighbour_tsne <- FindNeighbors(haircells, dims = 1:2, reduction = "tsne", l2.norm = TRUE, k.param = 30, return.neighbor = TRUE)
neighbour_umap <- FindNeighbors(haircells, dims = 1:2, reduction = "umap", l2.norm = TRUE, k.param = 30,return.neighbor = TRUE)

#Local distortion
neighbourhood_sim <- data.frame(pca_var = cbind(neighbour_ambient@nn.idx, neighbour_pca@neighbors$RNA.nn@nn.idx) %>% apply(1, jaccard, 30),
                                var_tsne = cbind(neighbour_ambient@nn.idx, neighbour_tsne@neighbors$RNA.nn@nn.idx) %>% apply(1, jaccard, 30),
                                var_umap = cbind(neighbour_ambient@nn.idx, neighbour_umap@neighbors$RNA.nn@nn.idx) %>% apply(1, jaccard, 30)
                                )
ggplot(neighbourhood_sim) + geom_density(aes(pca_var)) + geom_density(aes(var_tsne), col = "blue") + geom_density(aes(var_umap), col = "red")


#Knn accuracy
#Cluster identity for the cells and each neighbour
neighbour_ambient_clusterid <- apply(neighbour_ambient@nn.idx, c(1,2), function(x){haircells[[]]$clusters[x]})
neighbour_pca_clusterid <- apply(neighbour_pca@neighbors$RNA.nn@nn.idx, c(1,2), function(x){haircells[[]]$clusters[x]})
neighbour_tsne_clusterid <- apply(neighbour_tsne@neighbors$RNA.nn@nn.idx, c(1,2), function(x){haircells[[]]$clusters[x]})
neighbour_umap_clusterid <- apply(neighbour_umap@neighbors$RNA.nn@nn.idx, c(1,2), function(x){haircells[[]]$clusters[x]})

#function to give the majority vote from a vector
majority <- function(arreglo){
  conteos<- arreglo %>% table %>% as.data.frame
  return(conteos$.[1])
} 

knn_accuracy <- data.frame(embedding = rep("ambient", length(Cells(haircells))), cluster = haircells[[]]$clusters, knnvote = apply(neighbour_ambient_clusterid[,2:30], 1, majority))

knn_accuracy <- rbind(knn_accuracy,
                      data.frame(embedding = rep("pca", length(Cells(haircells))), cluster = haircells[[]]$clusters, knnvote = apply(neighbour_pca_clusterid[,2:30], 1, majority)),
                      data.frame(embedding = rep("tsne", length(Cells(haircells))), cluster = haircells[[]]$clusters, knnvote = apply(neighbour_tsne_clusterid[,2:30], 1, majority)),
                      data.frame(embedding = rep("umap", length(Cells(haircells))), cluster = haircells[[]]$clusters, knnvote = apply(neighbour_umap_clusterid[,2:30], 1, majority))
                      )
knn_accuracy %>% mutate(acertado = cluster == knnvote) %>% group_by(embedding) %>% summarise(knnac = sum(acertado)/n())

supportcells <- LoadSeuratRds("./9_qcfilters/supportcells10X_clase9.RDS")
SCRB <- LoadSeuratRds("./9_qcfilters/haircells1SCRB_clase9.RDS")
