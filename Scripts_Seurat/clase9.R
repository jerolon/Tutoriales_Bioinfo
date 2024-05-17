#Funcion recursiva para quitar una libreria cargada junto con todos
#sus dependientes
detachWithParents <- function(ns){
ns <- asNamespace(ns, base.OK = FALSE)
nsname <- getNamespaceName(ns)
users <- getNamespaceUsers(ns)
setwd("/mnt/atgc-d2/bioa03124/jmiranda")
if (length(users)){
  print(gettextf("Erasing %s with its dependencies %s", 
           sQuote(getNamespaceName(ns)), paste(sQuote(users), 
                                               collapse = ", ")))
  sapply(users, detachWithParents)                                           
}
pos <- match(paste0("package:", nsname), search())
if (!is.na(pos)){
  detach(pos = pos)
}
print(paste0("Borrando: ", nsname))
unloadNamespace(nsname)
}
###############Loading libs#######################
#Cargar versiones especificas de librerias       #
#Seurat solo funciona con ggplot2 v3.4 no con 3.5#
##################Loading libs####################
while(packageVersion("ggplot2") >= 3.5){
detachWithParents("ggplot2")
try(library(ggplot2, lib.loc = "/cm/shared/apps/r/4.3.1-studio/lib64/R/library"), TRUE)
}
packageVersion("ggplot2")
#BiocManager::install("plger/scDblFinder")
library(scDblFinder)
#dos lineas para saber si esta funcionando scDbl
sce <- mockDoubletSCE()
sce <- scDblFinder(sce)
#Cargamos la matriz que salió de alevinfry en formato SingleCellExp de Bioconductor
library(fishpond)
#Solo requiere el path a donde esta la carpeta alevin y el archivo quant.json
scexp_obj <- loadFry("./haircells_10X")

#libreria para manipular el objeto ya se cargo
#library(SingleCellExperiment)
#la matriz de conteos "raw" se obtiene con assay()
matriz <- assay(scexp_obj, "counts")
library(Seurat)
packageVersion("Seurat")
#make sure it is v5
#otherwise add /mnt/atgc-d2/bioa03124/.rlibraries/4.3/ to R_LIBS in terminal

#hacer un objeto Seurat con la matriz. Filtra los features que salgan en menos de 3 celulas y las celulas con menos de 200 genes
haircells <- CreateSeuratObject(counts = matriz, project = "kozak2020", min.cells = 3, min.features = 200)

#Features es cualquier caracteristica, conteo, propiedad de las celulas
#añadir propiedades al df de metadata de celulas
haircells[["percent.mt"]] <- PercentageFeatureSet(haircells, pattern = "^mt-")
haircells[["percent.ribo"]] <- PercentageFeatureSet(haircells, pattern = "^rp[sl]")
#hemoglobinas, pero hbp es un factor de transcripcion
haircells[["percent.hemo"]] <- PercentageFeatureSet(haircells, pattern = "^hb[^p]")
haircells[["EGFP"]] <- FetchData(haircells, vars = "transEGFP")

#Graficas pre-filtrado
VlnPlot(haircells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(haircells, features = c("percent.ribo", "percent.hemo", "EGFP"), ncol = 3)
plot1 <- FeatureScatter(haircells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(haircells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(haircells, feature1 = "EGFP", feature2 = "percent.hemo")
plot1  + plot2 + plot3


# Filtering based on MAD --------------------------------------------------
library(dplyr)
#funcion para clasificar como outliers
is_outlier <- function(seurobj, metric, nmads){
  M <- seurobj[[metric]][,1]
  #between nos da las observaciones ques estan entre nmads desviaciones de la mediana
  outlier <- !between(M, median(M) - nmads*mad(M), median(M) + nmads*mad(M))
  return(outlier)
}

#transform to log1p for linearity
haircells[["log1p_nCount"]] <- log1p(haircells[["nCount_RNA"]])
haircells[["log1p_nGenes"]] <- log1p(haircells[["nFeature_RNA"]])
FeatureScatter(haircells, feature1 = "log1p_nCount", feature2 = "log1p_nGenes")
#make outlier based on deviations
haircells[["outlier"]] <- is_outlier(haircells, "log1p_nCount", 5) | 
                          is_outlier(haircells, "log1p_nGenes", 5)|
                          is_outlier(haircells, "percent.hemo", 5)
haircells[["mt.outlier"]] <- is_outlier(haircells, "percent.mt", 3) | haircells[["percent.mt"]] > 8
ggplot(haircells@meta.data) + geom_point(aes(log1p_nCount, log1p_nGenes, col = outlier)) + theme_bw()
ggplot(haircells@meta.data) + geom_point(aes(log1p_nCount, log1p_nGenes, col = mt.outlier)) + theme_bw()

#filtering with the outlier columns
haircells <- subset(haircells, subset = outlier == FALSE & mt.outlier == FALSE)

# Substracting the soup ---------------------------------------------------
library(SoupX)
#Get quick clusters for the soup
haircellpp <- SCTransform(haircells, vars.to.regress = "percent.mt")
haircellpp <- RunPCA(haircellpp)
haircellpp <- FindNeighbors(haircellpp)
haircellpp <- FindClusters(haircellpp)
haircells[["clusters"]] <- haircellpp[["seurat_clusters"]]
#We also copy this even if we dont understand it yet for vizualization
haircells[["PC1"]] <- FetchData(haircellpp, layer = "pca", vars = c("PC_1"))
haircells[["PC2"]] <- FetchData(haircellpp, layer = "pca", vars = c("PC_2"))
rm(haircellpp)

#get the matrix with the original unfiltered number of GEMs, but the same number of genes
soupCells <- CreateSeuratObject(counts = matriz, project = "kozak2020", min.cells = 3)
soupCells <- subset(soupCells, features = Features(haircells))
matriz <- soupCells[["RNA"]]$counts

sc <- SoupX::SoupChannel(matriz, haircells[["RNA"]]$counts, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
sc <- setClusters(sc, setNames(haircells[["clusters"]]$clusters, row.names(haircells[["clusters"]])))
sc <- autoEstCont(sc)
out <- adjustCounts(sc, roundToInt = TRUE)
#Overwrite the count matrix with the corrected counts
haircells[["RNA"]]$counts <- out


# Simulate and identify droplets ------------------------------------------
#use the out object since we have not modified it 
sce <- scDblFinder(SingleCellExperiment(list(counts=out)))
haircells[["scDblFinder_score"]] <- sce$scDblFinder.score
haircells[["scDblFinder_class"]] <- sce$scDblFinder.class
table(haircells[["scDblFinder_class"]])
SaveSeuratRds(haircells, "haircells10X_clase9.RDS")

###Team SCRB####################################################
#Team SCRB                                                     #            
#Team SCRB######################################################
scexp_obj <- loadFry("./mcSCRB")
#la matriz de conteos "raw" se obtiene con assay()
matriz <- assay(scexp_obj, "counts")
#We dont really want features that are zero for all cells
hcSCRB <- CreateSeuratObject(counts = matriz, project = "kozak2020", min.cells = 1, min.features = 200)

#Some  metadata is the same
hcSCRB[["percent.mt"]] <- PercentageFeatureSet(hcSCRB, pattern = "^mt-")
hcSCRB[["percent.ribo"]] <- PercentageFeatureSet(hcSCRB, pattern = "^rp[sl]")
#hemoglobinas, pero hbp es un factor de transcripcion
hcSCRB[["percent.hemo"]] <- PercentageFeatureSet(hcSCRB, pattern = "^hb[^p]")
hcSCRB[["EGFP"]] <- FetchData(hcSCRB, vars = "transEGFP")

#Microwells can be enriched with extra cell features
library(readr)
barcodes <- read_csv(file = "kozak2020_cellIdent.csv")
barcodeOrig <- tibble(barcode = colnames(x = hcSCRB))
#left join on the barcode columns
#just a trick to order the barcodes as in the seurat object
barcodeOrig <- dplyr::left_join(barcodeOrig, barcodes)
hcSCRB <- AddMetaData(object = hcSCRB, metadata = barcodeOrig)
#Do all barcodes correspond?
all(colnames(hcSCRB) == hcSCRB[["barcode"]])

#Graficas pre-filtrado
VlnPlot(hcSCRB, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hcSCRB, features = c("percent.ribo", "percent.hemo", "EGFP"), ncol = 3)
plot1 <- FeatureScatter(hcSCRB, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hcSCRB, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(hcSCRB, feature1 = "EGFP", feature2 = "percent.hemo")
plot1  + plot2 + plot3

hcSCRB[["log1p_nCount"]] <- log1p(hcSCRB[["nCount_RNA"]])
hcSCRB[["log1p_nGenes"]] <- log1p(hcSCRB[["nFeature_RNA"]])
FeatureScatter(hcSCRB, feature1 = "log1p_nCount", feature2 = "log1p_nGenes")
#¿Como se comparan los features y counts entre técnicas?
#Y los conteos de GFP?

#make outlier based on deviations
hcSCRB[["outlier"]] <- is_outlier(hcSCRB, "log1p_nCount", 5) | 
  is_outlier(hcSCRB, "log1p_nGenes", 5)|
  is_outlier(hcSCRB, "percent.hemo", 5)
hcSCRB[["mt.outlier"]] <- is_outlier(hcSCRB, "percent.mt", 3) | hcSCRB[["percent.mt"]] > 8
ggplot(hcSCRB@meta.data) + geom_point(aes(log1p_nCount, log1p_nGenes, col = outlier)) + theme_bw()
ggplot(hcSCRB@meta.data) + geom_point(aes(log1p_nCount, log1p_nGenes, col = mt.outlier)) + theme_bw()

#Do we filter?
hcSCRB <- subset(hcSCRB, subset = outlier == FALSE & mt.outlier == FALSE)

# Soup SCRB---------------------------------------------------
#Get quick clusters for the soup
haircellpp <- SCTransform(hcSCRB, vars.to.regress = "percent.mt")
haircellpp <- RunPCA(haircellpp)
haircellpp <- FindNeighbors(haircellpp)
haircellpp <- FindClusters(haircellpp)
hcSCRB[["clusters"]] <- haircellpp[["seurat_clusters"]]
#We also copy this even if we dont understand it yet for vizualization
hcSCRB[["PC1"]] <- FetchData(haircellpp, layer = "pca", vars = c("PC_1"))
hcSCRB[["PC2"]] <- FetchData(haircellpp, layer = "pca", vars = c("PC_2"))
rm(haircellpp)

#The cells are physically separated during library construction
#Therefore, it does not make sense to calculate the soup
matriz <- hcSCRB[["RNA"]]$counts

sc <- SoupX::SoupChannel(matriz, matriz, calcSoupProfile = FALSE)
#This will fail, since all cells have enough UMIs to not be considered soup
sc <- estimateSoup(sc)
sc$metaData %>% arrange(nUMIs) %>% head

#SCRB droplets?  ------------------------------------------
#use the out object since we have not modified it 
sce <- scDblFinder(SingleCellExperiment(list(counts=matriz)))
hcSCRB[["scDblFinder_score"]] <- sce$scDblFinder.score
hcSCRB[["scDblFinder_class"]] <- sce$scDblFinder.class
table(hcSCRB[["scDblFinder_class"]])
#¿More evidence to delete the outliers?
ggplot() + geom_point(data = hcSCRB[[]], aes(log1p_nCount,log1p_nGenes), col = "lightgray") + geom_point(data = filter(hcSCRB[[]], scDblFinder_class == "doublet"),  aes(log1p_nCount,log1p_nGenes), col = "red") + theme_bw()
SaveSeuratRds(hcSCRB, "haircellsSCRB_clase9.RDS")


# Back to 10X -------------------------------------------------------------
#Solo requiere el path a donde esta la carpeta alevin y el archivo quant.json
scexp_obj <- loadFry("./supportcells_10X/")

#la matriz de conteos "raw" se obtiene con assay()
matriz <- assay(scexp_obj, "counts")

#hacer un objeto Seurat con la matriz. Filtra los features que salgan en menos de 3 celulas y las celulas con menos de 200 genes
supportcells <- CreateSeuratObject(counts = matriz, project = "kozak2020", min.cells = 3, min.features = 200)

#Features es cualquier caracteristica, conteo, propiedad de las celulas
#añadir propiedades al df de metadata de celulas
supportcells[["percent.mt"]] <- PercentageFeatureSet(supportcells, pattern = "^mt-")
supportcells[["percent.ribo"]] <- PercentageFeatureSet(supportcells, pattern = "^rp[sl]")
#hemoglobinas, pero hbp es un factor de transcripcion
supportcells[["percent.hemo"]] <- PercentageFeatureSet(supportcells, pattern = "^hb[^p]")
supportcells[["EGFP"]] <- FetchData(supportcells, vars = "transEGFP")

#Graficas pre-filtrado
VlnPlot(supportcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(supportcells, features = c("percent.ribo", "percent.hemo", "EGFP"), ncol = 3)
plot1 <- FeatureScatter(supportcells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(supportcells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(supportcells, feature1 = "EGFP", feature2 = "percent.hemo")
plot1  + plot2 + plot3
FeatureScatter(supportcells, feature1 = "percent.ribo", feature2 = "percent.hemo")

# Filtering based on MAD --------------------------------------------------
#transform to log1p for linearity
supportcells[["log1p_nCount"]] <- log1p(supportcells[["nCount_RNA"]])
supportcells[["log1p_nGenes"]] <- log1p(supportcells[["nFeature_RNA"]])
FeatureScatter(supportcells, feature1 = "log1p_nCount", feature2 = "log1p_nGenes")
#make outlier based on deviations
supportcells[["outlier"]] <- is_outlier(supportcells, "log1p_nCount", 5) | 
  is_outlier(supportcells, "log1p_nGenes", 5)|
  is_outlier(supportcells, "percent.hemo", 5)
supportcells[["mt.outlier"]] <- is_outlier(supportcells, "percent.mt", 3) | supportcells[["percent.mt"]] > 8
ggplot(supportcells@meta.data) + geom_point(aes(log1p_nCount, log1p_nGenes, col = outlier)) + theme_bw()
ggplot(supportcells@meta.data) + geom_point(aes(log1p_nCount, log1p_nGenes, col = mt.outlier)) + theme_bw()

#filtering with the outlier columns
supportcells <- subset(supportcells, subset = outlier == FALSE & mt.outlier == FALSE)

# Substracting the soup ---------------------------------------------------
#Get quick clusters for the soup
soupport <- SCTransform(supportcells, vars.to.regress = "percent.mt")
soupport <- RunPCA(soupport)
soupport <- FindNeighbors(soupport)
soupport <- FindClusters(soupport)
supportcells[["clusters"]] <- soupport[["seurat_clusters"]]
#We also copy this even if we dont understand it yet for vizualization
supportcells[["PC1"]] <- FetchData(soupport, layer = "pca", vars = c("PC_1"))
supportcells[["PC2"]] <- FetchData(soupport, layer = "pca", vars = c("PC_2"))
rm(soupport)

#get the matrix with the original unfiltered number of GEMs, but the same number of genes
soupCells <- CreateSeuratObject(counts = matriz, project = "kozak2020", min.cells = 3)
soupCells <- subset(soupCells, features = Features(supportcells))
matriz <- soupCells[["RNA"]]$counts

sc <- SoupX::SoupChannel(matriz, supportcells[["RNA"]]$counts, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
sc <- setClusters(sc, setNames(supportcells[["clusters"]]$clusters, row.names(supportcells[["clusters"]])))
sc <- autoEstCont(sc)
out <- adjustCounts(sc, roundToInt = TRUE)
#Overwrite the count matrix with the corrected counts
supportcells[["RNA"]]$counts <- out


# Simulate and identify droplets ------------------------------------------
#use the out object since we have not modified it 
sce <- scDblFinder(SingleCellExperiment(list(counts=out)))
supportcells[["scDblFinder_score"]] <- sce$scDblFinder.score
supportcells[["scDblFinder_class"]] <- sce$scDblFinder.class
table(supportcells[["scDblFinder_class"]])
SaveSeuratRds(supportcells, "supportcells10X_clase9.RDS")
ggplot() + geom_point(data = supportcells[[]], aes(PC1, PC2), col = "lightgray") + geom_point(data = filter(supportcells[[]], scDblFinder_class == "doublet"), aes(PC1, PC2), col = "red") + theme_bw()
ggplot() + geom_point(data = supportcells[[]], aes(log1p_nCount, log1p_nGenes), col = "lightgray") + geom_point(data = filter(supportcells[[]], scDblFinder_class == "doublet"), aes(log1p_nCount, log1p_nGenes), col = "red") + theme_bw()


haircells[["RNA"]]@meta.data

haircells$nCount_data <- apply(haircells[["RNA"]]$data, 2, sum)
> haircells$nCount_scale <- apply(haircells[["RNA"]]$scale.data, 2, sum)
> VlnPlot(haircells, features = c("nCount_RNA", "nCount_data", "nCount_scale"))
> VlnPlot(haircells, features = c("nCount_RNA", "nCount_data", "nCount_scale"), group.by = "orig.ident")