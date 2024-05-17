setwd("/mnt/atgc-d2/bioa03124/jmiranda/")
library(ggplot2, lib.loc = "/cm/shared/apps/r/4.3.1-studio/lib64/R/library")
library(Seurat)
library(patchwork)
library(dplyr)

brain <- readRDS("visium_brain_clase12.RDS")
allen_reference <- readRDS("clase15/allen_cortex.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance

allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

SpatialDimPlot(brain, label.size = 7, label = TRUE)
cortex <- subset(brain, idents=c(2,3,4,6,7))
#Correr las veces que haga falta
SpatialDimPlot(cortex, cells.highlight = WhichCells(cortex, expression=anterior1_imagerow > 400))
cortex <- subset(cortex, anterior1_imagerow > 400, invert = TRUE)
SpatialDimPlot(cortex, label.size = 7, label = TRUE)
SpatialDimPlot(cortex, cells.highlight = WhichCells(cortex, expression=anterior1_imagecol < 100))
cortex <- subset(cortex, anterior1_imagecol < 100, invert = TRUE)

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)

cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi",
                                        features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- SVFInfo(cortex, method = "moransi", status = TRUE) %>%
  dplyr::filter(!is.na(variable)) %>% 
  arrange(rank) %>% head(n=4) %>% rownames()
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)

SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
                                        "L6b", "Oligo"), pt.size.factor = 1, ncol = 5, crop = FALSE, alpha = c(0.1, 1))

library(SeuratData)
slide.seq <- LoadData("ssHippo")

plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)

SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, 
                                                              idents = c(1,6, 13)), facet.highlight = TRUE)

#Referencia single cell
ref <- readRDS("clase15/mouse_hippocampus_reference.rds")
ref <- UpdateSeuratObject(ref)

anchors <- FindTransferAnchors(reference = ref, query = slide.seq, normalization.method = "SCT",
                               npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype, prediction.assay = TRUE,
                                  weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay

DefaultAssay(slide.seq) <- "predictions"
SpatialFeaturePlot(slide.seq, features = c("Dentate Principal cells", "CA3 Principal cells",
                                           "Entorhinal cortex","Endothelial tip", "Ependymal", 
                                           "Oligodendrocyte"), alpha = c(0.1, 1))

slide.seq$predicted.id <- GetTransferPredictions(slide.seq)
Idents(slide.seq) <- "predicted.id"
SpatialDimPlot(slide.seq, cells.highlight = 
                 CellsByIdentities(object = slide.seq, idents = c("CA3 Principal cells",
                                                                  "Dentate Principal cells", 
                                                                  "Endothelial tip")), facet.highlight = TRUE)


slide.seq$predicted.id <- GetTransferPredictions(slide.seq)
Idents(slide.seq) <- "predicted.id"
SpatialDimPlot(slide.seq, cells.highlight = 
                 CellsByIdentities(object = slide.seq, idents = c("CA3 Principal cells",
                                                                  "Dentate Principal cells", 
                                                                  "Endothelial tip")), facet.highlight = TRUE)

DefaultAssay(slide.seq) <- "SCT"
slide.seq <- FindSpatiallyVariableFeatures(slide.seq, assay = "SCT", slot = "scale.data", features = VariableFeatures(slide.seq)[1:1000],
                                           selection.method = "moransi", x.cuts = 100, y.cuts = 100)

top.clusters <- SVFInfo(slide.seq, method = "moransi", status = TRUE) %>%
  dplyr::filter(!is.na(variable)) %>% 
  arrange(rank) %>% head(n=6) %>% rownames()
SpatialFeaturePlot(slide.seq, features = 
                     top.clusters, ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")

#Spot deconvolution
library(spacexr)
Idents(ref) <- "celltype"

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$celltype)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

counts <- slide.seq[["Spatial"]]$counts
coords <- GetTissueCoordinates(slide.seq)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
slide.seq <- AddMetaData(slide.seq, metadata = )
library(readr)
RCTD <- readRDS("RCTD.rds")
slide.seq <- AddMetaData(slide.seq, metadata = RCTD@results$results_df)
p1 <- SpatialDimPlot(slide.seq, group.by = "first_type")
p2 <- SpatialDimPlot(slide.seq, group.by = "second_type")
p1 | p2

