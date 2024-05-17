library(Seurat)
library(matrixStats)
library(sctransform)
library(dplyr)
library(SeuratData)
library(ggplot2)
library(patchwork)
#######Visium
#InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
plot1 + plot2

brain@images$anterior1 %>% GetTissueCoordinates()
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(brain, interactive = TRUE)
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)

#Identificacion de variables espaciales
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

#Genes variables sin anotacion:
brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
    selection.method = "moransi")


svf <- SVFInfo(brain, method = "moransi", status = TRUE)

top.features <-  svf$rank %>% arrange(svf$rank) %>% head(n = 6) %>% rownames()
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
saveRDS(brain, "visium_brain_clase12.RDS")

##########Vizgen. Microscopy based spatial data
setwd("/mnt/atgc-d2/bioa03124/download.brainimagelibrary.org/63/05/63053886eb9f9771/660410/1231122825")
vizgen.obj <- ReadVizgen(data.dir = "merfish_output/cellpose_cyto2_nuclei/", type = "centroids")

cents <- CreateCentroids(vizgen.obj$centroids)
seg.data <- list(centroids = cents)
coords <- CreateFOV(coords = seg.data, 
                    type = "centroids", 
                    molecules = vizgen.obj$microns,
                    assay = "Vizgen")


viz.obj <- CreateSeuratObject(counts = vizgen.obj$transcripts, assay = "Vizgen", min.features=10)
coords <- subset(x = coords, cells = intersect(x = Cells(x = coords[["centroids"]]),
                                               y = Cells(x = viz.obj)))
viz.obj[["ace-dip-vat"]] <- coords
#The object is ready for visualization
p1 <- ImageFeaturePlot(viz.obj, features = "nCount_Vizgen")
p2 <- VlnPlot(viz.obj, features = c("nFeature_Vizgen", "nCount_Vizgen"), ncol = 2, pt.size = 0)
p1 + p2
#Intentar:
ImageFeaturePlot(viz.obj, features = "nCount_Vizgen", size = 5)

viz.obj@meta.data$log_umi <- log1p(viz.obj@meta.data$nCount_Vizgen)

#Te puedes brincar todo esto hasta el PCA con:
#viz.obj <- SCTransform(viz.obj, assay = "Vizgen", clip.range = c(-10, 10))
#En lugar de SCTrandform, llamamos a vst directamente
ref.model <- vst(viz.obj[["Vizgen"]]$counts,
                 cell_attr = viz.obj@meta.data,
                 verbosity = TRUE, n_cells = 5000)

residual.feature.mat <- get_residuals(vst_out = ref.model, 
                                      umi = viz.obj[["Vizgen"]]$counts,
                                      verbosity = 2, cell_attr = viz.obj@meta.data)
ref.model$y <- residual.feature.mat
ref.model$gene_attr$residual_mean <- NA_real_
ref.model$gene_attr$residual_variance <- NA_real_
ref.model$gene_attr$residual_mean <- rowMeans2(x = ref.model$y)
ref.model$gene_attr$residual_variance <- rowVars(x = ref.model$y)
ref.model

scale.data <- ref.model$y
ref.model$y <- ScaleData(scale.data, features = NULL, model.use = "linear", use.umi = FALSE, 
                         do.scale = TRUE, do.center = TRUE, scale.max = Inf, 
                         block.size = 750, min.cells.to.block = 3000, verbose = TRUE)

assay.out <- CreateSCTAssayObject(viz.obj[["Vizgen"]]$counts, scale.data = ref.model$y)
Misc(object = assay.out, slot = "vst.out") <- ref.model
viz.obj[["SCT"]] <- assay.out
DefaultAssay(object = viz.obj) <- "SCT"


viz.obj <- RunPCA(viz.obj, npcs = 30, features = rownames(viz.obj))
viz.obj <- RunUMAP(viz.obj, dims = 1:30)
viz.obj <- FindNeighbors(viz.obj, reduction = "pca", dims = 1:30)
viz.obj <- FindClusters(viz.obj, resolution = 0.3)                                                      


DimPlot(viz.obj, reduction = "umap")
ImageDimPlot(viz.obj, axes = TRUE, cols = "polychrome")

#Los genes variables estan contenidos en el modelo de SCT
ref.model$gene_attr %>% filter(detection_rate > 0.5) %>% arrange(desc(residual_variance)) %>% head(n = 20)
ImageFeaturePlot(viz.obj, features = c("Sox2ot", "Tcf7l2", "Enpp2", "Sox10", "Nrgn", "Opalin"))
#Clusters especificos
ImageDimPlot(viz.obj, cols = "red", cells = WhichCells(viz.obj, idents = 30))
#Las coordenadas de estos experimentos se encuentran en 
viz.obj@images$ace.dip.vat$centroids 
library(sf)
library(future)
options(future.globals.maxSize = 8000 * 1024^2)
cropped.coords <- Crop(viz.obj@images$ace.dip.vat, x = c(3132.945,3237.376), y = c(3196.334,3302.879), coords = "plot")
viz.obj[["zoom"]] <- cropped.coords
DefaultBoundary(viz.obj[["zoom"]]) <- "centroids"
ImageDimPlot(viz.obj, axes = TRUE, cols = "polychrome", fov = "zoom", size = 2, molecules = c("Gad1", "Sst", "Npy2r", "Pvalb", "Nrn1"), border.size = 0.1, border.color = "white")