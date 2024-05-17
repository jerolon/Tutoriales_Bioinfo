library(SPATA2)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(skimr)
#Leer el visium de la clase pasada
brain <- readRDS("visium_brain_clase12.RDS")
spata_brain <- asSPATA2(brain, spatial_method = "Visium",
                        assay_name = "SCT", image_name = "anterior1",
                        sample_name = "stxBrain")
#Misma grafica que con Seurat
SPATA2::plotSurface(spata_brain)
#Gene sets
skim(spata_brain@used_genesets)
printGeneSetOverview(spata_brain)

#Manipulacion de imagenes
plotSurface(spata_brain, pt_alpha = 0.2)

flipImage(spata_brain, axis = "horizontal") %>% 
plotSurface(pt_alpha = 0.5) +
  ggpLayerThemeCoords(unit = "px")

# Anotacion manual --------------------------------------------------------
#un circulo añadido como anotacion
theta <- seq(from= 0, to = 360, by = 20) * pi /180
radio <- 75
centro <- c(350,210)
circulo <- data.frame(x = radio*cos(theta) + centro[1],
                      y = radio * sin(theta) + centro[2])

spata_brain <- addImageAnnotation(spata_brain, 
                                  tags = "circulo_hechiso",
                                  area = list(outer = circulo))

layer <- ggpLayerImgAnnOutline(spata_brain, ids = "img_ann_1")

plotSurface(spata_brain, color_by = "nCount_Spatial", alpha = 0.2, display_image = TRUE) + layer


###Obtener info de los barcodes----
bcsp_dist <- getBarcodeSpotDistances(spata_brain) %>% 
  dplyr::filter(bc_origin != bc_destination)
#Cuntas micras hay en un pixel
pix2mic <- getPixelScaleFactor(spata_brain, unit = "um", switch = TRUE)

bcsp_dist_neighbors <- dplyr::filter(bcsp_dist, distance < 101*pix2mic) %>%
  group_by(bc_origin) %>% mutate(n_neighbors = n())

six_n <- bcsp_dist_neighbors %>% dplyr::filter(n_neighbors == 6) %>%
  ungroup() %>% slice_head(n=18)
two_n <- bcsp_dist_neighbors %>% dplyr::filter(n_neighbors == 2)
levels_loc <- c("Two", "Six", "Neigh", "None")

getCoordsDf(spata_brain) %>% mutate(
  category = case_when(
    barcodes %in% six_n$bc_origin ~"Six",
    barcodes %in% six_n$bc_destination ~ "Neigh",
    barcodes %in% two_n$bc_origin ~ "Two",
    barcodes %in% two_n$bc_destination ~ "Neigh",
    TRUE ~ "None"
  ), category = factor(x = category, levels = levels_loc)
) %>% 
  plotSurface(color_by = "category")

#Segmentacion basada en expresion
plotSurface(spata_brain, color_by = "seurat_clusters") + ggpLayerThemeCoords(unit = "px")

#Segmentacion basada en localizacion de los barcodes-----
getImageDims(spata_brain)
punto <- c(600, 0)
feature_df <- getCoordsDf(spata_brain) %>% mutate(distancia = 
  sqrt((punto[1] - x)^2 + (punto[2] - y)^2),
  grupo = as.factor(ntile(distancia, 10))
) %>% select(barcodes, grupo)

spata_brain <- addFeatures(spata_brain, feature_df = feature_df)
plotSurface(spata_brain, color_by = "grupo")
#Segmentación puede hacer análisis de expresion diferencial
spata_brain <- runDeAnalysis(spata_brain, across = "grupo")
plotDeaVolcano(spata_brain, across = "grupo", label_size = 7)
getDeaResultsDf(spata_brain, across = "grupo") %>% 
  arrange(desc(abs(avg_log2FC)))
genes_comp <- c("Lgr6","Penk","Meis2","Slc17a7","Olfm1","Ppp1r1b")
plotSurfaceComparison(spata_brain, color_by = genes_comp, nrow = 2)
#Comparacion de clusters contra histologia
plotBarchart(
  object = spata_brain, 
  grouping_variables = "grupo",
  across = "seurat_clusters"
)
#Anotacion a segmentacion-------
spata_brain <- imageAnnotationToSegmentation(spata_brain, 
                                             ids = "img_ann_1",
                                             segmentation_name = "circulo", 
                                             inside = "intra_circulo", outside = "extra_circulo",
                                             overwrite = TRUE
                                             )
plotSurface(spata_brain, color_by = "circulo", pt_clrp = "npg")
spata_brain <- runDeAnalysis(spata_brain, across = "circulo")
plotDeaVolcano(spata_brain, across = "circulo", label_size = 7)
genes_comp <- c("Dlgap1", "Stmn1", "Gnas", "Adora2a", "Gpr6", "Drd2")
plotSurfaceComparison(spata_brain, color_by = genes_comp, nrow = 2)

#Spatial trajectories
spata_brain <- 
  addSpatialTrajectory(
    object = spata_brain,
    id = "antero_posterior",
    start = c("500", "250") ,
    end = c("10", "250") ,
    width = "1mm",
    overwrite = TRUE
  )
spata_brain <- 
  addSpatialTrajectory(
    object = spata_brain,
    id = "dorso_ventral",
    start = c("400", "600") ,
    end = c("400", "0") ,
    width = "1mm",
    overwrite = TRUE
  )

#Solo buscamos genes que ya son espacialmente variables
#Si no se guardo de la clase pasada:
#spata_brain <- FindSpatiallyVariableFeatures(spata_brain, assay = "SCT", features = VariableFeatures(spata_brain)[1:1000],
#                                       selection.method = "moransi")

svf <- SVFInfo(spata_brain, method = "moransi", status = TRUE)
genes_toscreen <- dplyr::filter(svf, variable) %>% rownames
showModels()
modelos <- create_model_df(input = 10) %>% colnames()
spatial_traj_screen <- spatialTrajectoryScreening(spata_brain,
                                                  id = "antero_posterior",
                                                  variables = genes_toscreen)
plotOverview(spatial_traj_screen,
             label_vars = 4, 
             label_size = 5
)

library(patchwork)
trajectory <- 
  ggpLayerTrajectories(
    spata_brain, 
    ids = "antero_posterior",
    size = 1
  )

genes <- c("Nxph3", "Rgs2", "Fth1", "Shank1")
gene_colors <- color_vector(clrp = "npg", names = genes)
plist <- 
  imap(
    .x = gene_colors, 
    .f = function(color, gene){
      plotSurface(spata_brain, color_by = gene) + 
        scale_color_gradient(low = alpha("white", 0), high = color) + 
        trajectory
      
    })

wrap_plots(plist, ncol = 2)

plotTrajectoryLineplot(
  spata_brain,
  unit = "px",
  id = "antero_posterior",
  variables = genes, 
  smooth_se = TRUE, 
  clrp_adjust = gene_colors
)


shank_traj_df <- getStsDf(spata_brain,
                          id = "antero_posterior",
                          variables = "Shank1",
                          n_bins = 100, 
                          format = "wide")



# create data.frame with model variables
model_df <- 
  create_model_df(
    input = 1:100, # = 50
    var_order = "trajectory_order",
    model_subset = c("sharp_peak", "linear_descending", "sinus"), # only use models that contain these catchphrases
    model_remove = c("immediate", "two") # remove models that contain these catchphrases
  )

# show results 
model_df

# combine gene and model variables to a new data.frame
combined_df <- 
  cbind(
    shank_traj_df[,c("proj_length_binned", "trajectory_order","Shank1")],
    model_df[,c("sinus", "sinus_rev", "linear_descending","sharp_peak")]
  ) %>% 
  as_tibble()

# show results 
combined_df

p1 <- ggplot(data = combined_df, mapping = aes(x = trajectory_order, y = Shank1)) + 
  geom_point() + 
  geom_smooth(method = "loess", span = 0.25) + 
  theme_bw() + 
  labs(subtitle = "a) Shank1 - Course")

p2 <- ggplot(data = combined_df, mapping = aes(x = trajectory_order, y= sinus)) + 
  geom_point() + 
  geom_smooth(method = "loess", span = 0.25) + 
  theme_bw() + 
  labs(subtitle = "b) sinus")

p3 <- ggplot(data = combined_df, mapping = aes(x = trajectory_order, y = sinus_rev)) + 
  geom_point() + 
  geom_smooth(method = "loess", span = 0.25) + 
  theme_bw() + 
  labs(subtitle = "c) sinus_rev")

p4 <- ggplot(data = combined_df, mapping = aes(x = trajectory_order, y = linear_descending)) + 
  geom_point() + 
  geom_smooth(method = "loess", span = 0.25) + 
  theme_bw() + 
  labs(subtitle = "d) linear_descending")

#Correlation -------
p1 <- ggplot(data = combined_df, mapping = aes(x = Shank1, y = Shank1)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw() + 
  labs(subtitle = "a) Shank1 - Shank1")

p2 <- ggplot(data = combined_df, mapping = aes(x = sinus, y = Shank1)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw()  + 
  labs(subtitle = "b) Shank1  - sinus")

p3 <- ggplot(data = combined_df, mapping = aes(x = sinus_rev, y = Shank1)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw() +
  labs(subtitle = "c) Shank1  - sinus rev")

p4 <- ggplot(data = combined_df, mapping = aes(x = linear_descending, y= Shank1)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw() +
  labs(subtitle = "c) Shank1 - Linear descending")

p1 + p2 + p3 + p4


cor.test(x = combined_df[["Shank1"]], y = combined_df[["sinus"]])