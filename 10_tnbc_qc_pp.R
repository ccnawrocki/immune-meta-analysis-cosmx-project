rm(list = ls())
.rs.restartR(clean = T)

### Setup ----------------------------------------------------------------------

library(Matrix)
adata <- anndataR::read_h5ad("tnbc.h5ad", mode = "r+")
adata$obs[,c("slide", "corecode")] <- (adata$obs_names |> stringr::str_split(pattern = "_", simplify = T))[,2:3]

library(ggplot2)
ggplot() + 
  geom_point(data = adata$obs, mapping = aes(x = centroid_x, y = centroid_y, color = slide), shape = ".") + 
  coord_fixed()


### Fix spatial coordinates ----------------------------------------------------

# Need the min x and min y for each core
library(dplyr)
offsets <- adata$obs |> group_by(slide, corecode) |> summarise(min_x = min(centroid_x), min_y = min(centroid_y))
# dir.create("10_tnbc_qc_pp_outs")
# write.csv(x = offsets, file = "10_tnbc_qc_pp_outs/TNBC_offsets.csv", row.names = F)

# Creating "standard" coordinates 
x_standard <- adata$obs$centroid_x - as.numeric(plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$corecode), from = paste0(offsets$slide, offsets$corecode), to = offsets$min_x))
y_standard <- adata$obs$centroid_y - as.numeric(plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$corecode), from = paste0(offsets$slide, offsets$corecode), to = offsets$min_y))
plot(x_standard, y_standard, pch = ".")

# Need a key for where to send each core on our artificial canvas
offsets$origin_x <- case_when(
  offsets$slide %in% c("TNBC4") ~ 0, 
  offsets$slide %in% c("TNBC1", "TNBC3") ~ 3000, 
  offsets$slide %in% c("TNBC2") ~ 6000, 
)
offsets$origin_y <- case_when(
  
  offsets$slide == "TNBC4" & offsets$corecode == "coreD3" ~ 0, 
  offsets$slide == "TNBC4" & offsets$corecode == "coreE4" ~ 3000, 
  offsets$slide == "TNBC4" & offsets$corecode == "coreF2" ~ 6000, 
  offsets$slide == "TNBC4" & offsets$corecode == "coreF4" ~ 9000, 
  
  offsets$slide == "TNBC1" & offsets$corecode == "coreF1" ~ 0, 
  offsets$slide == "TNBC1" & offsets$corecode == "coreG1" ~ 3000, 
  offsets$slide == "TNBC3" & offsets$corecode == "coreD4" ~ 6000, 
  offsets$slide == "TNBC3" & offsets$corecode == "coreE2" ~ 9000, 
  
  offsets$slide == "TNBC2" & offsets$corecode == "coreC4" ~ 0, 
  offsets$slide == "TNBC2" & offsets$corecode == "coreD2" ~ 3000, 
  offsets$slide == "TNBC2" & offsets$corecode == "coreD4" ~ 6000, 
  
)

# Creating "global" coordinates
x_global <- x_standard + as.numeric(plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$corecode), from = paste0(offsets$slide, offsets$corecode), to = offsets$origin_x))
y_global <- y_standard + as.numeric(plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$corecode), from = paste0(offsets$slide, offsets$corecode), to = offsets$origin_y))
plot(x_global, y_global, pch = ".", asp = 1)

# Adding to the adata
adata$obs[,c("x_standard", "y_standard", "x_global", "y_global")] <- cbind(x_standard, y_standard, x_global, y_global)


### Add patient identifiers ----------------------------------------------------
tnbc_samples <- read.csv("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/sample-level-metadata/1kplex/core-map.csv", row.names = 1)
tnbc_samples$tma <- gsub(pattern = "TMA", replacement = "TNBC", x = tnbc_samples$tma)
tnbc_samples$corecode <- paste0("core", tnbc_samples$corecode)
adata$obs$patient <- plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$corecode), from = paste0(tnbc_samples$tma, tnbc_samples$corecode), to = tnbc_samples$patient)
adata$obs$patient <- gsub(pattern = "TN", replacement = "", x = adata$obs$patient) |> sprintf(fmt = "%03s")
adata$obs$patient <- paste0("bp", adata$obs$patient) # "bp" means "breast patient"


ggplot() + 
  geom_point(data = adata$obs, mapping = aes(x = x_global, y = y_global, color = patient), shape = ".") +
  scale_color_manual(values = InSituType::colorCellTypes(names = adata$obs$patient |> unique())) +
  coord_fixed() + 
  guides(color = guide_legend(override.aes = list(size = 5, shape = 16)))
ggplot() + 
  geom_point(data = adata$obs, mapping = aes(x = x_global, y = y_global, color = slide), shape = ".") +
  scale_color_manual(values = InSituType::colorCellTypes(names = adata$obs$slide |> unique())) +
  coord_fixed() + 
  guides(color = guide_legend(override.aes = list(size = 5, shape = 16)))

coremap <- openxlsx::read.xlsx("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/sample-level-metadata/1kplex/sgroi-tma-maps.xlsx", sheet = 1, startRow = 11, cols = 1:3)
adata$obs$corecode <- gsub(pattern = "core", replacement = "", x = adata$obs$corecode)
adata$obs$corenumber <- plyr::mapvalues(x = adata$obs$corecode, from = coremap$corecode, to = coremap$corenum)
adata$obs$corenumber <- paste0("core", sprintf("%02s", adata$obs$corenumber))

ggplot() + 
  geom_point(data = adata$obs, mapping = aes(x = x_global, y = y_global, color = corenumber), shape = ".") +
  scale_color_manual(values = InSituType::colorCellTypes(names = adata$obs$corenumber |> unique())) +
  coord_fixed() + 
  guides(color = guide_legend(override.aes = list(size = 5, shape = 16)))

adata$obs <- dplyr::relocate(adata$obs, corenumber, .after = slide)
adata$obs <- dplyr::relocate(adata$obs, corecode, .after = y_global)
adata$obs <- dplyr::rename(adata$obs, core = corecode)
glimpse(adata$obs)

anndataR::write_h5ad(object = adata, "tnbc.h5ad", mode = "w")


### Organizing the segmentation ------------------------------------------------
segfiles <- list.files(path = "02_proseg_outs", pattern = "shapes.parquet", full.names = T, recursive = T)
segfiles <- segfiles[grepl(pattern = "TNBC\\d+", x = segfiles)]
names(segfiles) <- stringr::str_split(string = segfiles, pattern = "/", simplify = T)[,2:3] |> 
  apply(MARGIN = 1, FUN = as.list) |> 
  lapply(FUN = do.call, what = paste0) |>
  gsub(pattern = "_outs", replacement = "", x = _)

seglist <- lapply(X = segfiles, FUN = arrow::read_parquet)
seglist <- lapply(X = seglist, FUN = sf::st_as_sf)
seglist$TNBC3coreE2 |> plot()

offsets$id <- paste0(offsets$slide, offsets$corecode)

translate_core <- function(core_seg, core_id) {
  local_translation <- core_seg
  sf::st_geometry(local_translation) <- sf::st_geometry(local_translation) + -offsets[offsets$id == core_id, c("min_x", "min_y"),] |> as.matrix() |> as.vector()
  global_translation <- local_translation
  sf::st_geometry(global_translation) <- sf::st_geometry(global_translation) + offsets[offsets$id == core_id, c("origin_x", "origin_y"),] |> as.matrix() |> as.vector()
  return(list("local" = local_translation, "global" = global_translation))
}

seglist <- purrr::map2(.x = seglist, .y = names(seglist), .f = translate_core)

local_seg <- sapply(X = seglist, FUN = "[[", "local", simplify = F)
local_seg <- dplyr::bind_rows(local_seg, .id = "sample") 
sample_ids <- stringr::str_split(string = local_seg$sample, pattern = "core", simplify = T)
local_seg$cell_ID <- paste0(local_seg$cell, "_", sample_ids[,1], "_", "core", sample_ids[,2])
local_seg <- select(local_seg, cell_ID, geometry)
(local_seg$cell_ID %in% adata$obs_names) |> mean() # 1

# global_seg <- sapply(X = seglist, FUN = "[[", "global", simplify = F)
# global_seg <- dplyr::bind_rows(global_seg, .id = "sample") 
# sample_ids <- stringr::str_split(string = global_seg$sample, pattern = "core", simplify = T)
# global_seg$cell_ID <- paste0(global_seg$cell, "_", sample_ids[,1], "_", "core", sample_ids[,2])
# global_seg <- select(global_seg, cell_ID, geometry)
# (global_seg$cell_ID %in% adata$obs_names) |> mean() # 1

sf::write_sf(local_seg, "tnbc_segmentation_local_coordinates.geojson", append = F)
# sf::write_sf(global_seg, "tnbc_segmentation_global_coordinates.geojson", append = F)


### Outliers -------------------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

library(Matrix)

adata <- anndataR::read_h5ad(path = "tnbc.h5ad", mode = "r+")
adata$obs$nCount_RNA <- rowSums(adata$X)
adata$obs$nFeature_RNA <- rowSums(adata$X > 0)
adata$obs$sample_ID <- paste(adata$obs$slide, adata$obs$core, sep = "_")

library(scuttle)

# First, identify cells with less than 15 counts
lowcounts <- (adata$obs$nCount_RNA < 15)

# Counts outliers
cts_outliers <- isOutlier(adata$obs$nCount_RNA, 
                          log = T, 
                          type = "lower", 
                          batch = adata$obs$sample_ID, 
                          nmads = 2, 
                          share.medians = F, share.mads = F
)
summary(cts_outliers)
#    Mode   FALSE    TRUE 
# logical  150545   10834 

# Some cores have thresholds that are much too low
(attributes(cts_outliers)$thresholds["lower",]) |> sort() |> knitr::kable()

# Thus, cells will be deemed low-quality if they have outlier counts OR <15 counts
adata$obs$counts_outlier <- (cts_outliers | lowcounts)
adata$obs |> dplyr::group_by(sample_ID) |> dplyr::summarise(prop_outlier = mean(counts_outlier))
# # A tibble: 11 × 2
# sample_ID prop_outlier
#   <chr>            <dbl>
# 1 TNBC1_F1        0.149 --> eval --> keep
# 2 TNBC1_G1        0.130 
# 3 TNBC2_C4        0.0754
# 4 TNBC2_D2        0.0461
# 5 TNBC2_D4        0.0781
# 6 TNBC3_D4        0.0537
# 7 TNBC3_E2        0.111 
# 8 TNBC4_D3        0.107 
# 9 TNBC4_E4        0.0721
# 10 TNBC4_F2       0.142 --> eval --> keep
# 11 TNBC4_F4       0.249 --> eval --> filter out

# I think that it makes sense generally for any core with >15% of cells that are
# low-quality to be excluded.

# Looking at the highly-affected cores 
tinyplot::plt(y_standard ~ x_standard | counts_outlier, 
              data = adata$obs |> dplyr::filter(sample_ID == "TNBC4_F4"),
              pch = 16, cex = 0.4,
              legend = legend(pt.cex = 2), 
              asp = 1
)
tinyplot::plt(y_standard ~ x_standard | counts_outlier, 
              data = adata$obs |> dplyr::filter(sample_ID == "TNBC1_F1"),
              pch = 16, cex = 0.4,
              legend = legend(pt.cex = 2), 
              asp = 1
)
tinyplot::plt(y_standard ~ x_standard | counts_outlier, 
              data = adata$obs |> dplyr::filter(sample_ID == "TNBC4_F2"),
              pch = 16, cex = 0.4,
              legend = legend(pt.cex = 2), 
              asp = 1
)

# Size outliers
vols <- adata$obs$volume
par(mar=c(4, 4, 4, 4))
hist(x = log(vols), breaks = 50)
ol <- isOutlier(log(vols), log = F, type = "both", nmads = 3)
(th <- attr(ol, "threshold")[1:2])
abline(v = th, col="blue")

adata$obs$volume_outlier <- as.vector(ol)


### QC -------------------------------------------------------------------------

tokeep <- !((adata$obs$counts_outlier) | (adata$obs$volume_outlier) | (adata$obs$sample_ID == "TNBC4_F4"))
mean(tokeep) # 0.8611839

# Need to use reticulate here
ad <- reticulate::import("anndata")
adata$write_h5ad(path = "tnbc.h5ad", mode = "w")

adata <- ad$read_h5ad("tnbc.h5ad")
adata <- adata[tokeep, ]
adata$write_h5ad("tnbc.h5ad")


### Processing -----------------------------------------------------------------

rm(list = ls())
.rs.restartR(clean = T)

## I will use reticulate to run scVI
# It should work already because we set it up back in script 04
ad <- reticulate::import("anndata")
SCVI <- reticulate::import("scvi")
sc <- reticulate::import("scanpy")
pd <- reticulate::import("pandas")

adata <- sc$read_h5ad("tnbc.h5ad")
adata$layers["counts"] <- adata$X

SCVI$model$LinearSCVI$setup_anndata(
  adata = adata, 
  layer = "counts", 
  batch_key = "sample_ID"
)

model <- SCVI$model$LinearSCVI(
  adata, 
  gene_likelihood = "nb", 
  n_latent = 30L
)

model$train(
  check_val_every_n_epoch = 4L,
  batch_size = 4096L,
  max_epochs = 400L,
  early_stopping = T,
  early_stopping_patience = 20L,
  early_stopping_monitor="elbo_validation",
  accelerator = "mps"
)
# RUNTIME: 06:19

elbo <- cbind(
  reticulate::py_to_r(model$history$elbo_train), 
  reticulate::py_to_r(model$history$elbo_validation), 
  epoch = 0:399
)

library(magrittr)
elbo$elbo_train %<>% as.numeric()
elbo$elbo_validation %<>% as.numeric()

library(ggplot2)
ggplot() + 
  geom_line(data = elbo |> tidyr::pivot_longer(cols = c("elbo_train", "elbo_validation"), values_to = "value", names_to = "type"), mapping = aes(x = epoch, y = value, color = type)) +
  scale_color_manual(values = c("elbo_train" = "steelblue", "elbo_validation" = "orange2")) +
  labs(title = "ELBO Plot for scVI Modeling") + 
  ggthemes::theme_hc()
# ggsave(filename = "10_tnbc_qc_pp_outs/ELBO_plot.pdf", height = 5, width = 5)
# model$save("10_tnbc_qc_pp_outs/scvi_model", overwrite = TRUE)

## Dimensionality reduction
SCVI_LATENT_KEY <- "X_scVI"
latent <- model$get_latent_representation()
adata$obsm[SCVI_LATENT_KEY] <- latent

adata$write_h5ad("tnbc.h5ad")

## UMAP
rm(list = ls())
.rs.restartR(clean = T)

library(RcppHNSW)

adata <- anndataR::read_h5ad("tnbc.h5ad", mode = "r+")
UM <- uwot::umap(X = adata$obsm$X_scVI, n_neighbors = 50, min_dist = 0.05, metric = "cosine", nn_method = "hnsw", spread = 1, 
                 fast_sgd = T, n_sgd_threads = 40, verbose = T)
adata$obsm$X_umap <- UM

library(BPCells)
plot_embedding(embedding = UM, source = adata$obs$sample_ID, rasterize = T, labels_discrete = F)

## Leiden clustering
knn <- knn_hnsw(adata$obsm$X_scVI, k = 50, metric = "cosine", ef = 1500)
knng <- Matrix::sparseMatrix(
  i = rep.int(seq_len(nrow(knn$idx)), ncol(knn$idx)),
  j = knn$idx |> as.vector(),
  x = 1
)
knng <- knng + Matrix::t(knng)  # Add transpose
knng <- (knng > 0) * 1  # Convert back to binary (0 or 1)
knng <- Matrix::drop0(knng)

iterative_clustering <- list()
for (i in seq(0.05, 2, 0.05)) {
  iterative_clustering[[as.character(i)]] <- cluster_graph_leiden(knng, resolution = i)
}
iterative_clustering <- dplyr::bind_cols(iterative_clustering)
colnames(iterative_clustering) <- paste("leiden_res", colnames(iterative_clustering), sep = "_")

# pdf(file = "10_tnbc_qc_pp_outs/umap_by_clusters.pdf", width = 6, height = 6)
for (nm in colnames(iterative_clustering)) {
  p <- plot_embedding(iterative_clustering[[nm]], UM, rasterize = T, labels_discrete = F) + 
    ggplot2::labs(title = nm)
  print(p)
}
# dev.off()

iterative_clustering <- as.data.frame(iterative_clustering)
rownames(iterative_clustering) <- adata$obs_names
adata$uns$iterative_clustering <- iterative_clustering

adata$obsp$knn <- knng

anndataR::write_h5ad(object = adata, path = "tnbc.h5ad", mode = "w")

## Normalization
rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("tnbc.h5ad", mode = "r+")

cts <- adata$layers$counts |>  as("CsparseMatrix")
scaling_factor <- 1000

norm_factors <- Matrix::Diagonal(x = scaling_factor/adata$obs$nCount_RNA, names = rownames(adata$layers$counts))
norm <- ((norm_factors %*% adata$layers$counts) |> log1p())/log(2)
adata$layers$lognorm <- norm

## Last small fix
adata$obsm$spatial <- dplyr::select(adata$obs, x_global, y_global) |> magrittr::set_rownames(NULL) |> magrittr::set_colnames(NULL) |> as.matrix()


### Final save -----------------------------------------------------------------
adata
# AnnData object with n_obs × n_vars = 138977 × 956
# obs: 'cell', 'original_cell_id', 'centroid_x', 'centroid_y', 'centroid_z', 'component', 'volume', 'surface_area', 'scale', 'region', 'slide', 'corenumber', 'x_standard', 'y_standard', 'x_global', 'y_global', 'core', 'patient', 'nCount_RNA', 'nFeature_RNA', 'sample_ID', 'counts_outlier', 'volume_outlier', '_scvi_batch', '_scvi_labels'
# uns: '_scvi_manager_uuid', '_scvi_uuid', 'iterative_clustering'
# obsm: 'X_scVI', 'X_umap', 'spatial'
# layers: 'counts', 'lognorm'
# obsp: 'knn'

anndataR::write_h5ad(object = adata, path = "tnbc.h5ad", mode = "w")

