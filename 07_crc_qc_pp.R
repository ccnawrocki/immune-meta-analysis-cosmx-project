rm(list = ls())
.rs.restartR(clean = T)

### Setup ----------------------------------------------------------------------

library(Matrix)
adata <- anndataR::read_h5ad("crc.h5ad", mode = "r+")
adata$obs[,c("slide", "core")] <- (adata$obs_names |> stringr::str_split(pattern = "_", simplify = T))[,2:3]
adata$obs$core <- gsub(pattern = "core", replacement = "", x = adata$obs$core)

library(ggplot2)
ggplot() + 
  geom_point(data = adata$obs, mapping = aes(x = centroid_x, y = centroid_y, color = slide), shape = ".") + 
  coord_fixed()


### Fix spatial coordinates ----------------------------------------------------

# Need the min x and min y for each core
library(dplyr)
offsets <- adata$obs |> group_by(slide, core) |> summarise(min_x = min(centroid_x), min_y = min(centroid_y))
# dir.create("07_crc_qc_pp_outs")
# write.csv(x = offsets, file = "07_crc_qc_pp_outs/CRC_offsets.csv", row.names = F)

# Creating "standard" coordinates 
x_standard <- adata$obs$centroid_x - as.numeric(plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$core), from = paste0(offsets$slide, offsets$core), to = offsets$min_x))
y_standard <- adata$obs$centroid_y - as.numeric(plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$core), from = paste0(offsets$slide, offsets$core), to = offsets$min_y))
plot(x_standard, y_standard, pch = ".")

# Need a key for where to send each core on our artificial canvas
offsets$origin_x <- case_when(
  (offsets$slide == "CRC1" & offsets$core %in% c("B1", "B2", "E2", "F3", "C4")) ~ 0, 
  (offsets$slide == "CRC2" & offsets$core %in% c("A4", "C3", "D1")) ~ 4000, 
  (offsets$slide == "CRC2" & offsets$core %in% c("D2", "D3", "D4", "E2", "E4")) ~ 8000, 
)
offsets$origin_y <- c(seq(0, 16000, 4000), seq(0, 8000, 4000), seq(0, 16000, 4000))

# Creating "global" coordinates
x_global <- x_standard + as.numeric(plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$core), from = paste0(offsets$slide, offsets$core), to = offsets$origin_x))
y_global <- y_standard + as.numeric(plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$core), from = paste0(offsets$slide, offsets$core), to = offsets$origin_y))
plot(x_global, y_global, pch = ".", asp = 1)

# Adding to the adata
adata$obs[,c("x_standard", "y_standard", "x_global", "y_global")] <- cbind(x_standard, y_standard, x_global, y_global)


### Add patient identifiers ----------------------------------------------------
crc_samples <- openxlsx::read.xlsx(xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-crc-data/sample_map.xlsx", sheet = "core map")
crc_samples <- crc_samples[crc_samples$tma |> grepl(pattern = "prim"),]
crc_samples$tma <- gsub(pattern = "prim", replacement = "CRC", x = crc_samples$tma)
adata$obs$patient <- plyr::mapvalues(x = paste0(adata$obs$slide, adata$obs$core), from = paste0(crc_samples$tma, crc_samples$core), to = crc_samples$patient)
adata$obs$patient <- gsub(pattern = "p", replacement = "", x = adata$obs$patient) |> sprintf(fmt = "%03s")
adata$obs$patient <- paste0("cp", adata$obs$patient) # "cp" means "colon patient"

# renv::install("NanoString-Biostats/InSituType")
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

anndataR::write_h5ad(object = adata, "crc.h5ad", mode = "w")


### Organizing the segmentation ------------------------------------------------
segfiles <- list.files(path = "02_proseg_outs", pattern = "shapes.parquet", full.names = T, recursive = T)
segfiles <- segfiles[grepl(pattern = "CRC\\d+", x = segfiles)]
names(segfiles) <- stringr::str_split(string = segfiles, pattern = "/", simplify = T)[,2:3] |> 
  apply(MARGIN = 1, FUN = as.list) |> 
  lapply(FUN = do.call, what = paste0) |>
  gsub(pattern = "_outs", replacement = "", x = _)

seglist <- lapply(X = segfiles, FUN = arrow::read_parquet)
seglist <- lapply(X = seglist, FUN = sf::st_as_sf)
seglist$CRC2coreE4 |> plot()

offsets$id <- paste0(offsets$slide, "core", offsets$core)

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

sf::write_sf(local_seg, "crc_segmentation_local_coordinates.geojson", append = F)
# sf::write_sf(global_seg, "crc_segmentation_global_coordinates.geojson", append = F)


### Outliers -------------------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

library(Matrix)

adata <- anndataR::read_h5ad(path = "crc.h5ad", mode = "r+")
adata$obs$nCount_RNA <- rowSums(adata$X)
adata$obs$nFeature_RNA <- rowSums(adata$X > 0)
adata$obs$sample_ID <- paste0(adata$obs$slide, "_core", adata$obs$core)

# First, identify cells with less than 30 counts
lowcounts <- (adata$obs$nCount_RNA < 30)

library(scuttle)

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
# logical  410273   21732 

# Some cores have thresholds that are much too low
(attributes(cts_outliers)$thresholds["lower",]) |> sort() |> knitr::kable()

# Thus, cells will be deemed low-quality if they have outlier counts OR <30 counts
adata$obs$counts_outlier <- (cts_outliers | lowcounts)
adata$obs |> dplyr::group_by(sample_ID) |> dplyr::summarise(prop_outlier = mean(counts_outlier))
# # A tibble: 13 × 2
#   sample_ID   prop_outlier
#   <chr>              <dbl>
# 1 CRC1_coreB1       0.0600
# 2 CRC1_coreB2       0.0655
# 3 CRC1_coreC4       0.0955 --> eval --> seems okay to keep
# 4 CRC1_coreE2       0.0601
# 5 CRC1_coreF3       0.0491
# 6 CRC2_coreA4       0.0495
# 7 CRC2_coreC3       0.0292
# 8 CRC2_coreD1       0.185 --> eval --> sadly, I think it must be filtered out
# 9 CRC2_coreD2       0.0533
# 10 CRC2_coreD3      0.0501
# 11 CRC2_coreD4      0.0899
# 12 CRC2_coreE2      0.0512
# 13 CRC2_coreE4      0.107 --> eval --> seems okay to keep

# I think that it makes sense generally for any core with >15% of cells that are
# low-quality to be excluded.

# Looking at the highly-affected cores 
tinyplot::plt(y_standard ~ x_standard | counts_outlier, 
              data = adata$obs |> dplyr::filter(sample_ID == "CRC2_coreD1"),
              pch = 16, cex = 0.4,
              legend = legend(pt.cex = 2), 
              asp = 1
)
tinyplot::plt(y_standard ~ x_standard | counts_outlier, 
              data = adata$obs |> dplyr::filter(sample_ID == "CRC2_coreE4"),
              pch = 16, cex = 0.4,
              legend = legend(pt.cex = 2), 
              asp = 1
)
tinyplot::plt(y_standard ~ x_standard | counts_outlier, 
              data = adata$obs |> dplyr::filter(sample_ID == "CRC1_coreC4"),
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

tokeep <- !((adata$obs$counts_outlier) | (adata$obs$volume_outlier) | (adata$obs$sample_ID == "CRC2_coreD1"))
mean(tokeep) # 0.8459416

# Need to use reticulate here
ad <- reticulate::import("anndata")
adata$write_h5ad(path = "crc.h5ad", mode = "w")

adata <- ad$read_h5ad("crc.h5ad")
adata <- adata[tokeep, ]
adata$write_h5ad("crc.h5ad")


### Processing -----------------------------------------------------------------

## Finding HVGs
rm(list = ls())
.rs.restartR(clean = T)

# renv::install("nebula")

# NOTE: nebula has some issues with serialization in RStudio. One may have to 
# make a fresh environment as follows and run this piece in that environment:
# renv::upgrade()
# renv::init(bare = T)
# renv::install(packages = list("bioc::SingleCellExperiment", "scverse/anndataR@084a187", "nebula@1.5.6"))
# renv::install("bioc::rhdf5")
# renv::snapshot()

adata <- anndataR::read_h5ad(path = "crc.h5ad", mode = "r+")

# # Setting up the data 
# cts <- adata$X
# dimnames(cts) <- list(adata$obs_names, adata$var_names)
# cts <- Matrix::t(cts) |> as("CsparseMatrix")
# 
# meta <- adata$obs
# idx <- meta |> dplyr::arrange(sample_ID) |> rownames()
# cts <- cts[,idx]
# meta <- meta[idx,]
# mm <- model.matrix(~1, data = meta)
# eff <- meta$nCount_RNA
# ids <- meta$sample_ID
# 
# # Fitting the nebula model
# nfit <- nebula::nebula(count = cts, id = ids, pred = mm, offset = eff,
#                        covariance = F, cpc = 0.001, ncore = 8)
# saveRDS(object = nfit, file = "07_crc_qc_pp_outs/hvg_nfit.RDS")
nfit <- readRDS(file = "07_crc_qc_pp_outs/hvg_nfit.RDS")

# Plotting the results
tokeep <- nfit$convergence >= -10
hvg_data <- data.frame(mean_expr = nfit$summary$`logFC_(Intercept)`,
                       se_expr = nfit$summary$`se_(Intercept)`,
                       overdispersion = nfit$overdispersion$Cell, 
                       gene = nfit$summary$gene)
hvg_data <- hvg_data[tokeep,]
hvgs <- dplyr::arrange(hvg_data, desc(overdispersion))[1:2500, "gene"]

plot(x = hvg_data$mean_expr, hvg_data$overdispersion, pch = 16, cex = 0.5,
     col = ifelse(hvg_data$gene %in% hvgs, yes = "red", no = "black"))

plot(x = hvg_data$mean_expr, y = sqrt(hvg_data$se_expr), pch = 16, cex = 0.5)
fitted <- loess(sqrt(hvg_data$se_expr) ~ hvg_data$mean_expr, span = 0.3)
lines(x = fitted$x[order(fitted$x)], y = fitted$fitted[order(fitted$x)], col = "red", lwd = 4)

# Adding these results to the AnnData
adata$var[["overdispersion"]] <- NA
adata$var[hvg_data$gene, "overdispersion"] <- hvg_data$overdispersion
adata$var[["log_cpc"]] <- NA
adata$var[hvg_data$gene, "log_cpc"] <- hvg_data$mean_expr
adata$var[["log_cpc_se"]] <- NA
adata$var[hvg_data$gene, "log_cpc_se"] <- hvg_data$se_expr

# Specifying HVGs in the AnnData
hvgs <- (adata$var |> dplyr::filter(log_cpc > quantile(log_cpc, 0.025, na.rm = T)) |> dplyr::arrange(desc(overdispersion)) |> rownames())[1:5500]
adata$var[["highly_variable"]] <- (adata$var_names %in% hvgs)
plot(x = adata$var$log_cpc, y = adata$var$overdispersion, pch = 16, cex = 0.5,
     col = ifelse(adata$var$highly_variable, yes = "red", no = "black"))

# Even after subsetting to HVGs, every cell still has >=15 counts.
sum((magrittr::set_colnames(adata$X, rownames(adata$var))[,hvgs] |> Matrix::rowSums()) < 15) # 0

anndataR::write_h5ad(adata, path = "crc.h5ad", mode = "w")

## I will use reticulate to run scVI
rm(list = ls())
.rs.restartR(clean = T)

# It should work already because we set it up back in script 04
ad <- reticulate::import("anndata")
SCVI <- reticulate::import("scvi")
sc <- reticulate::import("scanpy")
pd <- reticulate::import("pandas")

adata <- sc$read_h5ad("crc.h5ad")
adata$layers["counts"] <- adata$X

adata_hvgs <- adata[,adata$var$highly_variable]$copy()

SCVI$model$LinearSCVI$setup_anndata(
  adata = adata_hvgs, 
  layer = "counts", 
  batch_key = "sample_ID"
)

model <- SCVI$model$LinearSCVI(
  adata_hvgs, 
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
# RUNTIME: 1:15:26

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
# ggsave(filename = "07_crc_qc_pp_outs/ELBO_plot.pdf", height = 5, width = 5)
# model$save("07_crc_qc_pp_outs/scvi_model", overwrite = TRUE)

## Dimensionality reduction
SCVI_LATENT_KEY <- "X_scVI"
latent <- model$get_latent_representation()
adata$obsm[SCVI_LATENT_KEY] <- latent

adata$write_h5ad("crc.h5ad")

## UMAP
rm(list = ls())
.rs.restartR(clean = T)

library(RcppHNSW)
adata <- anndataR::read_h5ad("crc.h5ad", mode = "r+")
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

# pdf(file = "07_crc_qc_pp_outs/umap_by_clusters.pdf", width = 6, height = 6)
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

anndataR::write_h5ad(object = adata, path = "crc.h5ad", mode = "w")

## Normalization
rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("crc.h5ad", mode = "r+")

cts <- adata$layers$counts |>  as("CsparseMatrix")
scaling_factor <- 1000

norm_factors <- Matrix::Diagonal(x = scaling_factor/adata$obs$nCount_RNA, names = rownames(adata$layers$counts))
norm <- ((norm_factors %*% adata$layers$counts) |> log1p())/log(2)
adata$layers$lognorm <- norm

## Last small fix
adata$obsm$spatial <- dplyr::select(adata$obs, x_global, y_global) |> magrittr::set_rownames(NULL) |> magrittr::set_colnames(NULL) |> as.matrix()


### Final save -----------------------------------------------------------------
adata
# AnnData object with n_obs × n_vars = 365451 × 6180
# obs: 'cell', 'original_cell_id', 'centroid_x', 'centroid_y', 'centroid_z', 'component', 'volume', 'surface_area', 'scale', 'region', 'slide', 'core', 'x_standard', 'y_standard', 'x_global', 'y_global', 'patient', 'nCount_RNA', 'nFeature_RNA', 'sample_ID', 'counts_outlier', 'volume_outlier'
# var: 'overdispersion', 'log_cpc', 'log_cpc_se', 'highly_variable'
# uns: 'iterative_clustering'
# obsm: 'X_scVI', 'X_umap', 'spatial'
# layers: 'counts', 'lognorm'
# obsp: 'knn'

anndataR::write_h5ad(object = adata, path = "crc.h5ad", mode = "w")

