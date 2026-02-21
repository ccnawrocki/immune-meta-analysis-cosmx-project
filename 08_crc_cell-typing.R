rm(list = ls())
.rs.restartR(clean = T)

library(Matrix)
library(BPCells)
library(ggplot2)

# dir.create(path = "08_crc_cell-typing_outs")

## scVI normalized expression  -------------------------------------------------
# Note that this is how one can get it, but I am not going to use it.
SCVI <- reticulate::import("scvi")
ad <- reticulate::import("anndata")

adata <- ad$read_h5ad("crc.h5ad")

# model <- SCVI$model$SCVI$load("07_crc_qc_pp_outs/scvi_model")
# genes_oi <- c("EPCAM",
#               "COL1A1",
#               "CD68",
#               "CD74",
#               "PECAM1", 
#               "CD2", 
#               "CD4", 
#               "CD8A", 
#               "CD3E", 
#               "JCHAIN",
#               "IGHG1/2",
#               "NKG7"
#               )
# adata <- adata[,adata$var$highly_variable]$copy()
# denoised <- model$get_normalized_expression(adata = adata, library_size = 1000, gene_list = genes_oi)


## Plotting some canonical markers ---------------------------------------------
# pdf(file = "08_crc_cell-typing_outs/umap_by_canonical_markers.pdf", width = 12, height = 10)
plot_embedding(
  source = t(adata$layers["lognorm"]) |> magrittr::set_rownames(rownames(adata$var)) |> as("CsparseMatrix"),
  embedding = adata$obsm["X_umap"],
  features = c("EPCAM",
               "COL1A1",
               "CD68",
               "CD74",
               "PECAM1", 
               "CD2"
  ),
  rasterize = T, 
  colors_continuous = viridis::viridis(n = 71), 
  quantile_range = c(0.01, 0.99)
) & 
  labs(x = "UM1", y = "UM2")

plot_embedding(
  source = t(adata$layers["lognorm"]) |> magrittr::set_rownames(rownames(adata$var)) |> as("CsparseMatrix"),
  embedding = adata$obsm["X_umap"],
  features = c("CD4", 
               "CD8A", 
               "CD3E", 
               "JCHAIN",
               "IGHG1/2",
               "NKG7"
  ), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71), 
  quantile_range = c(0.01, 0.99)
) & 
  labs(x = "UM1", y = "UM2")

plot_embedding(
  source = t(adata$layers["lognorm"]) |> magrittr::set_rownames(rownames(adata$var)) |> as("CsparseMatrix"),
  embedding = adata$obsm["X_umap"],
  features = c("PIGR", 
               "ACTA2", 
               "MYL9", 
               "ELANE", 
               "KRT19", 
               "ANXA4"
  ), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71), 
  quantile_range = c(0.01, 0.99)
) & 
  labs(x = "UM1", y = "UM2")
# dev.off()


## Plotting some other variables -----------------------------------------------
pdf(file = "08_crc_cell-typing_outs/umap_by_other_variables.pdf", width = 6, height = 5)
plot_embedding(
  source = adata$obs$slide,
  embedding = adata$obsm["X_umap"],
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "CRC: slide")
plot_embedding(
  source = adata$obs$patient,
  embedding = adata$obsm["X_umap"],
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "CRC: patient")
plot_embedding(
  source = adata$obs$sample_ID,
  embedding = adata$obsm["X_umap"],
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "CRC: sample ID") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, shape = 15)))
dev.off()


## Clusters that I think look good ---------------------------------------------
# pdf(file = "08_crc_cell-typing_outs/umap_by_leiden_clusters.pdf", width = 6, height = 5)
plot_embedding(
  source = adata$uns$iterative_clustering$leiden_res_0.7,
  embedding = adata$obsm["X_umap"],
  rasterize = T, labels_discrete = T
) + 
  labs(x = "UM1", y = "UM2", title = "CRC: leiden clustering (res = 0.7)") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, shape = 15)))
# dev.off()


## Autofluorescence score ------------------------------------------------------
adata <- anndataR::read_h5ad("crc.h5ad", mode = "r+")
probes <- readRDS(file = url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/refs/heads/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))
probes <- probes$Hs_6k
probes[probes$barcode |> grepl(pattern = ".*BB.*BB.*BB.*") |> which(),]
bluetargets <- probes$gene[probes$barcode |> grepl(pattern = ".*BB.*BB.*BB.*") |> which()] |> intersect(adata$var |> rownames())

# Calculating the meta score can be memory-intensive, so we will just visualize for now.
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  groups = adata$uns$iterative_clustering$leiden_res_0.7, 
                  features = bluetargets
) + 
  scale_color_viridis_c(limits = c(-1, 2), oob = scales::squish, option = "B", direction = -1) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_blank()) + 
  labs(title = "Blue Targets", y = "Cluster", x = "Gene")

# Tough to say, but cluster 11 is probably just RBCs.
# Clusters 1, 8, and 9 are all epithelial. 1 and 9 are malignant.


## Finding cluster markers with presto GLMM ------------------------------------
# Pseudo-bulking
library(data.table)
idx <- !(adata$uns$iterative_clustering$leiden_res_0.7 %in% c("11", "1", "9"))
psb <- presto::collapse_counts(counts_mat = (Matrix::t(adata$layers$counts[idx,]) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names)),
                               meta_data = data.frame("sample_ID" = adata$obs$sample_ID[idx],
                                                      "cluster" = adata$uns$iterative_clustering$leiden_res_0.7[idx, drop = T]),
                               get_norm = F,
                               varnames = c("cluster", "sample_ID"), 
                               min_cells_per_group = 20
)
psb$meta_data$logUMI <- log(psb$counts_mat |> colSums())

# Modeling
library(presto)
library(lme4)
library(purrr)
library(dplyr)
presto_res <- presto.presto(
  y ~ 1 + (1|cluster) + (1|sample_ID) + (1|cluster:sample_ID) + offset(logUMI),
  psb$meta_data,
  psb$counts_mat,
  size_varname = "logUMI",
  effects_cov = c("cluster"),
  ncore = 8,
  min_sigma = 0.05,
  family = "poisson",
  nsim = 1000
)
saveRDS(object = presto_res, file = "08_crc_cell-typing_outs/presto_model_psb.RDS")
presto_res <- readRDS(file = "08_crc_cell-typing_outs/presto_model_psb.RDS")

# Visulizing results
contrasts_mat <- make_contrast.presto(presto_res, "cluster")
effects_marginal <- contrasts.presto(presto_res, contrasts_mat, one_tailed = T) |> 
  dplyr::mutate(cluster = contrast) |> 
  dplyr::mutate(
    logFC = sign(beta) * log2(exp(abs(beta))), # convert stats to log2
    SD = log2(exp(sigma)),
    zscore = logFC / SD
  ) |> 
  dplyr::select(cluster, feature, logFC, SD, zscore, pvalue) |> 
  dplyr::arrange(pvalue)
effects_marginal$fdr <- p.adjust(effects_marginal$pvalue, method = "BH")

library(ggplot2)
p <- effects_marginal |> 
  ggplot() + 
  scattermore::geom_scattermore(mapping = aes(x = logFC, y = -log10(fdr)), pointsize = 3) + 
  geom_hline(yintercept = -log10(0.01)) + 
  geom_vline(xintercept = 1.5) + 
  ggrepel::geom_text_repel(data = effects_marginal |> filter(fdr < 0.01 & logFC > 1.5), 
                           mapping = aes(x = logFC, y = -log10(fdr), label = feature), 
                           color = "red", size = 3, max.overlaps = 25, box.padding = 0.5) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 10, color = "black"), 
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(size = 12, color = "black"), 
        axis.title.y = element_text(size = 12, color = "black"), 
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5)
  ) + 
  labs(title = "Leiden Cluster Markers", x = "log2FC", y = "-log10(Adj. P-value)") +
  facet_wrap(.~cluster, ncol = 3, scales = "free")
p
# ggsave(filename = "08_crc_cell-typing_outs/cluster_markers_by_presto_glmm_volcano_plots.pdf", height = 10, width = 12)

top_marks <- dplyr::group_by(effects_marginal, cluster) |> 
  dplyr::filter(fdr < 0.01 & logFC > 1.5) |>
  dplyr::mutate(rank = order(zscore, decreasing = T)) |> 
  dplyr::group_by(cluster) |> 
  dplyr::top_n(n = -10, wt = rank)
top_marks <- tidyr::pivot_wider(data = top_marks, id_cols = rank, names_from = cluster, values_from = feature) |> 
  dplyr::arrange(rank)
top_marks

BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm[idx,]) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  groups = adata$uns$iterative_clustering$leiden_res_0.7[idx, drop = T], 
                  group_order = colnames(top_marks)[-1],
                  features = top_marks[,-1] |> as.matrix() |> as.vector() |> unique()
) + 
  scale_color_viridis_c(limits = c(-1, 2), oob = scales::squish, option = "B", direction = -1) + 
  ggthemes::theme_hc() +
  theme(axis.text.x = element_text(size = 6)) +
  labs(title = "DE Markers", y = "Cluster", x = "Gene")
# ggsave(filename = "08_crc_cell-typing_outs/cluster_markers_by_presto_glmm_bubble_plot.pdf", height = 5, width = 16)

# Looking at some cores: 
d <- adata$obs
d$cluster <- adata$uns$iterative_clustering$leiden_res_0.7

# pdf(file = "08_crc_cell-typing_outs/clusters_in_space.pdf", height = 10, width = 10)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "CRC1_coreF3"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "CRC2_coreD3"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "CRC2_coreA4"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "CRC2_coreC3"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "CRC1_coreE2"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "CRC2_coreD2"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
# dev.off()


## Annotating ------------------------------------------------------------------
adata$obs$celltype_level1 <- dplyr::case_when(
  adata$uns$iterative_clustering$leiden_res_0.7 == 1 ~ "Epithelial",
  adata$uns$iterative_clustering$leiden_res_0.7 == 2 ~ "Immune.Myeloid", 
  adata$uns$iterative_clustering$leiden_res_0.7 == 3 ~ "Fibroblast", 
  adata$uns$iterative_clustering$leiden_res_0.7 == 4 ~ "Endothelial",
  adata$uns$iterative_clustering$leiden_res_0.7 == 5 ~ "B.cell.Lineage",
  adata$uns$iterative_clustering$leiden_res_0.7 == 6 ~ "Neutrophil",
  adata$uns$iterative_clustering$leiden_res_0.7 == 7 ~ "T.and.NK",
  adata$uns$iterative_clustering$leiden_res_0.7 == 8 ~ "Epithelial",
  adata$uns$iterative_clustering$leiden_res_0.7 == 9 ~ "Epithelial",
  adata$uns$iterative_clustering$leiden_res_0.7 == 10 ~ "Myocyte",
  adata$uns$iterative_clustering$leiden_res_0.7 == 11 ~ "Erythrocyte", 
  adata$uns$iterative_clustering$leiden_res_0.7 == 12 ~ "Neural"
)

# pdf(file = "08_crc_cell-typing_outs/umap_by_celltype_level1.pdf", width = 6, height = 5)
plot_embedding(
  source = adata$obs$celltype_level1,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "CRC: cell types (level 1)") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, shape = 15)))
# dev.off()

# Saving
anndataR::write_h5ad(object = adata, path = "crc.h5ad", mode = "w")


## Validation ------------------------------------------------------------------

# I want to check on cluster 6...
SCVI <- reticulate::import("scvi")
ad <- reticulate::import("anndata")

adata <- ad$read_h5ad("crc.h5ad")

plot_embedding(
  source = adata$obs[,c("nCount_RNA", "nFeature_RNA"), drop = F] |> log2() |> t(),
  embedding = adata$obsm["X_umap"],
  rasterize = T, 
  features = c("nCount_RNA", "nFeature_RNA"),
  colors_continuous = viridis::viridis(n = 71), 
  quantile_range = c(0.01, 0.99)
) & 
  labs(x = "UM1", y = "UM2")

d <- adata$obs
d$neutrophil <- ifelse(test = d$celltype_level1 == "Neutrophil", yes = T, no = F)
tinyplot::plt(y_global ~ x_global | neutrophil,
              data = d |> dplyr::filter(sample_ID == "CRC2_coreC3"), 
              col = BPCells::discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)

# I also want to check on cluster 11...
d <- adata$obs
tinyplot::plt(y_global ~ x_global | celltype_level1, 
              data = d |> dplyr::filter(sample_ID == "CRC2_coreE2"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)

# The fact that cluster 6 has very low counts is evidence that it is likely 
# neutrophils. Cluster 11 does indeed look like RBCs on the image.


## InSituType ------------------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("crc.h5ad", mode = "r+")

# Tricky part: need the mean negative probe expression in each cell
library(DBI)
library(duckdb)
library(dplyr)

crc_stow <- "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-crc-data"
txpaths <- list.files(pattern = "tx_file.csv.gz", 
                      recursive = T, 
                      path = crc_stow)
txpaths <- txpaths[grepl(pattern = "(ColonPrimaryTMA1)|(PrimaryColonTMA2)", x = txpaths)]

# Getting all the negative probe locations for our samples
negprobe_list <- list()
for (p in txpaths) {
  con <- dbConnect(duckdb(), dbdir = ":memory:")
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  P <- file.path(crc_stow, p)
  tx <- dbGetQuery(
    con, 
    sprintf("
            FROM read_csv('%s')
            WHERE target IN ('Negative16', 'Negative8', 'Negative14', 'Negative11', 'Negative5', 'Negative20', 'Negative7', 'Negative3', 'Negative4', 'Negative10', 'Negative12', 'Negative13', 'Negative17', 'Negative2', 'Negative15', 'Negative6', 'Negative18', 'Negative1', 'Negative19', 'Negative9')
            ",
            P
    )
  )
  dbDisconnect(con, shutdown = T)
  negprobe_list[[p]] <- tx
}
names(negprobe_list) <- stringr::str_split(string = names(negprobe_list), pattern = "/", simplify = T)[,1]
negprobes <- dplyr::bind_rows(negprobe_list, .id = "slide")
negprobes <- dplyr::select(negprobes, slide, fov, target, x_global_px, y_global_px)
crc_samples <- openxlsx::read.xlsx(file.path(crc_stow, "sample_map.xlsx"), sheet = "fov map")
negprobes$core <- plyr::mapvalues(x = paste(negprobes$slide, negprobes$fov, sep = "_"), from = paste(crc_samples$flowcell, crc_samples$fov, sep = "_"), to = crc_samples$core)
negprobes$slide <- ifelse(test = grepl(pattern = "1", x = negprobes$slide), yes = "CRC1", no = "CRC2")
negprobes$sample_ID <- paste0(negprobes$slide, "_core", negprobes$core)
negprobes <- negprobes[negprobes$sample_ID %in% unique(adata$obs$sample_ID),]

# Converting to microns
negprobes$x_global_um <- negprobes$x_global_px*0.12028
negprobes$y_global_um <- negprobes$y_global_px*0.12028

# Creating separate sf objects for each sample
negprobes <- split(x = negprobes, f = negprobes$sample_ID)
negprobes <- lapply(X = negprobes, FUN = sf::st_as_sf, coords = c("x_global_um", "y_global_um"))
negprobes[["CRC2_coreE2"]] |> plot() # Nice

# Transforming the coordinates to standard scale, as we did with the segmentations.
offsets <- read.csv(file = "07_crc_qc_pp_outs/CRC_offsets.csv")
offsets$sample_ID <- paste0(offsets$slide, "_core", offsets$core)

translate_tx <- function(tx_df, sampleid) {
  local_translation <- tx_df
  sf::st_geometry(local_translation) <- sf::st_geometry(local_translation) + -offsets[offsets$sample_ID == sampleid, c("min_x", "min_y"),] |> as.matrix() |> as.vector()
  return(local_translation)
}
negprobes <- purrr::map2(.x = negprobes, .y = names(negprobes), .f = translate_tx)

# Counting each negative probe in every cell
segs <- sf::read_sf("crc_segmentation_local_coordinates.geojson")
segs <- sf::st_set_crs(segs, NA)
segs <- segs[segs$cell_ID %in% adata$obs_names,]
segs$sample_ID <- plyr::mapvalues(x = segs$cell_ID, from = adata$obs_names, to = as.character(adata$obs$sample_ID))
segs <- split(x = segs, f = segs$sample_ID)

library(sf)
library(ggplot2)
ggplot() + 
  geom_sf(data = segs$CRC2_coreD2) + 
  geom_sf(data = negprobes$CRC2_coreD2, mapping = aes(color = target), size = 0.5) + 
  scale_color_manual(values = BPCells::discrete_palette(name = "stallion")) + 
  theme_void() # Nice

makeNPmeans <- function(.sampleid) {
  
  seg <- segs[[.sampleid]] |> dplyr::select(cell_ID, geometry)
  txdf <- negprobes[[.sampleid]] |> dplyr::select(target, geometry)
  out <- sf::st_join(txdf, seg)
  out$counting <- ifelse(test = is.na(out$cell_ID), yes = 0, no = 1)
  out <- out |> dplyr::select(target, cell_ID, counting) |> sf::st_drop_geometry()
  mf <- tidyr::pivot_wider(data = out, names_from = cell_ID, values_from = counting, values_fill = 0, values_fn = sum, names_sep = "_")
  mf <- dplyr::select(mf, -`NA`)
  mf <- tibble::column_to_rownames(.data = mf, var = "target")
  negmeans <- colMeans(mf)
  cat(.sampleid, "\n") 
  return(negmeans) 
  
}

negmeans_list <- lapply(X = names(segs), FUN = makeNPmeans)
negmeans <- unlist(negmeans_list)
mean(names(negmeans) %in% adata$obs_names) # 1

# Adding to the adata 
adata$obs$negmean <- 0
adata$obs[names(negmeans),]$negmean <- negmeans
anndataR::write_h5ad(object = adata, path = "crc.h5ad", mode = "w")

# Now we can actually run InSituType. We will do it in a supervised manner.
# Need a proper reference profile first
ref <- read.csv(file = url("https://github.com/Nanostring-Biostats/CosMx-Cell-Profiles/raw/refs/heads/main/Human/ColonCRC_6k/ColonCRC_6k.profiles.csv"), row.names = 1)

# Running the insitutype function for each slide separately due to memory constraints
# We want to exclude the erythrocytes
library(Matrix)
cts <- adata$layers$counts |> as("CsparseMatrix") |> magrittr::set_rownames(adata$obs_names) |> magrittr::set_colnames(adata$var_names)
cidx <- adata$obs_names[(adata$obs$celltype_level1 != "Erythrocyte")]

reslist <-  list()
for (slideid in unique(adata$obs$slide)) {
  set.seed(2001)
  idx <- intersect(adata$obs_names[adata$obs$slide == slideid], cidx)
  reslist[[slideid]] <- InSituType::insitutype(x = cts[idx,],
                                                neg = adata$obs[idx,]$negmean,
                                                assay_type = "rna", 
                                                n_clusts=0, 
                                                anchors = NULL, 
                                                reference_profiles = ref,
                                                update_reference_profiles = T, 
                                                rescale = T, 
                                                refit = T
  )
  cat(slideid, "\n")
}
clusts <- sapply(X = reslist, FUN = "[[", "clust")
clusts <- lapply(X = clusts, FUN = as.data.frame) |> 
  dplyr::bind_rows(.id = "slide") |> 
  magrittr::set_colnames(c("slide", "clust"))
celltypes <- data.frame(cell_ID = adata$obs_names, celltype_level1 = adata$obs$celltype_level1)
rownames(celltypes) <- celltypes$cell_ID
celltypes$sup_celltype <- celltypes$celltype_level1
celltypes[rownames(clusts),] <- clusts$clust

# pdf(file = "08_crc_cell-typing_outs/umap_by_celltype_insitutype.pdf", width = 8, height = 5)
BPCells::plot_embedding(
  source =celltypes$sup_celltype,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "CRC: cell types (InSituType)") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, shape = 15)))
# dev.off()

# Adding to the adata
adata$obs$celltype_supervised <- celltypes$sup_celltype
anndataR::write_h5ad(object = adata, path = "crc.h5ad", mode = "w")


## Plotting --------------------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

library(sf)
library(dplyr)

adata <- anndataR::read_h5ad("crc.h5ad", mode = "r+")
seg <- sf::read_sf("crc_segmentation_local_coordinates.geojson")
st_crs(seg) <- NA

seg <- seg[seg$cell_ID %in% adata$obs_names,]
seg$sample_ID <- plyr::mapvalues(x = seg$cell_ID, from = adata$obs_names, to = as.character(adata$obs$sample_ID))
seg$celltype <- plyr::mapvalues(x = seg$cell_ID, from = adata$obs_names, to = adata$obs$celltype_level1)

library(ggplot2)
seg$sample_ID |> unique()
ggplot() + 
  geom_sf(data = seg[seg$sample_ID == "CRC1_coreE2",], mapping = aes(geometry = geometry, fill = celltype), color = "black", size = 0.2) + 
  scale_fill_manual(values = BPCells::discrete_palette("stallion")) + 
  theme_void()
ggplot() + 
  geom_sf(data = seg[seg$sample_ID == "CRC1_coreC4",], mapping = aes(geometry = geometry, fill = celltype), color = "black", size = 0.2) + 
  scale_fill_manual(values = BPCells::discrete_palette("stallion")) + 
  theme_void()
ggplot() + 
  geom_sf(data = seg[seg$sample_ID == "CRC2_coreA4",], mapping = aes(geometry = geometry, fill = celltype), color = "black", size = 0.2) + 
  scale_fill_manual(values = BPCells::discrete_palette("stallion")) + 
  theme_void()
ggplot() + 
  geom_sf(data = seg[seg$sample_ID == "CRC2_coreC3",], mapping = aes(geometry = geometry, fill = celltype), color = "black", size = 0.2) + 
  scale_fill_manual(values = BPCells::discrete_palette("stallion")) + 
  theme_void()
ggplot() + 
  geom_sf(data = seg[seg$sample_ID == "CRC2_coreE4",], mapping = aes(geometry = geometry, fill = celltype), color = "black", size = 0.2) + 
  scale_fill_manual(values = BPCells::discrete_palette("stallion")) + 
  theme_void()
ggplot() + 
  geom_sf(data = seg[seg$sample_ID == "CRC2_coreD3",], mapping = aes(geometry = geometry, fill = celltype), color = "black", size = 0.2) + 
  scale_fill_manual(values = BPCells::discrete_palette("stallion")) + 
  theme_void()


## HieraType -------------------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("crc.h5ad", mode = "r+")

# T and NK cells
tcell_ids <- adata$obs[adata$obs$celltype_level1 == "T.and.NK",] |> rownames()

library(HieraType)
cts <- (adata$layers$counts) |> magrittr::set_colnames(rownames(adata$var)) |> magrittr::set_rownames(rownames(adata$obs)) |> as("CsparseMatrix")
cts <- cts[tcell_ids,]
gr <- (adata$obsp$knn) |> magrittr::set_colnames(rownames(adata$obs)) |> magrittr::set_rownames(rownames(adata$obs))

# Retrieving our pipeline
pipeline_lymphoid <- readRDS("hieratype_custom_lymphoid_pipeline.RDS")

# Finally, we can run HieraType, using this pipeline
lymphoid_typing <- 
  run_pipeline(pipeline = pipeline_lymphoid, 
               counts_matrix = cts, 
               adjacency_matrix = gr[tcell_ids,tcell_ids], 
               celltype_call_threshold = 0.5
  )

# Visualizing
norm <- (adata$layers$lognorm) |> magrittr::set_colnames(rownames(adata$var)) |> magrittr::set_rownames(rownames(adata$obs)) |> as("CsparseMatrix")
norm <- norm[tcell_ids,]
fct <- clusterwise_foldchange_metrics(normed = Matrix::t(norm),
                                      metadata = lymphoid_typing$post_probs$lymphoidmajor,
                                      cluster_column = "celltype_granular")
idxs <- unique(unname(c(unlist(lapply(pipeline_lymphoid$markerslists$tcd4minor, "[[", "index_marker")), 
                        unlist(lapply(pipeline_lymphoid$markerslists$tcd8minor, "[[", "index_marker")), 
                        unlist(lapply(pipeline_lymphoid$markerslists$lymphoidmajor, "[[", "index_marker")))))
nms <- c("NK", names(pipeline_lymphoid$markerslists$tcd4minor), names(pipeline_lymphoid$markerslists$tcd8minor))
hmsubtype <- marker_heatmap(fct, featsuse = c("CD8A", "CD8B", "CD4", "CD27", idxs[is.element(idxs, fct$gene)]), 
                            clusterorder = nms, orient_diagonal = T) + 
  ggplot2::scale_fill_viridis_c(option = "C") + 
  ggplot2::labs(fill = "Mean expression in group", title = "CRC: Cell Sub-type Markers")

# Looks good
print(hmsubtype)
ggplot2::ggsave("08_crc_cell-typing_outs/crc_T_and_NK_subtyping.pdf", height = 4, width = 12)

# Contingencies are reasonable
table(lymphoid_typing$post_probs$lymphoidmajor$celltype_granular)
#   NK  T.CD4.effector    T.CD4.memory     T.CD4.naive T.CD8.cytotoxic T.CD8.exhausted    T.CD8.memory     T.CD8.naive            Treg 
# 8510            1072             513             343            1248            1321             802             603            2061 

# Adding our results to the AnnData object
adata$obs$celltype_granular <- NA
adata$obs[lymphoid_typing$post_probs$lymphoidmajor$cell_ID,]$celltype_granular <- lymphoid_typing$post_probs$lymphoidmajor$celltype_granular

# Saving
saveRDS(lymphoid_typing, file = "08_crc_cell-typing_outs/hieratype_custom_lymphoid_pipeline_results.RDS")
anndataR::write_h5ad(object = adata, path = "crc.h5ad", mode = "w")


## Finalizing ------------------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

# For now, we will rename the "Neutrophil" cluster "Neutrophils.and.Necrosis"
# since this cluster is still somewhat ambiguous.

adata <- anndataR::read_h5ad("crc.h5ad", mode = "r+")
adata$obs$celltype_level1 <- ifelse(test = adata$obs$celltype_level1 == "Neutrophil", 
                                    yes = "Neutrophils.and.Necrosis", 
                                    no = adata$obs$celltype_level1)

# We will also create a "celltype_final" column in the metadata
adata$obs$celltype_final <- dplyr::case_when(
  is.na(adata$obs$celltype_granular) ~ adata$obs$celltype_level1, 
  T ~ adata$obs$celltype_granular
)

# Creating celltype colors for viz downstream
cat(unique(adata$obs$celltype_final), sep = "' = \n'")

celltype_cols <- c(
  'Epithelial' = 'olivedrab',
  'Immune.Myeloid' = 'orange',
  'Endothelial' = 'firebrick',
  'B.cell.Lineage' = 'purple',
  'Fibroblast' = 'dodgerblue',
  'Neutrophils.and.Necrosis' = 'deeppink',
  'T.CD8.memory' = 'thistle',
  'NK' = 'khaki',
  'Myocyte' = 'darkgrey',
  'T.CD8.exhausted' = 'darkmagenta',
  'Treg' = 'darkblue',
  'T.CD4.effector' = 'cyan',
  'T.CD8.naive' = 'plum',
  'T.CD8.cytotoxic' = 'magenta',
  'T.CD4.memory' = 'deepskyblue',
  'T.CD4.naive' = 'lightblue',
  'Erythrocyte' = 'lavenderblush',
  'Neural' = 'black'
)

library(BPCells)
library(ggplot2)
BPCells::plot_embedding(
  source = adata$obs$celltype_final,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  scale_color_manual(values = celltype_cols) +
  labs(x = "UM1", y = "UM2", title = "CRC: final cell types") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, shape = 15)))
# ggsave("08_crc_cell-typing_outs/umap_by_celltype_final.pdf", width = 8, height = 5)

adata$uns$celltype_colors <- celltype_cols

# Validation of cell types with DE
library(data.table)
idx <- !(adata$obs$celltype_final == "Erythrocyte")
psb <- presto::collapse_counts(counts_mat = (Matrix::t(adata$layers$counts[idx,]) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names)),
                               meta_data = data.frame("sample_ID" = adata$obs$sample_ID[idx],
                                                      "celltype" = adata$obs$celltype_final[idx, drop = T]),
                               get_norm = F,
                               varnames = c("celltype", "sample_ID"), 
                               min_cells_per_group = 20
)
psb$meta_data$logUMI <- log(psb$counts_mat |> colSums())

# Modeling
library(presto)
library(lme4)
library(purrr)
library(dplyr)
# presto_res <- presto.presto(
#   y ~ 1 + (1|celltype) + (1|sample_ID) + (1|celltype:sample_ID) + offset(logUMI),
#   psb$meta_data,
#   psb$counts_mat,
#   size_varname = "logUMI",
#   effects_cov = c("celltype"),
#   ncore = 8,
#   min_sigma = 0.05,
#   family = "poisson",
#   nsim = 1000
# )
# saveRDS(object = presto_res, file = "08_crc_cell-typing_outs/presto_model_psb_validation.RDS")
presto_res <- readRDS(file = "08_crc_cell-typing_outs/presto_model_psb_validation.RDS")

# Visulizing results
contrasts_mat <- make_contrast.presto(presto_res, "celltype")
effects_marginal <- contrasts.presto(presto_res, contrasts_mat, one_tailed = T) |> 
  dplyr::mutate(celltype = contrast) |> 
  dplyr::mutate(
    logFC = sign(beta) * log2(exp(abs(beta))), # convert stats to log2
    SD = log2(exp(sigma)),
    zscore = logFC / SD
  ) |> 
  dplyr::select(celltype, feature, logFC, SD, zscore, pvalue) |> 
  dplyr::arrange(pvalue)
effects_marginal$fdr <- p.adjust(effects_marginal$pvalue, method = "BH")

library(ggplot2)
p <- effects_marginal |> 
  ggplot() + 
  scattermore::geom_scattermore(mapping = aes(x = logFC, y = -log10(fdr)), pointsize = 3) + 
  geom_hline(yintercept = -log10(0.01)) + 
  geom_vline(xintercept = 1.5) + 
  ggrepel::geom_text_repel(data = effects_marginal |> filter(fdr < 0.01 & logFC > 1.5), 
                           mapping = aes(x = logFC, y = -log10(fdr), label = feature), 
                           color = "red", size = 3, max.overlaps = 25, box.padding = 0.5) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(hjust = 0.5, size = 10, color = "black"), 
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(size = 12, color = "black"), 
        axis.title.y = element_text(size = 12, color = "black"), 
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5)
  ) + 
  labs(title = "Cell Type Markers", x = "log2FC", y = "-log10(Adj. P-value)") +
  facet_wrap(.~celltype, ncol = 4, scales = "free")
p
ggsave(filename = "08_crc_cell-typing_outs/celltype_markers_by_presto_glmm_volcano_plots.pdf", height = 15, width = 12)

top_marks <- dplyr::group_by(effects_marginal, celltype) |> 
  dplyr::filter(fdr < 0.01 & logFC > 1.5) |>
  dplyr::mutate(rank = order(zscore, decreasing = T)) |> 
  dplyr::group_by(celltype) |> 
  dplyr::top_n(n = -10, wt = rank)
top_marks <- tidyr::pivot_wider(data = top_marks, id_cols = rank, names_from = celltype, values_from = feature) |> 
  dplyr::arrange(rank)
top_marks

BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm[idx,]) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  groups = adata$obs$celltype_final[idx, drop = T], 
                  group_order = colnames(top_marks)[-1],
                  features = top_marks[,-1] |> as.matrix() |> as.vector() |> unique()
) + 
  scale_color_viridis_c(limits = c(-1, 2), oob = scales::squish, option = "B", direction = -1) + 
  ggthemes::theme_hc() +
  theme(axis.text.x = element_text(size = 6)) +
  labs(title = "DE Markers", y = "Cell Type", x = "Gene")
ggsave(filename = "08_crc_cell-typing_outs/celltype_markers_by_presto_glmm_bubble_plot.pdf", height = 5, width = 16)

# Saving
adata$obs$celltype_final <- factor(adata$obs$celltype_final, levels = names(adata$uns$celltype_colors))
anndataR::write_h5ad(object = adata, path = "crc.h5ad", mode = "w")

# Normal colon cells? 
rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("crc.h5ad", mode = "r+")

palette(value = adata$uns$celltype_colors)
tinyplot::plt(y_global ~ x_global | celltype_final, 
              data = adata$obs |> dplyr::filter(sample_ID == "CRC1_coreE2"), 
              col = (adata$uns$celltype_colors),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_global ~ x_global | celltype_final, 
              data = adata$obs |> dplyr::filter(sample_ID == "CRC2_coreD2"), 
              col = (adata$uns$celltype_colors),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)

library(ggplot2)
library(sf)
seg <- sf::read_sf("crc_segmentation_local_coordinates.geojson")
st_crs(seg) <- NA

seg <- seg[seg$cell_ID %in% adata$obs_names,]
seg$sample_ID <- plyr::mapvalues(x = seg$cell_ID, from = adata$obs_names, to = as.character(adata$obs$sample_ID))
seg$cluster <- plyr::mapvalues(x = seg$cell_ID, from = adata$obs_names, to = as.character(adata$uns$iterative_clustering$leiden_res_0.65))

ggplot() + 
  geom_sf(data = seg[seg$sample_ID == "CRC2_coreD2",], mapping = aes(geometry = geometry, fill = cluster), color = "black", size = 0.2) + 
  scale_fill_manual(values = BPCells::discrete_palette("stallion")) + 
  theme_void()
ggplot() + 
  geom_sf(data = seg[seg$sample_ID == "CRC1_coreE2",], mapping = aes(geometry = geometry, fill = cluster), color = "black", size = 0.2) + 
  scale_fill_manual(values = BPCells::discrete_palette("stallion")) + 
  theme_void()

# It seems pretty clear that cluster 13 from resolution 0.65 is normal colon, based on morphology. 
# I will update the data accordingly.

adata$uns$celltype_colors <- c(adata$uns$celltype_colors, "limegreen")
newlevels <- c(gsub(pattern = "Epithelial", replacement = "Epithelial.Malignant", x = levels(adata$obs$celltype_final)), "Epithelial.Normal")

adata$obs$celltype_final <- dplyr::case_when(
  adata$uns$iterative_clustering$leiden_res_0.65 == "13" ~ "Epithelial.Normal", 
  T ~ adata$obs$celltype_final
)
adata$obs$celltype_final <- ifelse(test = adata$obs$celltype_final == "Epithelial", 
                                   yes = "Epithelial.Malignant", 
                                   no = adata$obs$celltype_final)
adata$obs$celltype_final <- factor(x = adata$obs$celltype_final, levels = newlevels)

library(BPCells)
BPCells::plot_embedding(
  source = adata$obs$celltype_final,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  scale_color_manual(values = adata$uns$celltype_colors) +
  labs(x = "UM1", y = "UM2", title = "CRC: final cell types") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, shape = 15)))
# ggsave("08_crc_cell-typing_outs/umap_by_celltype_final.pdf", width = 8, height = 5)

