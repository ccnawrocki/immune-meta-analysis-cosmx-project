rm(list = ls())
.rs.restartR(clean = T)

library(Matrix)
library(BPCells)
library(ggplot2)

ad <- reticulate::import("anndata")
adata <- ad$read_h5ad("tnbc.h5ad")

## Plotting some canonical markers ---------------------------------------------
# dir.create(path = "11_tnbc_cell-typing_outs")
# pdf(file = "11_tnbc_cell-typing_outs/umap_by_canonical_markers.pdf", width = 12, height = 10)
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
               "IGHG1",
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
  features = c("VWF", 
               "ERBB2", 
               "MYL9", 
               "ELANE", 
               "KRT19", 
               "FPR1"
  ), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71), 
  quantile_range = c(0.01, 0.99)
) & 
  labs(x = "UM1", y = "UM2")
# dev.off()


## Plotting some other variables -----------------------------------------------
# pdf(file = "11_tnbc_cell-typing_outs/umap_by_other_variables.pdf", width = 6, height = 5)
plot_embedding(
  source = adata$obs$slide,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "TNBC: slide")
plot_embedding(
  source = adata$obs$patient,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "TNBC: patient")
plot_embedding(
  source = adata$obs$sample_ID,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "TNBC: sample ID") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, shape = 15)))
# dev.off()


## Clusters that I think look good ---------------------------------------------
pdf(file = "11_tnbc_cell-typing_outs/umap_by_leiden_clusters.pdf", width = 6, height = 5)
plot_embedding(
  source = adata$uns$iterative_clustering$leiden_res_0.9,
  embedding = adata$obsm["X_umap"],
  rasterize = T, labels_discrete = T
) + 
  labs(x = "UM1", y = "UM2", title = "TNBC: leiden clustering (res = 0.9)") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, shape = 15)))
dev.off()


## Autofluorescence score ------------------------------------------------------
adata <- anndataR::read_h5ad("tnbc.h5ad", mode = "r+")
probes <- readRDS(file = url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/refs/heads/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))
probes <- probes$Hs_UCC
probes$gene[probes$barcode |> grepl(pattern = ".*BB.*BB.*BB.*") |> which()]
bluetargets <- probes$gene[probes$barcode |> grepl(pattern = ".*BB.*BB.*BB.*") |> which()] |> intersect(adata$var |> rownames())

afscore <- Seurat::AddModuleScore(Seurat::CreateAssayObject(data = t(adata$layers$lognorm) |> magrittr::set_rownames(rownames(adata$var)) |> magrittr::set_colnames(rownames(adata$obs)) |> as("CsparseMatrix") |> Seurat::ScaleData()), 
                                  features = list(bluetargets), nbin = 7, name = "AF_score")

afscore_emb <- adata$obsm["X_umap"] |> as.data.frame() |> magrittr::set_rownames(adata$obs |> rownames()) |> magrittr::set_colnames(c("UM1", "UM2"))
afscore_emb$AF_score <- plyr::mapvalues(x = rownames(afscore_emb), from = rownames(afscore), to = afscore$AF_score1) |> as.numeric()
ggplot() + 
  scattermore::geom_scattermore(data = afscore_emb, mapping = aes(x = UM1, y = UM2, colour = AF_score)) + 
  scale_color_viridis_c(option = "G") + 
  theme_classic()
# ggsave(filename = "11_tnbc_cell-typing_outs/umap_by_AF_score.pdf", width = 6, height = 5)

# Don't really see any RBCs.


## Finding cluster markers with presto GLMM ------------------------------------
# Pseudo-bulking
library(data.table)
idx <- !(adata$uns$iterative_clustering$leiden_res_0.9 %in% c("2", "3", "6", "10", "9"))
psb <- presto::collapse_counts(counts_mat = (Matrix::t(adata$layers$counts[idx,]) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names)),
                               meta_data = data.frame("sample_ID" = adata$obs$sample_ID[idx],
                                                      "cluster" = adata$uns$iterative_clustering$leiden_res_0.9[idx, drop = T]),
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
saveRDS(object = presto_res, file = "11_tnbc_cell-typing_outs/presto_model_psb.RDS")
presto_res <- readRDS(file = "11_tnbc_cell-typing_outs/presto_model_psb.RDS")

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
# ggsave(filename = "11_tnbc_cell-typing_outs/cluster_markers_by_presto_glmm_volcano_plots.pdf", height = 10, width = 12)

top_marks <- dplyr::group_by(effects_marginal, cluster) |> 
  dplyr::filter(fdr < 0.01 & logFC > 1.5) |>
  dplyr::mutate(rank = order(zscore, decreasing = T)) |> 
  dplyr::group_by(cluster) |> 
  dplyr::top_n(n = -10, wt = rank)
top_marks <- tidyr::pivot_wider(data = top_marks, id_cols = rank, names_from = cluster, values_from = feature) |> 
  dplyr::arrange(rank)
top_marks

library(BPCells)
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm[idx,]) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  groups = adata$uns$iterative_clustering$leiden_res_0.9[idx, drop = T], 
                  group_order = colnames(top_marks)[-1],
                  features = top_marks[,-1] |> as.matrix() |> as.vector() |> unique()
) + 
  scale_color_viridis_c(limits = c(-1, 2), oob = scales::squish, option = "B", direction = -1) + 
  ggthemes::theme_hc() +
  theme(axis.text.x = element_text(size = 8)) +
  labs(title = "DE Markers", y = "Cluster", x = "Gene")
# ggsave(filename = "11_tnbc_cell-typing_outs/cluster_markers_by_presto_glmm_bubble_plot.pdf", height = 5, width = 16)

# Looking at some cores: 
d <- adata$obs
d$cluster <- adata$uns$iterative_clustering$leiden_res_0.9

# pdf(file = "11_tnbc_cell-typing_outs/clusters_in_space.pdf", height = 10, width = 10)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "TNBC1_G1"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "TNBC2_D4"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "TNBC3_D4"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_global ~ x_global | cluster, 
              data = d |> dplyr::filter(sample_ID == "TNBC4_F2"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
# dev.off()


## Annotating ------------------------------------------------------------------
adata$obs$celltype_level1 <- dplyr::case_when(
  adata$uns$iterative_clustering$leiden_res_0.9 == 1 ~ "Endothelial",
  adata$uns$iterative_clustering$leiden_res_0.9 == 2 ~ "Epithelial.Malignant", 
  adata$uns$iterative_clustering$leiden_res_0.9 == 3 ~ "Epithelial.Malignant",
  adata$uns$iterative_clustering$leiden_res_0.9 == 4 ~ "B.cell",
  adata$uns$iterative_clustering$leiden_res_0.9 == 5 ~ "Immune.Myeloid",
  adata$uns$iterative_clustering$leiden_res_0.9 == 6 ~ "Epithelial.Malignant",
  adata$uns$iterative_clustering$leiden_res_0.9 == 7 ~ "Neutrophil",
  adata$uns$iterative_clustering$leiden_res_0.9 == 8 ~ "Fibroblast",
  adata$uns$iterative_clustering$leiden_res_0.9 == 9 ~ "Epithelial.Malignant",
  adata$uns$iterative_clustering$leiden_res_0.9 == 10 ~ "Epithelial.Malignant",
  adata$uns$iterative_clustering$leiden_res_0.9 == 11 ~ "T.and.NK", 
  adata$uns$iterative_clustering$leiden_res_0.9 == 12 ~ "Plasma"
)

# pdf(file = "11_tnbc_cell-typing_outs/umap_by_celltype_level1.pdf", width = 6, height = 5)
plot_embedding(
  source = adata$obs$celltype_level1,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "TNBC: cell types (level 1)") + 
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3, shape = 15)))
# dev.off()

# Saving
anndataR::write_h5ad(object = adata, path = "tnbc.h5ad", mode = "w")


## InSituType ------------------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("tnbc.h5ad", mode = "r+")

# Tricky part: need the mean negative probe expression in each cell
library(DBI)
library(duckdb)
library(dplyr)

tnbc_stow <- "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data"
txpaths <- list.files(pattern = "tx_file.csv.gz", 
                      recursive = T, 
                      path = tnbc_stow)
txpaths <- txpaths[grepl(pattern = "Sgroi_1K", x = txpaths)]

# Getting all the negative probe locations for our samples
negprobe_list <- list()
for (p in txpaths) {
  con <- dbConnect(duckdb(), dbdir = ":memory:")
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  P <- file.path(tnbc_stow, p)
  tx <- dbGetQuery(
    con, 
    sprintf("
            FROM read_csv('%s')
            WHERE target IN ('NegativeAdd', 'Negative8', 'Negative6', 'Negative2', 'Negative7', 'Negative1', 'Negative9', 'Negative10', 'Negative5', 'Negative4', 'Negative3')
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
tnbc_samples <- read.csv(file = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/sample-level-metadata/1kplex/fov-map.csv", row.names = 1)
tnbc_samples$slide <- paste0("Sgroi_1K_", tnbc_samples$tma)
negprobes$core <- plyr::mapvalues(x = paste(negprobes$slide, negprobes$fov, sep = "_"), from = paste(tnbc_samples$slide, tnbc_samples$fov, sep = "_"), to = tnbc_samples$corecode)
negprobes$sample_ID <- paste0("TNBC", gsub(pattern = "Sgroi_1K_TMA", replacement = "", x = negprobes$slide), "_", negprobes$core)
negprobes <- negprobes[negprobes$sample_ID %in% unique(adata$obs$sample_ID),]

# Converting to microns
negprobes$x_global_um <- negprobes$x_global_px*0.12028
negprobes$y_global_um <- negprobes$y_global_px*0.12028

# Creating separate sf objects for each sample
negprobes <- split(x = negprobes, f = negprobes$sample_ID)
negprobes <- lapply(X = negprobes, FUN = sf::st_as_sf, coords = c("x_global_um", "y_global_um"))
negprobes[["TNBC1_F1"]] |> plot() # Nice

# Transforming the coordinates to standard scale, as we did with the segmentations.
offsets <- read.csv(file = "10_tnbc_qc_pp_outs/TNBC_offsets.csv")
offsets$sample_ID <- paste(offsets$slide, gsub(pattern = "core", replacement = "", x = offsets$corecode), sep = "_")

translate_tx <- function(tx_df, sampleid) {
  local_translation <- tx_df
  sf::st_geometry(local_translation) <- sf::st_geometry(local_translation) + -offsets[offsets$sample_ID == sampleid, c("min_x", "min_y"),] |> as.matrix() |> as.vector()
  return(local_translation)
}
negprobes <- purrr::map2(.x = negprobes, .y = names(negprobes), .f = translate_tx)

# Counting each negative probe in every cell
segs <- sf::read_sf("tnbc_segmentation_local_coordinates.geojson")
segs <- sf::st_set_crs(segs, NA)
segs <- segs[segs$cell_ID %in% adata$obs_names,]
segs$sample_ID <- stringr::str_split(string = segs$cell_ID, pattern = "_", n = 2, simplify = T)[,2] |> gsub(pattern = "core", replacement = "", x = _)
segs <- split(x = segs, f = segs$sample_ID)

library(sf)
library(ggplot2)
ggplot() + 
  geom_sf(data = segs$TNBC1_G1) + 
  geom_sf(data = negprobes$TNBC1_G1, mapping = aes(color = target), size = 0.5) + 
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
anndataR::write_h5ad(object = adata, path = "tnbc.h5ad", mode = "w")

# Now we can actually run InSituType. We will do it in a supervised manner.
# Need a proper reference profile first
ref <- read.csv(file = url("https://github.com/Nanostring-Biostats/CosMx-Cell-Profiles/raw/refs/heads/main/Human/IO_1k/IO_1k.profiles.csv"), row.names = 1)

library(Matrix)
cts <- adata$layers$counts |> as("CsparseMatrix") |> magrittr::set_rownames(adata$obs_names) |> magrittr::set_colnames(adata$var_names)
cidx <- adata$obs_names[adata$obs$celltype_level1 != "Epithelial.Malignant"]

res <- InSituType::insitutype(x = cts[cidx,],
                              neg = adata$obs[cidx,]$negmean,
                              assay_type = "rna", 
                              n_clusts = 0, 
                              anchors = NULL, 
                              reference_profiles = ref,
                              update_reference_profiles = T, 
                              rescale = T, 
                              refit = T
)

celltypes <- data.frame(cell_ID = rownames(adata$obs), prior_clust = adata$obs$celltype_level1 |> as.character())
rownames(celltypes) <- celltypes$cell_ID
celltypes$sup_celltype <- celltypes$prior_clust
celltypes[names(res$clust),]$sup_celltype <- res$clust

# pdf(file = "11_tnbc_cell-typing_outs/umap_by_celltype_insitutype.pdf", width = 8, height = 5)
BPCells::plot_embedding(
  source = celltypes$sup_celltype,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  labs(x = "UM1", y = "UM2", title = "TNBC: cell types (InSituType)") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, shape = 15)))
# dev.off()

# Adding to the adata
adata$obs$celltype_supervised <- celltypes$sup_celltype
anndataR::write_h5ad(object = adata, path = "tnbc.h5ad", mode = "w")


## HieraType -------------------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("tnbc.h5ad", mode = "r+")

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
  ggplot2::labs(fill = "Mean expression in group", title = "TNBC: Cell Sub-type Markers")

# Looks good
print(hmsubtype)
ggplot2::ggsave("11_tnbc_cell-typing_outs/hcc_T_and_NK_subtyping.pdf", height = 4, width = 12)

# Contingencies are reasonable
table(lymphoid_typing$post_probs$lymphoidmajor$celltype_granular)
#   NK  T.CD4.effector    T.CD4.memory     T.CD4.naive T.CD8.cytotoxic T.CD8.exhausted    T.CD8.memory     T.CD8.naive            Treg 
# 4445             476             438             337             538             703             629             300             823 

# Adding our results to the AnnData object
adata$obs$celltype_granular <- NA
adata$obs[lymphoid_typing$post_probs$lymphoidmajor$cell_ID,]$celltype_granular <- lymphoid_typing$post_probs$lymphoidmajor$celltype_granular

# Saving
saveRDS(lymphoid_typing, file = "11_tnbc_cell-typing_outs/hieratype_custom_lymphoid_pipeline_results.RDS")
anndataR::write_h5ad(object = adata, path = "tnbc.h5ad", mode = "w")


## Finalizing ------------------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

# For now, we will rename the "Neutrophil" cluster "Neutrophils.and.Necrosis"
# since this cluster is still somewhat ambiguous.

adata <- anndataR::read_h5ad("tnbc.h5ad", mode = "r+")
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
  'Endothelial' = 'firebrick',
  'Epithelial.Malignant' = 'olivedrab',
  'B.cell' = 'purple',
  'Neutrophils.and.Necrosis' = 'deeppink',
  'Fibroblast' = 'dodgerblue',
  'Immune.Myeloid' = 'orange',
  'T.CD8.cytotoxic' = 'magenta',
  'T.CD4.memory' = 'deepskyblue',
  'Plasma' = 'red',
  'Treg' = 'darkblue',
  'NK' = 'khaki',
  'T.CD8.memory' = 'thistle',
  'T.CD8.exhausted' = 'darkmagenta',
  'T.CD8.naive' = 'plum',
  'T.CD4.effector' = 'cyan',
  'T.CD4.naive' = 'lightblue'
)

library(BPCells)
library(ggplot2)
BPCells::plot_embedding(
  source = adata$obs$celltype_final,
  embedding = adata$obsm$X_umap,
  rasterize = T, labels_discrete = F
) + 
  scale_color_manual(values = celltype_cols) +
  labs(x = "UM1", y = "UM2", title = "TNBC: final cell types") + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, shape = 15)))
# ggsave("11_tnbc_cell-typing_outs/umap_by_celltype_final.pdf", width = 8, height = 5)

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
# saveRDS(object = presto_res, file = "11_tnbc_cell-typing_outs/presto_model_psb_validation.RDS")
presto_res <- readRDS(file = "11_tnbc_cell-typing_outs/presto_model_psb_validation.RDS")

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
ggsave(filename = "11_tnbc_cell-typing_outs/celltype_markers_by_presto_glmm_volcano_plots.pdf", height = 15, width = 12)

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
ggsave(filename = "11_tnbc_cell-typing_outs/celltype_markers_by_presto_glmm_bubble_plot.pdf", height = 5, width = 16)

# Final save
adata$obs$celltype_final <- factor(adata$obs$celltype_final, levels = names(adata$uns$celltype_colors))
anndataR::write_h5ad(object = adata, path = "tnbc.h5ad", mode = "w")

