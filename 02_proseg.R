#!#!# NOTES #!#!#
# - ProSeg version 3.1.0 was used 
# - Device: MacBook Pro (M3 Pro)
# - Did a lot of copy + pasting below, since I wanted to check results as things ran
# - I would recommend doing this on the HPC in the future--especially for the 6K panel!

# - Eventually, I will probably delete/move the ProSeg outputs. Done on 2/10/2026.

rm(list = ls())
.rs.restartR()

# renv::install(packages = c("arrow", "openxlsx", "R.utils"), prompt = F)

library(dplyr)
dir.create("02_proseg_outs")

### PDAC ### -------------------------------------------------------------------
pdac_samples <- read.csv("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/sample-level-metadata/fov-to-core-master-map.csv", row.names = 1)

## PDAC12 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/PDAC12/flatFiles/PDAC12/PDAC12_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = pdac_samples[pdac_samples$tma == "PDAC12", ]$fov, 
                           to = pdac_samples[pdac_samples$tma == "PDAC12", ]$corenumber)
tx <- filter(tx, core %in% c(14, 21, 23, 43, 44, 45))

tx <- split(tx, tx$core)
names(tx) <- sprintf("%02s", names(tx))
dir.create("02_proseg_outs/PDAC12")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/PDAC12/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/PDAC12
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

# Looks good
test <- arrow::read_parquet("02_proseg_outs/PDAC12/core45_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/PDAC12", pattern = "\\.csv", full.names = T))

## PDAC8 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/PDAC8/flatFiles/PDAC8/PDAC8_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = pdac_samples[pdac_samples$tma == "PDAC8", ]$fov, 
                           to = pdac_samples[pdac_samples$tma == "PDAC8", ]$corenumber)
tx <- filter(tx, core %in% c(3, 4, 9, 12))

# Splitting by FOV
tx <- split(tx, tx$core)
names(tx) <- sprintf("%02s", names(tx))
dir.create("02_proseg_outs/PDAC8")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/PDAC8/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/PDAC8
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/PDAC8/core12_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/PDAC8", pattern = "\\.csv", full.names = T))

## PDAC10 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/PDAC10/flatFiles/PDAC10/PDAC10_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = pdac_samples[pdac_samples$tma == "PDAC10", ]$fov, 
                           to = pdac_samples[pdac_samples$tma == "PDAC10", ]$corenumber)
tx <- filter(tx, core %in% c(4))

# Splitting by FOV
tx <- split(tx, tx$core)
names(tx) <- sprintf("%02s", names(tx))
dir.create("02_proseg_outs/PDAC10")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/PDAC10/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/PDAC10
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/PDAC10/core04_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/PDAC10", pattern = "\\.csv", full.names = T))

## PDAC5 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/PDAC5/flatFiles/PDAC5/PDAC5_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = pdac_samples[pdac_samples$tma == "PDAC5", ]$fov, 
                           to = pdac_samples[pdac_samples$tma == "PDAC5", ]$corenumber)
tx <- filter(tx, core %in% c(5, 15, 27, 31))

# Splitting by FOV
tx <- split(tx, tx$core)
names(tx) <- sprintf("%02s", names(tx))
dir.create("02_proseg_outs/PDAC5")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/PDAC5/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/PDAC5
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/PDAC5/core31_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/PDAC5", pattern = "\\.csv", full.names = T))


### TNBC ### -------------------------------------------------------------------
tnbc_samples <- read.csv(file = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/sample-level-metadata/1kplex/fov-map.csv", row.names = 1)

## TMA 1 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/Sgroi_1K_TMA1/flatFiles/Sgroi_1K_TMA1/Sgroi_1K_TMA1_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = tnbc_samples[tnbc_samples$tma == "TMA1",]$fov, 
                           to = tnbc_samples[tnbc_samples$tma == "TMA1",]$corecode)
tx <- filter(tx, core %in% c("F1", "G1"))

tx <- split(tx, tx$core)
dir.create("02_proseg_outs/TNBC1")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/TNBC1/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/TNBC1
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/TNBC1/coreG1_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/TNBC1", pattern = "\\.csv", full.names = T))

## TMA 2 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/Sgroi_1K_TMA2/flatFiles/Sgroi_1K_TMA2/Sgroi_1K_TMA2_tx_file.csv")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = tnbc_samples[tnbc_samples$tma == "TMA2",]$fov, 
                           to = tnbc_samples[tnbc_samples$tma == "TMA2",]$corecode)
tx <- filter(tx, core %in% c("C4", "D2", "D4"))

tx <- split(tx, tx$core)
dir.create("02_proseg_outs/TNBC2")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/TNBC2/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/TNBC2
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/TNBC2/coreD4_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/TNBC2", pattern = "\\.csv", full.names = T))

## TMA3 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/Sgroi_1K_TMA3/flatFiles/Sgroi_1K_TMA3/Sgroi_1K_TMA3_tx_file.csv")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = tnbc_samples[tnbc_samples$tma == "TMA3",]$fov, 
                           to = tnbc_samples[tnbc_samples$tma == "TMA3",]$corecode)
tx <- filter(tx, core %in% c("D4", "E2"))

tx <- split(tx, tx$core)
dir.create("02_proseg_outs/TNBC3")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/TNBC3/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/TNBC3
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/TNBC3/coreE2_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/TNBC3", pattern = "\\.csv", full.names = T))

## TMA4 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/Sgroi_1K_TMA4/flatFiles/Sgroi_1K_TMA4/Sgroi_1K_TMA4_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = tnbc_samples[tnbc_samples$tma == "TMA4",]$fov, 
                           to = tnbc_samples[tnbc_samples$tma == "TMA4",]$corecode)
tx <- filter(tx, core %in% c("D3", "E4", "F2", "F4"))

tx <- split(tx, tx$core)
dir.create("02_proseg_outs/TNBC4")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/TNBC4/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/TNBC4
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/TNBC4/coreF4_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/TNBC4", pattern = "\\.csv", full.names = T))


### CRC ### --------------------------------------------------------------------
crc_samples <- openxlsx::read.xlsx(xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-crc-data/sample_map.xlsx", 
                                   sheet = "fov map")

## CRC1 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-crc-data/ColonPrimaryTMA1/flatFiles/ColonPrimaryTMA1/ColonPrimaryTMA1_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = crc_samples[crc_samples$tma == "prim1",]$fov, 
                           to = crc_samples[crc_samples$tma == "prim1",]$core)
tx <- filter(tx, core %in% c("B1", "B2", "E2", "F3", "C4"))

tx <- split(tx, tx$core)
dir.create("02_proseg_outs/CRC1")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/CRC1/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/CRC1
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/CRC1/coreB1_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/CRC1", pattern = "\\.csv", full.names = T))

## CRC2 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-crc-data/PrimaryColonTMA2/flatFiles/PrimaryColonTMA2/PrimaryColonTMA2_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = crc_samples[crc_samples$tma == "prim2",]$fov, 
                           to = crc_samples[crc_samples$tma == "prim2",]$core)
tx <- filter(tx, core %in% c("A4", "D4", "E4", "C3", "D3", "D2", "E2", "D1"))

tx <- split(tx, tx$core)
dir.create("02_proseg_outs/CRC2")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/CRC2/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/CRC2
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/CRC2/coreA4_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/CRC2", pattern = "\\.csv", full.names = T))


### HCC ### --------------------------------------------------------------------
hcc_samples <- lapply(X = c("tma1-pathology", "tma2-pathology"), 
                      FUN = openxlsx::read.xlsx, 
                      xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc-tmas/sample-level-metadata/hcc_tma_maps_de-identified.xlsx", 
                      cols = c(1, 3)) |> 
  dplyr::bind_rows(.id = "slide")
hcc_samples$slide <- paste0("HCC", hcc_samples$slide)
hcc_samples$core <- gsub(pattern = "c", replacement = "", x = hcc_samples$Comments) |> sprintf(fmt = "%02s")

## HCC1 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc-tmas/hcc_tma1/flatFiles/hcc_tma1/hcc_tma1_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = hcc_samples[hcc_samples$slide == "HCC1",]$FOV, 
                           to = hcc_samples[hcc_samples$slide == "HCC1",]$core)
tx <- filter(tx, core %in% c("05", "08", "09", "10", "11", "14", "16", "17", "19", "20", "21"))

tx <- split(tx, tx$core)
dir.create("02_proseg_outs/HCC1")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/HCC1/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/HCC1
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/HCC1/core11_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/HCC1", pattern = "\\.csv", full.names = T))

## HCC2 ##
tx <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc-tmas/hcc_tma2/flatFiles/hcc_tma2/hcc_tma2_tx_file.csv.gz")
head(tx)

tx$core <- plyr::mapvalues(x = tx$fov, 
                           from = hcc_samples[hcc_samples$slide == "HCC2",]$FOV, 
                           to = hcc_samples[hcc_samples$slide == "HCC2",]$core)
tx <- filter(tx, core %in% c("02", "03", "04", "06", "08", "09", "10", "11", "12", "15", "17", "18"))

tx <- split(tx, tx$core)
dir.create("02_proseg_outs/HCC2")
purrr::map2(.x = tx, .y = paste("02_proseg_outs/HCC2/core", names(tx), ".csv", sep = ""), 
            .f = data.table::fwrite)

# $ cd /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/analysis/immune-meta-analysis-project/02_proseg_outs/HCC2
# $ dirlist=($(ls ${F}*.csv))
# $ for f in $dirlist; do print ${f//.csv/}; done
# $ for f in $dirlist; do mkdir ${f//.csv/_outs}; done
# $ for f in $dirlist; do proseg $f --output-path ${f//.csv/_outs} --gene-column target --x-column x_global_px --y-column y_global_px --z-column z --compartment-column CellComp --compartment-nuclear Nuclear --fov-column fov --cell-id-column cell --cell-id-unassigned '' --cell-assignment-column cell_ID --cell-assignment-unassigned '0' --excluded-genes '^(SystemControl|Negative)' --coordinate-scale 0.12028; done

test <- arrow::read_parquet("02_proseg_outs/HCC2/core06_outs/proseg-output.zarr/shapes/cell_boundaries/shapes.parquet")
test <- sf::st_as_sf(test)
plot(test, asp = 1, col = "grey")

# Do not need these anymore
file.remove(list.files(path = "02_proseg_outs/HCC2", pattern = "\\.csv", full.names = T))

