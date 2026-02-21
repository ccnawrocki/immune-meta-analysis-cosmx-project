rm(list = ls())
.rs.restartR()

library(ggplot2)


### PDAC ### -------------------------------------------------------------------
pdac_samples <- read.csv("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/sample-level-metadata/fov-to-core-master-map.csv", row.names = 1)

# PDAC12 
pdac12 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/PDAC12/flatFiles/PDAC12/PDAC12_metadata_file.csv.gz")
pdac12$core <- plyr::mapvalues(x = pdac12$fov, 
                               from = pdac_samples[pdac_samples$tma == "PDAC12", ]$fov, 
                               to = pdac_samples[pdac_samples$tma == "PDAC12", ]$corenumber)

ggplot() + 
  geom_point(data = pdac12, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = pdac12 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

pdac12 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # | core|      cpc|       fpc|
  # |----:|--------:|---------:|
  # |    7| 308.0598| 107.72597|
  # |   12| 393.7331| 149.08333|
  # |   14| 208.3015|  93.77819| --> KEEP
  # |   21| 409.0807| 161.20354| --> KEEP
  # |   22| 325.7700| 141.61647|
  # |   23| 224.7895| 113.69850| --> KEEP
  # |   29| 166.9665|  79.66223|
  # |   38| 156.7085|  82.10004|
  # |   43| 336.5642| 121.72716| --> KEEP
  # |   44| 354.1252| 129.84997| --> KEEP
  # |   45| 347.4235| 132.06721| --> KEEP
  # |   47| 641.7995| 192.17759|

# PDAC8 
pdac8 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/PDAC8/flatFiles/PDAC8/PDAC8_metadata_file.csv.gz")
pdac8$core <- plyr::mapvalues(x = pdac8$fov, 
                              from = pdac_samples[pdac_samples$tma == "PDAC8", ]$fov, 
                              to = pdac_samples[pdac_samples$tma == "PDAC8", ]$corenumber)

ggplot() + 
  geom_point(data = pdac8, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = pdac8 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

pdac8 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # | core|       cpc|       fpc|
  # |----:|---------:|---------:|
  # |    3| 147.99815|  80.51374| --> KEEP
  # |    4| 180.21484|  73.35612| --> KEEP
  # |    6|  21.57403|  18.15152|
  # |    7|  29.75409|  25.55302|
  # |    9| 199.70228| 103.35627| --> KEEP
  # |   12| 148.98766|  61.31026| --> KEEP
  # |   13|  20.16456|  17.11821|
  # |   14|  16.91587|  14.60716|
  # |   15|  24.87481|  21.45595|
  # |   25|  98.14087|  55.38179|
  # |   37|  19.38392|  16.61608|
  # |   45|  15.79771|  13.73201|

# PDAC10 
pdac10 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/PDAC10/flatFiles/PDAC10/PDAC10_metadata_file.csv.gz")
pdac10$core <- plyr::mapvalues(x = pdac10$fov, 
                              from = pdac_samples[pdac_samples$tma == "PDAC10", ]$fov, 
                              to = pdac_samples[pdac_samples$tma == "PDAC10", ]$corenumber)

ggplot() + 
  geom_point(data = pdac10, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = pdac10 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

pdac10 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # | core|       cpc|      fpc|
  # |----:|---------:|--------:|
  # |    2|  96.34905| 44.46470|
  # |    4| 203.75767| 39.73871| --> KEEP
  # |   11|  90.05576| 48.42518|
  # |   12|  19.54631| 14.06776|
  # |   13|  51.28129| 35.19978|
  # |   19|  28.55530| 17.70880|
  # |   28| 101.48550| 45.35294|
  # |   34|  38.17088| 25.91429|
  # |   41|  22.73451| 18.05605|

# PDAC5 
pdac5 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-pdac-data/runs-0-1-2-3/PDAC5/flatFiles/PDAC5/PDAC5_metadata_file.csv.gz")
pdac5$core <- plyr::mapvalues(x = pdac5$fov, 
                               from = pdac_samples[pdac_samples$tma == "PDAC5", ]$fov, 
                               to = pdac_samples[pdac_samples$tma == "PDAC5", ]$corenumber)

ggplot() + 
  geom_point(data = pdac5, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = pdac5 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

pdac5 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # | core|       cpc|      fpc|
  # |----:|---------:|--------:|
  # |    3|  54.50754| 35.78395|
  # |    5| 103.11752| 54.98754| --> KEEP
  # |    6|  70.75115| 42.59388|
  # |   14|  66.98551| 37.46361|
  # |   15| 100.15334| 53.72850| --> KEEP
  # |   17|  92.52362| 48.96864|
  # |   27| 129.93362| 44.22192| --> KEEP
  # |   28| 111.91158| 40.31640|
  # |   31| 199.15090| 88.09516| --> KEEP
  # |   34|  68.88579| 42.39541|
  # |   36|  75.38112| 40.63565|
  # |   39| 121.41214| 62.85929|


### CRC ### --------------------------------------------------------------------

# ColonPrimaryTMA1
crc_samples <- openxlsx::read.xlsx(xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-crc-data/sample-level-metadata/CRC-to-liver-mets-tma-maps-de-identified.xlsx", 
                                   sheet = 1, cols = 5:6, startRow = 13)
crc1 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-crc-data/ColonPrimaryTMA1/flatFiles/ColonPrimaryTMA1/ColonPrimaryTMA1_metadata_file.csv.gz")
crc1$core <- plyr::mapvalues(x = crc1$fov, 
                              from = crc_samples$fov, 
                              to = crc_samples$corenumber)

ggplot() + 
  geom_point(data = crc1, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = crc1 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

crc1 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # | core|      cpc|      fpc|
  # |----:|--------:|--------:|
  # |    2| 685.3900| 441.9539| --> KEEP
  # |    8| 266.3454| 198.1063| --> KEEP
  # |   11| 419.5133| 294.8138| --> KEEP
  # |   18| 545.2673| 348.7544| --> KEEP
  # |   21| 590.1006| 384.6549| --> KEEP
  # |   23| 266.3015| 198.2730|
  # |   24| 210.6915| 156.7499|

# PrimaryColonTMA2
crc_samples <- openxlsx::read.xlsx(xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-crc-data/sample-level-metadata/CRC-to-liver-mets-tma-maps-de-identified.xlsx", 
                                   sheet = 3, cols = 5:6, startRow = 13)
crc2 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/ting-crc-data/PrimaryColonTMA2/flatFiles/PrimaryColonTMA2/PrimaryColonTMA2_metadata_file.csv.gz")
crc2$core <- plyr::mapvalues(x = crc2$fov, 
                             from = crc_samples$fov, 
                             to = crc_samples$corenumber)

ggplot() + 
  geom_point(data = crc2, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = crc2 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

crc2 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # | core|      cpc|       fpc|
  # |----:|--------:|---------:|
  # |    1| 992.3114| 559.60320| --> KEEP
  # |    2| 213.1678| 164.13439|
  # |    3| 184.1875| 147.95648|
  # |    4| 306.4864| 208.74545| --> KEEP
  # |    5| 411.5674| 285.98393| --> KEEP
  # |    7| 188.5445| 146.77904|
  # |    8| 217.8426| 160.17368|
  # |    9| 329.7836| 230.11543| --> KEEP
  # |   10| 575.8231| 390.37693| --> KEEP
  # |   14| 320.6555| 218.83598|
  # |   15| 109.4235|  74.30801|
  # |   16| 217.9461| 171.26088| --> KEEP
  # |   17| 414.5882| 288.46587| --> KEEP
  # |   21| 165.3263| 130.15954|
  # |   22| 148.5773| 117.91980| --> KEEP


### TNBC ### -------------------------------------------------------------------
tnbc_samples <- read.csv(file = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/sample-level-metadata/1kplex/fov-map.csv", row.names = 1)

# TMA1
tnbc1 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/Sgroi_1K_TMA1/flatFiles/Sgroi_1K_TMA1/Sgroi_1K_TMA1_metadata_file.csv.gz")
tnbc1$core <- plyr::mapvalues(x = tnbc1$fov, 
                              from = tnbc_samples[tnbc_samples$tma == "TMA1",]$fov, 
                              to = tnbc_samples[tnbc_samples$tma == "TMA1",]$corecode)

ggplot() + 
  geom_point(data = tnbc1, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = tnbc1 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

tnbc1 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # |core |       cpc|      fpc|
  # |:----|---------:|--------:|
  # |D1   |  44.65773| 27.17059|
  # |D2   |  21.85506| 16.85641|
  # |D3   |  35.82260| 21.57945|
  # |E1   |  18.52224| 14.68497|
  # |E2   |  75.73214| 46.13247|
  # |E3   |  54.18737| 36.07956|
  # |E4   |  31.00684| 22.01643|
  # |F1   | 119.64694| 56.87395| --> KEEP
  # |F2   |  65.67013| 35.15392|
  # |F3   |  44.50091| 29.96941|
  # |F4   |  43.87484| 29.02487|
  # |G1   | 108.86469| 43.20710| --> KEEP

# TMA2
tnbc2 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/Sgroi_1K_TMA2/flatFiles/Sgroi_1K_TMA2/Sgroi_1K_TMA2_metadata_file.csv.gz")
tnbc2$core <- plyr::mapvalues(x = tnbc2$fov, 
                             from = tnbc_samples[tnbc_samples$tma == "TMA2",]$fov, 
                             to = tnbc_samples[tnbc_samples$tma == "TMA2",]$corecode)

ggplot() + 
  geom_point(data = tnbc2, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = tnbc2 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

tnbc2 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # |core |       cpc|      fpc|
  # |:----|---------:|--------:|
  # |B1   |  67.29205| 35.61523|
  # |C1   |  65.18293| 33.56758|
  # |C2   | 102.06009| 39.33824|
  # |C3   |  91.51773| 50.78986|
  # |C4   | 156.93535| 77.82994| --> KEEP
  # |D1   |  82.08586| 52.37022|
  # |D2   | 170.63139| 62.87144| --> KEEP
  # |D3   |  79.52770| 39.36794|
  # |D4   | 156.91305| 77.01484| --> KEEP
  # |E1   |  50.45067| 32.72775|
  # |E2   |  56.33356| 34.76900|
  # |E3   |  69.21509| 50.31035|
  # |E4   |  61.89629| 30.62508|
  # |F1   |  89.16403| 45.81271|
  # |F2   |  48.75116| 31.14526|
  # |F3   |  71.04636| 52.30203|
  # |F4   |  73.52251| 36.86298|
  # |G1   |  95.46964| 54.62705|

# TMA3
tnbc3 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/Sgroi_1K_TMA3/flatFiles/Sgroi_1K_TMA3/Sgroi_1K_TMA3_metadata_file.csv.gz")
tnbc3$core <- plyr::mapvalues(x = tnbc3$fov, 
                              from = tnbc_samples[tnbc_samples$tma == "TMA3",]$fov, 
                              to = tnbc_samples[tnbc_samples$tma == "TMA3",]$corecode)

ggplot() + 
  geom_point(data = tnbc3, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = tnbc3 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

tnbc3 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # |core |       cpc|      fpc|
  # |:----|---------:|--------:|
  # |D1   |  61.34218| 41.20892|
  # |D2   |  61.89296| 41.47214|
  # |D3   |  56.82469| 41.15425|
  # |D4   | 232.54421| 64.35017| --> KEEP
  # |E1   |  86.11007| 52.63661|
  # |E2   | 117.09174| 65.98039| --> KEEP
  # |E3   |  23.63467| 17.34076|
  # |E4   | 102.30486| 51.53581|
  # |F1   | 101.70469| 57.97099|
  # |F2   | 104.96458| 62.21324|
  # |F3   |  61.28718| 40.37920|
  # |F4   | 100.24163| 44.94885|
  # |G1   |  81.19206| 49.25414|
  
# TMA4
tnbc4 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/sgroi-breast-data/Sgroi_1K_TMA4/flatFiles/Sgroi_1K_TMA4/Sgroi_1K_TMA4_metadata_file.csv.gz")
tnbc4$core <- plyr::mapvalues(x = tnbc4$fov, 
                              from = tnbc_samples[tnbc_samples$tma == "TMA4",]$fov, 
                              to = tnbc_samples[tnbc_samples$tma == "TMA4",]$corecode)

ggplot() + 
  geom_point(data = tnbc4, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = tnbc4 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

tnbc4 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # |core |       cpc|      fpc|
  # |:----|---------:|--------:|
  # |D1   |  80.10550| 47.35968|
  # |D2   |  32.85869| 17.92935|
  # |D3   | 107.31019| 45.21019| --> KEEP
  # |D4   |  41.92239| 30.92594|
  # |E1   |  93.45572| 51.63193|
  # |E2   |  67.10556| 38.32979|
  # |E3   | 111.51171| 54.83567|
  # |E4   | 117.92366| 52.43929| --> KEEP
  # |F1   |  63.37266| 39.91214|
  # |F2   | 137.30448| 34.69293| --> KEEP
  # |F3   |  92.83537| 49.43685|
  # |F4   | 107.49847| 55.95599| --> KEEP
  # |G1   |  57.72037| 27.81086|


### HCC ### --------------------------------------------------------------------
hcc_samples <- lapply(X = c("tma1-pathology", "tma2-pathology"), 
                      FUN = openxlsx::read.xlsx, 
                      xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc-tmas/sample-level-metadata/hcc_tma_maps_de-identified.xlsx", 
                      cols = c(1, 3)) |> 
  dplyr::bind_rows(.id = "slide")
hcc_samples$slide <- paste0("HCC", hcc_samples$slide)
hcc_samples$core <- gsub(pattern = "c", replacement = "", x = hcc_samples$Comments) |> sprintf(fmt = "%02s")

# HCC1
hcc1 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc-tmas/hcc_tma1/flatFiles/hcc_tma1/hcc_tma1_metadata_file.csv.gz")
hcc1$core <- plyr::mapvalues(x = hcc1$fov, 
                             from = hcc_samples[hcc_samples$slide == "HCC1",]$FOV, 
                             to = hcc_samples[hcc_samples$slide == "HCC1",]$core)

ggplot() + 
  geom_point(data = hcc1, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = hcc1 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

hcc1 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # |core |       cpc|       fpc|
  # |:----|---------:|---------:|
  # |01   |  56.42259|  41.15087|
  # |02   |  68.27745|  41.13577|
  # |05   | 101.78689|  46.33999| --> KEEP
  # |06   |  76.88983|  34.68503|
  # |07   |  58.95696|  37.13530|
  # |08   | 210.39499|  83.23478| --> KEEP
  # |09   | 126.52552|  43.31946| --> KEEP
  # |10   | 165.86113|  64.24787| --> KEEP
  # |11   | 155.39164|  56.10789| --> KEEP
  # |12   | 222.20899|  81.56622|
  # |13   |  94.77705|  45.20949|
  # |14   | 458.08250| 122.33379| --> KEEP
  # |16   | 175.10037|  79.10946| --> KEEP
  # |17   | 228.17451|  69.39004| --> KEEP
  # |19   | 103.00358|  46.10041| --> KEEP
  # |20   | 311.51952|  93.61689| --> KEEP
  # |21   | 237.94751|  81.18035| --> KEEP
  # |22   | 393.60488|  97.55024|
  # |23   |  49.26531|  34.61297|

# HCC2
hcc2 <- data.table::fread("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/cosmx/franses-hcc-data/hcc-tmas/hcc_tma2/flatFiles/hcc_tma2/hcc_tma2_metadata_file.csv.gz")
hcc2$core <- plyr::mapvalues(x = hcc2$fov, 
                             from = hcc_samples[hcc_samples$slide == "HCC2",]$FOV, 
                             to = hcc_samples[hcc_samples$slide == "HCC2",]$core)

ggplot() + 
  geom_point(data = hcc2, 
             mapping = aes(x = CenterX_global_px, y = CenterY_global_px, color = as.character(core)), 
             shape = ".") + 
  geom_text(data = hcc2 |> dplyr::group_by(core) |> dplyr::summarise(X = median(CenterX_global_px), Y = median(CenterY_global_px)), 
            mapping = aes(x = X, y = Y, label = core)) + 
  coord_fixed() +
  theme_classic() +
  Seurat::NoLegend()

hcc2 |> dplyr::group_by(core) |> 
  dplyr::summarise(cpc = mean(nCount_RNA), fpc = mean(nFeature_RNA)) |> 
  knitr::kable()
  # |core |        cpc|       fpc|
  # |:----|----------:|---------:|
  # |02   |  641.33628| 153.95924| --> KEEP
  # |03   |  522.12387| 156.23226| --> KEEP
  # |04   |  449.31446| 122.72438| --> KEEP
  # |06   |  484.98348| 117.90112| --> KEEP
  # |08   |  188.63032|  71.78650| --> KEEP
  # |09   |  701.17868| 120.42223| --> KEEP
  # |10   | 1272.49533| 167.77421| --> KEEP
  # |11   |  206.38315|  88.17505| --> KEEP
  # |12   |  118.37527|  55.56088| --> KEEP
  # |14   |  217.20711|  84.51883|
  # |15   |  287.34017|  98.49996| --> KEEP
  # |16   |  296.26660|  93.56439|
  # |17   |  359.61782|  94.25675| --> KEEP
  # |18   |  268.17544|  96.80140| --> KEEP
  # |22   |   68.46324|  38.25011|
  # |23   |   40.42257|  29.03980|

