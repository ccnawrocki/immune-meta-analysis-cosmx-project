# Did this in the R console:
# rm(list = ls())
# .rs.restartR(clean = T)

import spatialdata as sd
import anndata as ad
from pathlib import Path
import pandas as pd
import glob

zarr_paths = [p for p in Path(".").rglob("CRC*/*/*.zarr") if p.is_dir()]
zarr_paths = [str(p) for p in zarr_paths]

!mkdir 06_crc_concat_outs

for z in zarr_paths:
  nms = z.split("/")
  slide = nms[1]
  core = nms[2].replace("_outs", "")
  key = slide+"_"+core
  sd.read_zarr(z).tables["table"].write_h5ad("06_crc_concat_outs/"+key+".h5ad")

!ls 06_crc_concat_outs/*.h5ad
adata_paths = glob.glob("06_crc_concat_outs/*.h5ad")

adata = {}
for p in adata_paths:
  nm = p.split("/")[1].replace(".h5ad", "")
  adata[nm] = ad.read_h5ad(p)

adata = ad.concat(adata, index_unique = "_")
adata.shape
# (432005, 6180)

adata.write_h5ad("crc.h5ad")
globals().clear()

