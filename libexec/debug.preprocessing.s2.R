library(qs)
library(Seurat)
library(SeuratWrappers)
library(STdeconvolve)
library(spdep)
library(SpatialExperiment)
library(DESpace)
library(scuttle)
library(edgeR)
library(dplyr)
library(GSVA)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(CellChat)

set.seed(1234)
options(future.globals.maxSize = 1 * 1024^4) # 1TB

########
args <- commandArgs()
if (length(args) != 9) stop("miss parameters")

if (!dir.exists(args[7]))
  dir.create(args[7])

obj <- readRDS(file.path(args[7], "main.rds"))
pathseq <- readRDS(file.path(args[7], "pathseq.rds"))
########


#### GSVA ####
DefaultAssay(obj) <- "Spatial"

obj <- NormalizeData(obj)

db <- if (args[9] == "human") org.Hs.eg.db else if (args[9] == "mouse") org.Mm.eg.db else NULL

## total genes ##
print("GSVA (total genes)")

goannot <- AnnotationDbi::select(db, keys = rownames(obj), keytype = "SYMBOL", columns = "GO")
goannot <- goannot[goannot$ONTOLOGY %in% "BP", ]

genesbygo <- lapply(split(goannot$SYMBOL, goannot$GO), unique)

gsva_obj <- gsva(gsvaParam(GetAssayData(obj)[rownames(obj), ], genesbygo, minSize = 2), verbose = FALSE)

obj@assays$GSVA <- CreateAssayObject(gsva_obj, key = "GSVA_")
####

## hvgs ##
print("GSVA (hvgs)")

goannot <- AnnotationDbi::select(db, keys = obj@misc$HVGs, keytype = "SYMBOL", columns = "GO")
goannot <- goannot[goannot$ONTOLOGY %in% "BP", ]

genesbygo <- lapply(split(goannot$SYMBOL, goannot$GO), unique)

gsva_obj <- gsva(gsvaParam(GetAssayData(obj)[obj@misc$HVGs, ], genesbygo, minSize = 2), verbose = FALSE)

obj@assays$GSVA.hvgs <- CreateAssayObject(gsva_obj, key = "GSVAHVGs_")
####

## svgs ##
print("GSVA (svgs)")

goannot <- AnnotationDbi::select(db, keys = obj@misc$SVGs, keytype = "SYMBOL", columns = "GO")
goannot <- goannot[goannot$ONTOLOGY %in% "BP", ]

genesbygo <- lapply(split(goannot$SYMBOL, goannot$GO), unique)

gsva_obj <- gsva(gsvaParam(GetAssayData(obj)[obj@misc$SVGs, ], genesbygo, minSize = 2), verbose = FALSE)

obj@assays$GSVA.svgs <- CreateAssayObject(gsva_obj, key = "GSVASVGs_")
####

## degs.mes ##
print("GSVA (degs.mes)")

for (microbe in colnames(pathseq@meta.data)[4:ncol(pathseq@meta.data)]) {
  features <- union(rownames(obj@misc$DEGs.MES[[microbe]]$positive), rownames(obj@misc$DEGs.MES[[microbe]]$negative))

  if (length(features) > 2) {
    goannot <- AnnotationDbi::select(db, keys = features, keytype = "SYMBOL", columns = "GO")
    goannot <- goannot[goannot$ONTOLOGY %in% "BP", ]

    genesbygo <- lapply(split(goannot$SYMBOL, goannot$GO), unique)
    genesbygo <- genesbygo[sapply(genesbygo, \(x) length(x) >= 2)]

    if (length(genesbygo) > 0) {
      gsva_par <- gsvaParam(GetAssayData(obj)[features, ], genesbygo, minSize = 2)

      try({
        gsva_obj <- gsva(gsva_par, verbose = FALSE)
        obj@assays[[paste0("GSVA.degs.mes.", microbe)]] <- CreateAssayObject(gsva_obj, key = paste0("GSVADEGsMES", microbe, "_"))
      })
    }
  }
}
####

DefaultAssay(obj) <- "SCT"
########


########
saveRDS(obj, file.path(args[7], "main.rds"))
########