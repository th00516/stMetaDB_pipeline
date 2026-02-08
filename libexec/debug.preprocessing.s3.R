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


#### CellChat v2 ####
obj@meta.data$samples <- obj@meta.data$batch

cellchatdb <- if (args[9] == "human") CellChatDB.human else if (args[9] == "mouse") CellChatDB.mouse else NULL
cellchatdb <- subsetDB(cellchatdb, search = "Secreted Signaling", key = "annotation")

## seurat ##
print("CellChat v2 (seurat clusters)")

cellchat_seurat <- createCellChat(obj@assays$SCT$data, meta = obj@meta.data, group.by = "seurat_clusters", datatype = "spatial",
                                  coordinates = as.matrix(GetTissueCoordinates(obj)[, 1:2]),
                                  spatial.factors = data.frame(ratio = obj@images$slice1@scale.factors$lowres,
                                                               tol = obj@images$slice1@scale.factors$spot / 2))

cellchat_seurat@idents <- droplevels(cellchat_seurat@idents)

cellchat_seurat@DB <- cellchatdb

cellchat_seurat <- subsetData(cellchat_seurat)
cellchat_seurat <- identifyOverExpressedGenes(cellchat_seurat)
cellchat_seurat <- identifyOverExpressedInteractions(cellchat_seurat)

try(cellchat_seurat <- computeCommunProb(cellchat_seurat, contact.range = 100, contact.knn.k = 6, scale.distance = 10))
cellchat_seurat <- filterCommunication(cellchat_seurat, min.cells = 2)

try(cellchat_seurat <- computeCommunProbPathway(cellchat_seurat))

try(cellchat_seurat <- aggregateNet(cellchat_seurat))
####

## banksy ##
print("CellChat v2 (banksy clusters)")

cellchat_banksy <- createCellChat(obj@assays$BANKSY$data, meta = obj@meta.data, group.by = "banksy_clusters", datatype = "spatial",
                                  coordinates = as.matrix(GetTissueCoordinates(obj)[, 1:2]),
                                  spatial.factors = data.frame(ratio = obj@images$slice1@scale.factors$lowres,
                                                               tol = obj@images$slice1@scale.factors$spot / 2))

cellchat_banksy@idents <- droplevels(cellchat_banksy@idents)

cellchat_banksy@DB <- cellchatdb

cellchat_banksy <- subsetData(cellchat_banksy)
cellchat_banksy <- identifyOverExpressedGenes(cellchat_banksy)
cellchat_banksy <- identifyOverExpressedInteractions(cellchat_banksy)

try(cellchat_banksy <- computeCommunProb(cellchat_banksy, contact.range = 100, contact.knn.k = 6, scale.distance = 10))
cellchat_banksy <- filterCommunication(cellchat_banksy, min.cells = 2)

try(cellchat_banksy <- computeCommunProbPathway(cellchat_banksy))

try(cellchat_banksy <- aggregateNet(cellchat_banksy))
####

## spottypes_main ##
# cellchat_spottypes_main <- NULL

# if (length(levels(obj$spottypes.main)) > 1) {
#     cellchat_spottypes_main <- createCellChat(obj@assays$BANKSY$data, meta = obj@meta.data, group.by = "spottypes.main", datatype = "spatial",
#                                               coordinates = as.matrix(GetTissueCoordinates(obj)[, 1:2]),
#                                               spatial.factors = data.frame(ratio = obj@images$slice1@scale.factors$lowres,
#                                                                            tol = obj@images$slice1@scale.factors$spot / 2))

#     cellchat_spottypes_main@idents <- droplevels(cellchat_spottypes_main@idents)

#     cellchat_spottypes_main@DB <- cellchatdb

#     cellchat_spottypes_main <- subsetData(cellchat_spottypes_main)
#     cellchat_spottypes_main <- identifyOverExpressedGenes(cellchat_spottypes_main)
#     cellchat_spottypes_main <- identifyOverExpressedInteractions(cellchat_spottypes_main)

#     cellchat_spottypes_main <- computeCommunProb(cellchat_spottypes_main, contact.range = 100, contact.knn.k = 6, scale.distance = 10)
#     cellchat_spottypes_main <- filterCommunication(cellchat_spottypes_main, min.cells = 2)

#     try(cellchat_spottypes_main <- computeCommunProbPathway(cellchat_spottypes_main))

#     cellchat_spottypes_main <- aggregateNet(cellchat_spottypes_main)
# }
####

## mes vs. nonmes ##
print("CellChat v2 (mes vs. nonmes)")

cellchat_mes_vs_nonmes_list <- list()

for (microbe in colnames(pathseq@meta.data)[4:ncol(pathseq@meta.data)]) {
  print(microbe)

  obj@meta.data$mes_nonmes <- obj@misc$MES_nonMES[[microbe]]

  if (length(unique(obj@meta.data$mes_nonmes)) == 1) next

  if (length(colnames(obj)[obj@meta.data$mes_nonmes == "MES"]) < 2) next
  if (length(colnames(obj)[obj@meta.data$mes_nonmes == "nonMES"]) < 2) next
  
  cellchat_mes_nonmes <- createCellChat(obj@assays$SCT$data, meta = obj@meta.data, group.by = "mes_nonmes", datatype = "spatial",
                                        coordinates = as.matrix(GetTissueCoordinates(obj)[, 1:2]),
                                        spatial.factors = data.frame(ratio = obj@images$slice1@scale.factors$lowres,
                                                                     tol = obj@images$slice1@scale.factors$spot / 2))

  cellchat_mes_nonmes@idents <- droplevels(cellchat_mes_nonmes@idents)

  cellchat_mes_nonmes@DB <- cellchatdb

  cellchat_mes_nonmes <- subsetData(cellchat_mes_nonmes)
  cellchat_mes_nonmes <- identifyOverExpressedGenes(cellchat_mes_nonmes)
  cellchat_mes_nonmes <- identifyOverExpressedInteractions(cellchat_mes_nonmes)

  err <- try(cellchat_mes_nonmes <- computeCommunProb(cellchat_mes_nonmes, contact.range = 100, contact.knn.k = 6, scale.distance = 10))
  if (class(err) == "try-error") next

  cellchat_mes_nonmes <- filterCommunication(cellchat_mes_nonmes, min.cells = 2)

  err <- try(cellchat_mes_nonmes <- computeCommunProbPathway(cellchat_mes_nonmes))
  if (class(err) == "try-error") next

  try(cellchat_mes_nonmes <- aggregateNet(cellchat_mes_nonmes))

  cellchat_mes_vs_nonmes_list[[microbe]]$mes_nonmes <- cellchat_mes_nonmes
}
####

## internal mes or nonmes ##
print("CellChat v2 (internal mes or nonmes)")

cellchat_mes_nonmes_list <- list()

for (microbe in colnames(pathseq@meta.data)[4:ncol(pathseq@meta.data)]) {
  print(microbe)

  obj@meta.data$mes_nonmes <- obj@misc$MES_nonMES[[microbe]]

  if (length(unique(obj@meta.data$mes_nonmes)) == 1) next

  if (length(colnames(obj)[obj@meta.data$mes_nonmes == "MES"]) < 2) next
  if (length(colnames(obj)[obj@meta.data$mes_nonmes == "nonMES"]) < 2) next

  if (length(unique(obj@meta.data["banksy_clusters"][colnames(obj)[obj@meta.data$mes_nonmes == "MES"], ])) < 2) next
  if (length(unique(obj@meta.data["banksy_clusters"][colnames(obj)[obj@meta.data$mes_nonmes == "nonMES"], ])) < 2) next


  obj_mes <- CreateSeuratObject(obj@assays$BANKSY$data[, colnames(obj)[obj@meta.data$mes_nonmes == "MES"]],
                                meta.data = obj@meta.data,
                                assay = "BANKSY")

  obj_nonmes <- CreateSeuratObject(obj@assays$BANKSY$data[, colnames(obj)[obj@meta.data$mes_nonmes == "nonMES"]],
                                   meta.data = obj@meta.data,
                                   assay = "BANKSY")

  cellchat_mes <- createCellChat(obj_mes@assays$BANKSY$counts, meta = obj_mes@meta.data, group.by = "banksy_clusters",
                                 datatype = "spatial",
                                 coordinates = as.matrix(GetTissueCoordinates(obj)[colnames(obj)[obj@meta.data$mes_nonmes == "MES"], 1:2]),
                                 spatial.factors = data.frame(ratio = obj@images$slice1@scale.factors$lowres,
                                                              tol = obj@images$slice1@scale.factors$spot / 2))

  cellchat_nonmes <- createCellChat(obj_nonmes@assays$BANKSY$counts, meta = obj_nonmes@meta.data, group.by = "banksy_clusters",
                                    datatype = "spatial",
                                    coordinates = as.matrix(GetTissueCoordinates(obj)[colnames(obj)[obj@meta.data$mes_nonmes == "nonMES"], 1:2]),
                                    spatial.factors = data.frame(ratio = obj@images$slice1@scale.factors$lowres,
                                                                 tol = obj@images$slice1@scale.factors$spot / 2))


  ## cellchat_mes ##
  cellchat_mes@idents <- droplevels(cellchat_mes@idents)

  cellchat_mes@DB <- cellchatdb

  cellchat_mes <- subsetData(cellchat_mes)
  cellchat_mes <- identifyOverExpressedGenes(cellchat_mes)
  cellchat_mes <- identifyOverExpressedInteractions(cellchat_mes)

  err <- try(cellchat_mes <- computeCommunProb(cellchat_mes, contact.range = 100, contact.knn.k = 6, scale.distance = 10))
  if (class(err) == "try-error") next

  cellchat_mes <- filterCommunication(cellchat_mes, min.cells = 2)

  err <- try(cellchat_mes <- computeCommunProbPathway(cellchat_mes))
  if (class(err) == "try-error") next

  try(cellchat_mes <- aggregateNet(cellchat_mes))
  ####

  ## cellchat_nonmes ##
  cellchat_nonmes@idents <- droplevels(cellchat_nonmes@idents)

  cellchat_nonmes@DB <- cellchatdb

  cellchat_nonmes <- subsetData(cellchat_nonmes)
  cellchat_nonmes <- identifyOverExpressedGenes(cellchat_nonmes)
  cellchat_nonmes <- identifyOverExpressedInteractions(cellchat_nonmes)

  err <- try(cellchat_nonmes <- computeCommunProb(cellchat_nonmes, contact.range = 100, contact.knn.k = 6, scale.distance = 10))
  if (class(err) == "try-error") next

  cellchat_nonmes <- filterCommunication(cellchat_nonmes, min.cells = 2)

  err <- try(cellchat_nonmes <- computeCommunProbPathway(cellchat_nonmes))
  if (class(err) == "try-error") next

  try(cellchat_nonmes <- aggregateNet(cellchat_nonmes))
  ####

  cellchat_mes_nonmes_list[[microbe]]$mes <- cellchat_mes
  cellchat_mes_nonmes_list[[microbe]]$nonmes <- cellchat_nonmes
}
####


cellchat_all <- list(seurat_clusters = cellchat_seurat,
                     banksy_clusters = cellchat_banksy,
                    #  spottypes.main = cellchat_spottypes_main,
                     mes_vs_nonmes = cellchat_mes_vs_nonmes_list,
                     mes_nonmes = cellchat_mes_nonmes_list)

saveRDS(cellchat_all, file.path(args[7], "cellchat.rds"))
########