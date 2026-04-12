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

visium_path <- file.path(args[6], "outs")
obj <- Load10X_Spatial(visium_path)
obj <- obj[, obj$nFeature_Spatial > 0]
########


########
obj@misc$datatype <- "strna"

obj@meta.data$batch <- as.factor(args[7])
Idents(obj) <- obj@meta.data$batch

obj@meta.data$isEmptydrop <- factor("no", levels = c("yes", "no"))
obj@meta.data$isDoublet <- factor("no", levels = c("yes", "no"))
########


#### Basic ####
print("Basic preprocessing")

obj@meta.data$percent.mt <- if (args[9] == "human") {
  PercentageFeatureSet(obj, pattern = "^MT-")
} else if (args[9] == "mouse") {
  PercentageFeatureSet(obj, pattern = "^mt-")
} else {
  rep(0, ncol(obj))
}

obj <- subset(obj, subset = nCount_Spatial < 30 & nFeature_Spatial < 10)

#obj <- SCTransform(obj, vst.flavor = "v2", assay = "Spatial", vars.to.regress = "percent.mt")
obj <- SCTransform(obj, vst.flavor = "v2", assay = "Spatial")
obj@misc$HVGs <- VariableFeatures(obj)

## SVGs ##
obj <- FindSpatiallyVariableFeatures(obj, selection.method = "moransi", features = obj@misc$HVGs)
tmp <- obj@assays$SCT@meta.features[!is.na(obj@assays$SCT@meta.features$moransi.spatially.variable.rank), "moransi.spatially.variable.rank", drop = FALSE]
obj@misc$SVGs <- rownames(tmp[order(tmp$moransi.spatially.variable.rank), , drop = FALSE])[1:2000]
####

obj <- RunPCA(obj, features = obj@misc$SVGs)
obj <- RunUMAP(obj, dims = 1:10)

obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, algorithm = 4, resolution = 0.5)

## banksy clustering ##
obj <- RunBanksy(obj, assay = "SCT", slot = "data", lambda = 0.8, features = "variable", k_geom = 50)

DefaultAssay(obj) <- "BANKSY"

obj <- RunPCA(obj, assay = "BANKSY", reduction.name = "pca.banksy", features = obj@misc$SVGs)
obj <- RunUMAP(obj, dims = 1:10, reduction = "pca.banksy", reduction.name = "umap.banksy")

obj <- FindNeighbors(obj, dims = 1:10, reduction = "pca.banksy")
obj <- FindClusters(obj, algorithm = 4, resolution = 0.5, cluster.name = "banksy_clusters")
####

DefaultAssay(obj) <- "SCT"
obj@meta.data$seurat_clusters <- obj@meta.data$SCT_snn_res.0.5
########

levels(obj@meta.data$seurat_clusters) <- paste0("C", levels(obj@meta.data$seurat_clusters))
levels(obj@meta.data$banksy_clusters) <- paste0("C", levels(obj@meta.data$banksy_clusters))

## spottypes.main (SingleR version) ##
# load("singler_db/hpca.rda")
# ref_lst <- list(hpca)

# obj@meta.data$spottypes.main <- obj@meta.data$banksy_clusters
# ct <- SingleR(test = obj@assays$BANKSY$data,
#               ref = lapply(ref_lst, function(x) x$data),
#               labels = lapply(ref_lst, function(x) x$main_types),
#               clusters = obj@meta.data$banksy_clusters)
# ct$pruned.labels[is.na(ct$pruned.labels)] <- "Unknown"
# levels(obj@meta.data$spottypes.main) <- ct$pruned.labels
####

## spottypes.main (STdeconvolve version) ##
print("spottypes.main (STdeconvolve version)")

pos <- GetTissueCoordinates(obj)
pos <- pos[, 1:2]

DefaultAssay(obj) <- "Spatial"
counts <- GetAssay(obj)$counts
counts <- cleanCounts(counts, min.lib.size = 100)

corpus <- restrictCorpus(counts, removeAbove = 1.0, removeBelow = 0.05)

ldacorpus <- t(as.matrix(corpus))
ldacorpus <- ldacorpus[rowSums(ldacorpus) > 0, ]

ldas <- fitLDA(ldacorpus, Ks = seq(2, 20, by = 1))

optLDA <- optimalModel(models = ldas, opt = "min")

results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta

obj@misc$spottypes.main <- deconProp

topGenes <- topGenes(deconGexp, n = 20)
obj@misc$spottypes.main.topGenes <- topGenes
####

########


#### Process result of pathseq ####
pathseq <- NULL

if (length(grep("\\.rds$", args[8]))) {
  file.copy(args[8], file.path(args[7], "pathseq.rds"), overwrite = TRUE)
  pathseq <- readRDS(file.path(args[7], "pathseq.rds"))
} else if (length(grep("\\.qs$", args[8]))) {
  tmp <- qread(args[8])
  
  tmp <- UpdateSeuratObject(tmp)

  tmp$nCount_SCT <- NULL
  tmp$nFeature_SCT <- NULL
  tmp$SCT_snn_res.0.8 <- NULL
  tmp$seurat_clusters <- NULL
  
  saveRDS(tmp, file.path(args[7], "pathseq.rds"))
  pathseq <- readRDS(file.path(args[7], "pathseq.rds"))
}

if (is.null(pathseq)) stop("pathseq file format error")

## MES & non-MES based on local G statistic ##
print("MES & non-MES based on local G statistic")

coords <- GetTissueCoordinates(obj)
coords <- coords[, c(1,2)]
knn_nb <- knn2nb(knearneigh(coords, k = 6))
w <- nb2listw(knn_nb, style = "W", zero.policy = TRUE)
####

for (microbe in colnames(pathseq@meta.data)[4:ncol(pathseq@meta.data)]) {
  print(microbe)
  ## MES & non-MES ##
  # threshold <- quantile(pathseq[[microbe, drop = TRUE]], 0.75, na.rm = TRUE)
  # enriched_spots <- names(which(pathseq[[microbe, drop = TRUE]] > threshold))
  ####

  ## MES & non-MES based on local G statistic ##
  x <- pathseq@meta.data[rownames(coords), microbe]
  x[is.na(x)] <- 0
  local_g <- localG(x, w, zero.policy = TRUE)
  p_val <- 2 * pnorm(-abs(local_g))
  fdr <- p.adjust(p_val, method = "BH")

  local_g[is.na(local_g)] <- 0
  fdr[is.na(fdr)] <- 1

  enriched_spots <- rownames(coords)[(fdr < 0.05) & (local_g > 0)]
  ####

  obj@misc$MES_nonMES[[microbe]][colnames(obj)] <- "nonMES"
  obj@misc$MES_nonMES[[microbe]][enriched_spots] <- "MES"
  obj@misc$MES_nonMES[[microbe]] <- factor(obj@misc$MES_nonMES[[microbe]], levels = c("MES", "nonMES"))
  Idents(obj) <- obj@misc$MES_nonMES[[microbe]]

  obj@misc$DEGs.MES[[microbe]] <- list(positive = NULL, negative = NULL)

  obj$coord_x_scaled <- scale(coords[, 1])
  obj$coord_y_scaled <- scale(coords[, 2])

  if (length(which(Idents(obj) == "MES")) > 1) {
    # obj@misc$DEGs.MES[[microbe]]$positive <- FindMarkers(obj, ident.1 = "MES", ident.2 = "nonMES", test.use = "LR", latent.vars = c("nCount_Spatial", "coord_x_scaled", "coord_y_scaled"), only.pos = TRUE) %>% filter(p_val_adj < 0.05, avg_log2FC > 1)
    # obj@misc$DEGs.MES[[microbe]]$negative <- FindMarkers(obj, ident.1 = "nonMES", ident.2 = "MES", test.use = "LR", latent.vars = c("nCount_Spatial", "coord_x_scaled", "coord_y_scaled"), only.pos = TRUE) %>% filter(p_val_adj < 0.05, avg_log2FC > 1)
    
    ## Finding DE Markers based on DESpace ##
    spe <- SpatialExperiment(assays = list(counts = GetAssay(obj)$counts),
                             spatialCoords = as.matrix(GetTissueCoordinates(obj)[, 1:2]))
    rowData(spe)$gene_symbol <- rownames(spe)
    colData(spe)$regions <- obj@misc$MES_nonMES[[microbe]]

    spe <- spe[, colSums(counts(spe)) > 0]

    design_temp <- model.matrix(~ colData(spe)$regions)
    keep <- filterByExpr(counts(spe),
                         design = design_temp,
                         min.count = 1,
                         min.total.count = 2,
                         group = colData(spe)$regions)
    spe <- spe[keep, colSums(counts(spe)) > 0]

    spe <- logNormCounts(spe)

    svg_res <- svg_test(spe = spe,
                        cluster_col = "regions",
                        min_counts = 20,
                        min_non_zero_spots = 10,
                        replicates = FALSE,
                        verbose = TRUE)

    individual_res <- individual_svg(spe = spe,
                                     cluster_col = "regions",
                                     min_counts = 20,
                                     min_non_zero_spots = 10,
                                     edgeR_y = svg_res$estimated_y)

    top_regionA <- top_results(cluster_results = individual_res,
                               cluster = 'MES',
                               select = "both",
                               high_low = "high")
    top_regionB <- top_results(cluster_results = individual_res,
                               cluster = 'MES',
                               select = "both",
                               high_low = "low")
    

    obj@misc$DEGs.MES[[microbe]]$positive <- top_regionA$high_gene %>% filter(Cluster_FDR <= 0.05, Cluster_logFC >= 0.6)
    obj@misc$DEGs.MES[[microbe]]$negative <- top_regionB$low_gene %>% filter(Cluster_FDR <= 0.05, Cluster_logFC <= -0.6)
    ####
  }
}

Idents(obj) <- obj@meta.data$batch
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


######## ======== ########


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
