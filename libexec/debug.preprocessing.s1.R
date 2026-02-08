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

obj <- subset(obj, subset = percent.mt < 20)

obj <- SCTransform(obj, vst.flavor = "v2", assay = "Spatial", vars.to.regress = "percent.mt")
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


########
saveRDS(obj, file.path(args[7], "main.rds"))
########