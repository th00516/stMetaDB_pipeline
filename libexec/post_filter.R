library(Seurat)
library(SCRuB)
library(phyloseq)
library(spdep)


args <- commandArgs()

sample_id <- args[6]
samples_folder <- "/home/dataset-assist-0/PUB/yuhao/stmetadb_tmp/PMO/yuhao/DevHub/stmeta/workDir_new2/pipeline/2.gatk_pathseq/outputs/"
output_path <- "outputs/"

dir.create(output_path, recursive = TRUE)


obj <- readRDS(paste0(samples_folder, sample_id, ".raw_matrix.rds"))
microbe_counts <- obj@meta.data[4:(ncol(obj@meta.data) - 1), drop = FALSE]
microbe_counts[is.na(microbe_counts)] <- 0


sample_metadata <- obj@meta.data[1:3]

microbe_counts_total <- rowSums(microbe_counts)


#### 1. Filter contaminated microbe spots by SCRuB ####
scrub_counts <- microbe_counts[names(microbe_counts_total)[microbe_counts_total > 0], , drop = FALSE]
scrub_counts <- as.matrix(scrub_counts)

minTotalUMICounts <- 10
blank_spots <- names(microbe_counts_total)[
  (microbe_counts_total < minTotalUMICounts) & (microbe_counts_total > 0)
]
scrub_metadata <- data.frame(
  is_control = rep(FALSE, nrow(scrub_counts)),
  sample_type = rep("microbe spot", nrow(scrub_counts)),
  row.names = rownames(scrub_counts)
)
scrub_metadata[blank_spots, "is_control"] <- TRUE
scrub_metadata[blank_spots, "sample_type"] <- "blank spot"

err <- try(scr_out <- SCRuB(scrub_counts, metadata = scrub_metadata))

if (class(err) == "try-error") {
  decontaminated_microbe_counts <- microbe_counts
} else {
  decontaminated_microbe_counts <- scr_out$decontaminated_samples
}

#### 2. Filter microbe taxonomies with low abundance ####
ps <- phyloseq(otu_table(decontaminated_microbe_counts, taxa_are_rows = FALSE),
                         sample_data(sample_metadata[rownames(decontaminated_microbe_counts),]))

abundanceThreshold <- 0.0001
# prevalenceThreshold <- 0.05
# minUMICounts <- 5

# ps <- prune_taxa(taxa_sums(ps) >= minUMICounts, ps)
# if (sum(ps@otu_table > 0) == 0)
#   quit(status = 0)

ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps.rel@otu_table[is.na(ps.rel@otu_table)] <- 0

ps.filtered <- filter_taxa(ps.rel, function(x) sum(x) > abundanceThreshold, prune = TRUE)
if (sum(ps.filtered@otu_table > 0) == 0)
  quit(status = 0)

# prev_threshold <- round(prevalenceThreshold * nsamples(ps.filtered))
# ps.filtered <- filter_taxa(ps.filtered, function(x) sum(x > 0) > prev_threshold, prune = TRUE)
# if (sum(ps.filtered@otu_table > 0) == 0)
#   quit(status = 0)

cleaned_microbes <- taxa_names(ps.filtered)


#### 3. Filter microbe taxonomies with spatial autocorrelation ####
coords <- GetTissueCoordinates(obj)
coords <- coords[, c("x", "y")]
knn <- knn2nb(knearneigh(coords, k = 6))
w <- nb2listw(knn, style = "W")

refined_microbe_counts <- matrix(0, nrow = nrow(microbe_counts), ncol = ncol(microbe_counts))
refined_microbe_counts <- as.data.frame(refined_microbe_counts)
rownames(refined_microbe_counts) <- rownames(microbe_counts)
colnames(refined_microbe_counts) <- colnames(microbe_counts)

refined_microbe_counts[rownames(decontaminated_microbe_counts), ] <- decontaminated_microbe_counts
refined_microbe_counts <- refined_microbe_counts[, cleaned_microbes, drop = FALSE]

used_refined_microbe_counts <- refined_microbe_counts[rownames(coords), , drop = FALSE]

moran_results <- list()
for (microbe in cleaned_microbes) {
  moran_results[[microbe]] <- moran.test(used_refined_microbe_counts[, microbe], w)
}

moran_q <- p.adjust(sapply(moran_results, function(x) x$p.value), method = "BH")
significant_microbes <- names(moran_q)[moran_q < 0.05]

if (length(significant_microbes) == 0)
  quit(status = 0)


####
final_microbe_counts <- refined_microbe_counts[, significant_microbes, drop = FALSE]

obj@meta.data <- obj@meta.data[, c(colnames(obj@meta.data)[1:3],
                                   significant_microbes)]
obj@meta.data[, significant_microbes] <- final_microbe_counts
obj@meta.data[, 'Total'] <- rowSums(final_microbe_counts)


saveRDS(obj, file = paste0(output_path, sample_id, ".post_filtered.rds"))
