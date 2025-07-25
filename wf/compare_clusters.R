library(ArchR)
library(BiocManager)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(EnhancedVolcano)
library(ggrepel)
library(grid)
library(gridExtra)
library(hdf5r)
library(Matrix)
library(patchwork)
library(pheatmap)
library(purrr)
library(RColorBrewer)
library(Seurat)
library(SeuratObject)
library(stringr)
library(tibble)

multiple_conditions <- function(archrConditions, user) {
  # From a list of conditions, returns all containing string 'user'.

  match_cond <- c()
  user_lowercase <- tolower(user)
  pattern <- paste0("\\b", user_lowercase, "\\b")
  for (sample in archrConditions) {
    lowercase <- tolower(sample)
    result <- grepl(pattern, lowercase)
    if (result) {
      match_cond <- append(match_cond, sample)
    }
  }
  final <- paste(match_cond, collapse = ",")
  return(final)
}

args <- commandArgs(trailingOnly = TRUE)
print(args)
project_name <- args[1]
clusterA <- args[2]
conditionA <- args[3]
multipleA_flag <- args[4]
clusterB <- args[5]
conditionB <- args[6]
multipleB_flag <- args[7]
archr_path <- args[8]
genome <- args[9]
work_dir <- args[10]

addArchRGenome(genome)

setwd(work_dir)

gene_dir <- file.path(work_dir, "gene_results")
peak_dir <- file.path(work_dir, "peak_results")
motif_dir <- file.path(work_dir, "motif_results")
coverage_dir <- file.path(work_dir, "coverages")
for (dir in c(gene_dir, peak_dir, motif_dir, coverage_dir)) {
  if (!dir.exists(dir)) dir.create(dir)
}

proj_filter <- loadArchRProject(archr_path)

condition_values <- proj_filter@cellColData$Condition@values
if (multipleA_flag == "t") {  # Get all conditions containing "conditionA"
    conditionA <- multiple_conditions(condition_values, conditionA)
}
if (multipleB_flag == "t") {  # Get all conditions containing "conditionB"
    conditionB <- multiple_conditions(condition_values, conditionB)
}

combine_vec <- paste(unique(proj_filter$Clusters), collapse = ",")

clusterA_list <- unlist(strsplit(clusterA, ",", fixed = TRUE))
clusterB_list <- unlist(strsplit(clusterB, ",", fixed = TRUE))

conditionA_list <- unlist(strsplit(conditionA, ",", fixed = TRUE))
conditionB_list <- unlist(strsplit(conditionB, ",", fixed = TRUE))

# Create boolean index for Group A
if (length(conditionA_list) > 0 && length(clusterA_list) > 0) {
  subsetA <- which(
    proj_filter$Condition %in% conditionA_list & proj_filter$Clusters %in%
      clusterA_list,
  )

} else if (length(conditionA_list) == 0 && length(clusterA_list) > 0) {
  subsetA <- which(proj_filter$Clusters %in% clusterA_list)

} else if (length(conditionA_list) > 0 && length(clusterA_list) == 0) {
  subsetA <- which(proj_filter$Condition %in% conditionA_list)

} else {
  all_indexes <- length(proj_filter$Clusters)
  subsetA <- 1:all_indexes
}

# Create boolean index for Group B
if (length(conditionB_list) > 0 && length(clusterB_list) > 0) {
  subsetB <- which(
    proj_filter$Condition %in% conditionB_list & proj_filter$Clusters %in%
      clusterB_list,
  )
} else if (length(conditionB_list) == 0 && length(clusterB_list) > 0) {
  subsetB <- which(proj_filter$Clusters %in% clusterB_list)

} else if (length(conditionB) > 0 && length(clusterB_list) == 1) {
  subsetB <- which(proj_filter$Condition %in% conditionB_list)

} else {
  all_indexes <- length(proj_filter$Clusters)
  subsetB <- 1:all_indexes
}

# Label Cells in subsets with condition/group label
vector_value <- sapply(
  (1 : nrow(proj_filter@cellColData)),
  function(x) ifelse(
    x %in% subsetA,
    conditionA,
    ifelse(x %in% subsetB, conditionB, "NO")
  )
)
proj_filter$UpdateClustName <- vector_value

# Filter project to only Cells in groups subsets 
store_subsets <- sort(unique(c(subsetA, subsetB)))
project_select <- proj_filter[store_subsets]

### Calculate differential genes
select_genes <- getMarkerFeatures(
  ArchRProj = project_select,
  groupBy = "UpdateClustName",
  useMatrix = "GeneScoreMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest"
)

# Save stats for all genes
sample_gene_list <- getMarkers(select_genes, cutOff = "FDR <= 1 & Log2FC >= -Inf")
write.csv(
  sample_gene_list,
  file = file.path(gene_dir, "all_genes.csv"),
  row.names = FALSE
)

pairwise_genes <- rowData(select_genes)$name
log2FC <- assay(select_genes, "Log2FC")[, 1]
FDR <- assay(select_genes, "FDR")[, 1]
pvalue <- assay(select_genes, "Pval")[, 1]
pairwise_df <- data.frame(pairwise_genes, log2FC, pvalue, FDR)
pairwise_df <- na.omit(pairwise_df)

pairwise_df$Significance <- ifelse(
  pairwise_df$pvalue < 0.05 & abs(pairwise_df$log2FC) >= 0.4,
  ifelse(
    pairwise_df$log2FC > 0,
    colnames(assay(select_genes))[2],
    colnames(assay(select_genes))[1]
  ),
  "Not significant"
)
write.csv(
  pairwise_df,
  file = file.path(gene_dir, "marker_genes.csv"),
  row.names = FALSE
)

volcano <- EnhancedVolcano(
  pairwise_df,
  lab = pairwise_df$pairwise_genes,
  x = "log2FC",
  y = "pvalue",
  ylim = c(0, abs(min(log10(pairwise_df$pvalue)))),
  xlim =  c(-2.5, 2.5),
  title = paste0(
    colnames(assay(select_genes))[2],
    " vs ",
    colnames(assay(select_genes))[1]
  ),
  pCutoff = 0.05,
  FCcutoff = 0.4,
  pointSize = 1.0,
  labSize = 4.0
)

 
pdf(file.path(gene_dir, "volcano_gene.pdf"))
print(volcano)
dev.off()

############################# Compare Peaks ##################################

marker_test <- getMarkerFeatures(
  ArchRProj = project_select,
  groupBy = "UpdateClustName",
  useMatrix = "PeakMatrix",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

pma <- plotMarkers(
  seMarker = marker_test,
  name = conditionA,
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.4",
  plotAs = "MA"
)
pdf(file.path(peak_dir, "MA_peaks.pdf"))
print(pma)
dev.off()

marker_list <- getMarkers(marker_test, cutOff = "FDR <= 1 & Log2FC >= -Inf")

# Add annotations
peak_data <- data.frame(
  project_select@peakSet@ranges, project_select@peakSet@elementMetadata
)
total <- merge(peak_data, marker_list, by = c("start", "end"))

write.csv(
  total, file = file.path(peak_dir, "all_peaks.csv"), row.names = FALSE
)

############################# Compare Motifs #################################

motifs_up <- peakAnnoEnrichment(
  seMarker = marker_test,
  ArchRProj = project_select,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.1 & Log2FC > 0"
)
df <- data.frame(TF = rownames(motifs_up), mlog10Padj = assay(motifs_up)[, 1])
df <- df[order(df$mlog10Padj, decreasing = TRUE), ]
df$rank <- seq_len(nrow(df))

write.csv(
  df, file = file.path(motif_dir, "upRegulated_motifs.csv"), row.names = FALSE
)

gg_up <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ],
        aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) +
  theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

pdf(file.path(motif_dir, "upRegulated_motif_enrichment.pdf"))
print(gg_up)
dev.off()

motifs_do <- peakAnnoEnrichment(
  seMarker = marker_test,
  ArchRProj = project_select,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.1 & Log2FC < 0"
)
df2 <- data.frame(TF = rownames(motifs_do), mlog10Padj = assay(motifs_do)[, 1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE), ]
df2$rank <- seq_len(nrow(df2))

write.csv(
  df2, file = file.path(motif_dir, "downRegulated_motifs.csv"), row.names = FALSE
)

gg_do <- ggplot(df2, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df2[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) +
  theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

pdf(file.path(motif_dir, "downRegulated_motif_enrichment.pdf"))
print(gg_do)
dev.off()

markers_motifs <- getMarkerFeatures(
  ArchRProj = project_select,
  useMatrix = "MotifMatrix",
  groupBy = "UpdateClustName",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = "z",
  maxCells = 5000,
  normBy = "none"
)

# Save stats for all genes
motifs_list <- getMarkers(markers_motifs, cutOff = "FDR <= 1 & MeanDiff >= -Inf")
write.csv(
  motifs_list,
  file = file.path(motif_dir, "all_motifs.csv"),
  row.names = FALSE
)

pairwise_motifs <- rowData(markers_motifs)$name
mmean <- assay(markers_motifs, "MeanDiff")[, 1]
mFDR <- assay(markers_motifs, "FDR")[, 1]
mpvalue <- assay(markers_motifs, "Pval")[, 1]
pairwise_dfm <- data.frame(pairwise_motifs, mmean, mpvalue, mFDR)

pairwise_dfm$Significance <- ifelse(
  pairwise_dfm$mpvalue < 0.05 & abs(pairwise_dfm$mmean) >= 0.4,
  ifelse(
    pairwise_df$mpvalue > 0,
    colnames(assay(markers_motifs))[2],
    colnames(assay(markers_motifs))[1]
  ),
  "Not significant")
pairwise_dfm <- na.omit(pairwise_dfm)
write.csv(
  pairwise_dfm,
  file = file.path(motif_dir, "marker_motifs.csv"),
  row.names = FALSE
)

volcanom <- EnhancedVolcano(
  pairwise_dfm,
  lab = pairwise_dfm$pairwise_motifs,
  x = "mmean",
  y = "mpvalue",
  ylim = c(0, abs(min(log10(pairwise_dfm$mpvalue)))),
  xlim = c(-2.5, 2.5),
  xlab = bquote("MeanDiff"),
  title = paste0(
    colnames(assay(markers_motifs))[2],
    " vs ",
    colnames(assay(markers_motifs))[1]
  ),
  pCutoff = 0.05,
  FCcutoff = 0.4,
  pointSize = 1.0,
  labSize = 4.0
)

pdf(file.path(motif_dir, "volcano_motif.pdf"))
print(volcanom)
dev.off()

# Coverage files
file_names <- getGroupBW(
  ArchRProj = project_select,
  groupBy = "UpdateClustName",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
)

for (file_name in file_names) {
  file.copy(from = file_name, to = coverage_dir)  
}

saveArchRProject(ArchRProj = proj_filter)
