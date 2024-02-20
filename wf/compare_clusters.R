library(BiocManager)
library(ArchR)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(grid)
library(knitr)
library(patchwork)
library(gridExtra)
library(dplyr)
library(tibble)
library(hdf5r)
library(stringr)
library(rjson)
library(rmarkdown)
library(purrr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(EnhancedVolcano)


args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]
clusterA <- args[2]
conditionA <- args[3]
clusterB <- args[4]
conditionB <- args[5]
archr_path <- args[6]
genome <- args[7]
work_dir <- args[8]

# set genome to be used for gene and genome annotations to be mouse mm10 or
#  human hg38
addArchRGenome(genome)

setwd(work_dir)
proj_filter <- loadArchRProject(archr_path)

combine_vec <- paste(unique(proj_filter$Clusters), collapse = ",")
clusterA_list <- unlist(strsplit(clusterA, ",", fixed = TRUE))
clusterB_list <- unlist(strsplit(clusterB, ",", fixed = TRUE))
store_subsets <- c()
vector_length <- length(clusterA_list)

if (nchar(conditionA) > 0 && length(clusterA_list) > 0) {
  subsetA <- which(
    proj_filter$Condition == conditionA & proj_filter$Clusters %in%
      clusterA_list,
  )
  store_subsets <- append(store_subsets, subsetA)
} else if (nchar(conditionA) < 1 && length(clusterA_list) > 0) {
  subsetA <- which(proj_filter$Clusters %in% clusterA_list)
  store_subsets <- append(store_subsets, subsetA)
} else if (nchar(conditionA) > 0 && length(clusterA_list) < 1) {
  subsetA <- which(proj_filter$Condition == conditionA)
  store_subsets <- append(store_subsets, subsetA)
} else {
  all_indexes <- length(proj_filter$Clusters)
  subsetA <- 1:all_indexes
  store_subsets <- append(store_subsets, subsetA)
}
if (nchar(conditionB) > 0 && length(clusterB_list) > 0) {
  subsetB <- which(
    proj_filter$Condition == conditionB & proj_filter$Clusters %in%
      clusterB_list,
  )
  store_subsets <- append(store_subsets, subsetB)
} else if (nchar(conditionB) < 1 && length(clusterB_list) > 0) {
  subsetB <- which(proj_filter$Clusters %in% clusterB_list)
  store_subsets <- append(store_subsets, subsetB)
} else if (nchar(conditionB) > 0 && length(clusterB_list) < 1) {
  subsetB <- which(proj_filter$Condition == conditionB)
  store_subsets <- append(store_subsets, subsetB)
} else {
  all_indexes <- length(proj_filter$Clusters)
  subsetB <- 1:all_indexes
  store_subsets <- append(store_subsets, subsetB)
}

project_select <- proj_filter[store_subsets]

if (length(clusterA_list) < 1) {
  clusterA_list <- unlist(strsplit(combine_vec, ",", fixed = TRUE))
}
if (length(clusterB_list) < 1) {
  clusterB_list <- unlist(strsplit(combine_vec, ",", fixed = TRUE))
}

conditionA <- "ComparisonA"
conditionB <- "ComparisonB"
df <- data.frame(project_select@cellColData) %>%
  mutate(
    UpdateClustName = ifelse(
      Clusters %in% clusterB_list,
      conditionB,
      ifelse(
        Clusters %in% clusterA_list,
        conditionA,
        Clusters
      )
    )
  )

project_select$UpdateClustName <- df$UpdateClustName

groupcompare <- "UpdateClustName"


###Calculate differential genes

select_genes <- getMarkerFeatures(
  ArchRProj = project_select,
  useMatrix = "GeneScoreMatrix",
  groupBy = groupcompare,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

sample_gene_list <- getMarkers(select_genes, cutOff = "FDR <= 0.02")
write.csv(
  sample_gene_list,
  file = paste0(project_name, "_sample_gene_list.csv"),
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
  file = paste0(project_name, "_gene_markers.csv"),
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

pdf(paste0(project_name, "_", "volcano_gene.pdf"))
print(volcano)
dev.off()

marker_test <- getMarkerFeatures(
  ArchRProj = project_select,
  useMatrix = "PeakMatrix",
  groupBy = groupcompare,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = conditionA,
  bgdGroups = conditionB,
)

pma <- plotMarkers(
  seMarker = marker_test,
  name = conditionA,
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1",
  plotAs = "Volcano"
)
pdf(paste0(project_name, "_volcano_peak.pdf"))
print(pma)
dev.off()

marker_list <- getMarkers(marker_test, cutOff = "FDR <= 0.01 & Log2FC >= 1")

#Collect data with annotations
peak_data <- data.frame(
  project_select@peakSet@ranges, project_select@peakSet@elementMetadata
)
total <- merge(peak_data, marker_list, by = c("start", "end"))

write.csv(
  total, file = paste0(project_name, "_peak_markers.csv"), row.names = FALSE
)

motifs_up <- peakAnnoEnrichment(
  seMarker = marker_test,
  ArchRProj = project_select,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC > 0"
)
df <- data.frame(TF = rownames(motifs_up), mlog10Padj = assay(motifs_up)[, 1])
df <- df[order(df$mlog10Padj, decreasing = TRUE), ]
df$rank <- seq_len(nrow(df))

write.csv(df, file = paste0(project_name, "_motifsup.csv"), row.names = FALSE)

gg_up <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) +
  theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

pdf(paste0(project_name, "_UP", "motif_enrichment.pdf"))
print(gg_up)
dev.off()

motifs_do <- peakAnnoEnrichment(
    seMarker = marker_test,
    ArchRProj = project_select,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC < 0"
  )
df2 <- data.frame(TF = rownames(motifs_do), mlog10Padj = assay(motifs_do)[, 1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE), ]
df2$rank <- seq_len(nrow(df2))

write.csv(
  df2, file = paste0(project_name, "_motifsdown.csv"), row.names = FALSE
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

pdf(paste0(project_name, "_DOWN", "motif_enrichment.pdf"))
print(gg_do)
dev.off()

markers_motifs <- getMarkerFeatures(
  ArchRProj = project_select,
  useMatrix = "MotifMatrix",
  groupBy = groupcompare,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = "z",
  maxCells = 5000,
  normBy = "none"
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
  file = paste0(project_name, "_pairwise_motifs.csv"),
  row.names = FALSE
)

volcanom <- EnhancedVolcano(
  pairwise_dfm,
  lab = pairwise_dfm$pairwise_motifs,
  x = "mmean",
  y = "mpvalue",
  ylim = c(0, abs(min(log10(pairwise_dfm$mpvalue)))),
  xlim = c(-2.5, 2.5),
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

pdf(paste0(project_name, "_", "volcano_motif.pdf"))
print(volcanom)
dev.off()
