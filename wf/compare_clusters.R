library(BiocManager)
library(devtools)
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
library(dplyr)
library(ggpubr)
library(EnhancedVolcano)


# globals ---------------------------------------------------------------------
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

clusterA_list <- unlist(strsplit(clusterA, ",", fixed = TRUE))
clusterB_list <- unlist(strsplit(clusterB, ",", fixed = TRUE))
store_subsets <- c()
vector_length <- length(clusterA_list)

if (nchar(conditionA) > 1) {
  subsetA <- which(
    proj_filter$Condition == conditionA & proj_filter$Clusters %in% clusterA_list, 
  )
  store_subsets <- append(store_subsets, subsetA)
} else {
  subsetA <- which(proj_filter$Clusters %in% clusterA_list)
  store_subsets <- append(store_subsets, subsetA)
}

if (nchar(conditionB) > 1) {
  subsetB <- which(
    proj_filter$Condition == conditionB & proj_filter$Clusters %in% clusterB_list, 
  )
  store_subsets <- append(store_subsets, subsetB)
} else {
  subsetB <- which(proj_filter$Clusters %in% clusterB_list)
  store_subsets <- append(store_subsets, subsetB)
}

project_select <- proj_filter[store_subsets]

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

Select_genes <- getMarkerFeatures(
  ArchRProj = project_select,
  useMatrix = "GeneScoreMatrix",
  groupBy = groupcompare,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

SampleGeneList <- getMarkers(Select_genes, cutOff = "FDR <= 0.02")
write.csv(
  SampleGeneList,
  file = paste0(project_name, "_Sample_gene_list.csv"),
  row.names = FALSE
)

pairwise_Genes <- rowData(Select_genes)$name
log2FC <- assay(Select_genes, "Log2FC")[, 1]
FDR <- assay(Select_genes, "FDR")[, 1]
pvalue <- assay(Select_genes, "Pval")[, 1]
pairwise_df <- data.frame(pairwise_Genes, log2FC, pvalue, FDR)
pairwise_df <- na.omit(pairwise_df)
pairwise_df$Significance <- ifelse(
  pairwise_df$pvalue < 0.05 & abs(pairwise_df$log2FC) >= 0.4,
  ifelse(
    pairwise_df$log2FC > 0,
    colnames(assay(Select_genes))[2],
    colnames(assay(Select_genes))[1]
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
  lab = pairwise_df$pairwise_Genes,
  x = "log2FC",
  y = "pvalue",
  ylim = c(0, abs(min(log10(pairwise_df$pvalue)))),
  xlim =  c(-2.5, 2.5),
  title = paste0(
    colnames(assay(Select_genes))[2],
    " vs ",
    colnames(assay(Select_genes))[1]
  ),
  pCutoff = 0.05,
  FCcutoff = 0.4,
  pointSize = 1.0,
  labSize = 4.0
)

volcano
pdf(paste0(project_name, "_", "volcano_gene.pdf"))
print(volcano)
dev.off()

markerTest <- getMarkerFeatures(
  ArchRProj = project_select,
  useMatrix = "PeakMatrix",
  groupBy = groupcompare,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = conditionA,
  bgdGroups = conditionB,
)

pma <- plotMarkers(
  seMarker = markerTest,
  name = conditionA,
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1",
  plotAs = "Volcano"
)
pdf(paste0(project_name, "_volcano_peak.pdf"))
print(pma)
dev.off()

markerList <- getMarkers(markerTest, cutOff = "FDR <= 0.01 & Log2FC >= 1")

#Collect data with annotations
peak_data <- data.frame(
  project_select@peakSet@ranges, project_select@peakSet@elementMetadata
)
total <- merge(peak_data, markerList, by = c("start", "end"))

write.csv(
  total, file = paste0(project_name, "_peak_markers.csv"), row.names = FALSE
)

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = project_select,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC > 0"
)
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[, 1])
df <- df[order(df$mlog10Padj, decreasing = TRUE), ]
df$rank <- seq_len(nrow(df))

write.csv(df, file = paste0(project_name, "_motifsup.csv"), row.names = FALSE)

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
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

ggUp

pdf(paste0(project_name, "_UP", "motif_enrichment.pdf"))
print(ggUp)
dev.off()

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = project_select,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC < 0"
  )
df2 <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[, 1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE), ]
df2$rank <- seq_len(nrow(df2))

write.csv(
  df2, file = paste0(project_name, "_motifsdown.csv"), row.names = FALSE
)

ggDo <- ggplot(df2, aes(rank, mlog10Padj, color = mlog10Padj)) +
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

ggDo

pdf(paste0(project_name, "_DOWN", "motif_enrichment.pdf"))
print(ggDo)
dev.off()

markersMotifs <- getMarkerFeatures(
  ArchRProj = project_select,
  useMatrix = "MotifMatrix",
  groupBy = groupcompare,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = "z",
  maxCells = 5000,
  normBy = "none"
)

pairwise_motifs <- rowData(markersMotifs)$name
mmean <- assay(markersMotifs, "MeanDiff")[, 1]
mFDR <- assay(markersMotifs, "FDR")[, 1]
mpvalue <- assay(markersMotifs, "Pval")[, 1]
pairwise_dfm <- data.frame(pairwise_motifs, mmean, mpvalue, mFDR)

pairwise_dfm$Significance <- ifelse(
  pairwise_dfm$mpvalue < 0.05 & abs(pairwise_dfm$mmean) >= 0.4,
  ifelse(
    pairwise_df$mpvalue> 0,
    colnames(assay(markersMotifs))[2],
    colnames(assay(markersMotifs))[1]
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
    colnames(assay(markersMotifs))[2], " vs ", colnames(assay(markersMotifs))[1]
  ),
  pCutoff = 0.05,
  FCcutoff = 0.4,
  pointSize = 1.0,
  labSize = 4.0
)

volcanom
pdf(paste0(project_name, "_", "volcano_motif.pdf"))
print(volcanom)
dev.off()
