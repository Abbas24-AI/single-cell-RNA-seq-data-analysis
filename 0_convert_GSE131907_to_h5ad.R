suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(future)
  library(ggplot2)
  library(dplyr)
})

base_dir <- "/home/khan/nlsc_scrna"
seurat_rds <- file.path(base_dir, "NSCLC_Seurat.rds")

if (!file.exists(seurat_rds)) {
  stop("NSCLC_Seurat.rds not found at: ", seurat_rds)
}

# parallel (Linux)
n_cores <- max(1, parallel::detectCores() - 2)
plan(multicore, workers = n_cores)
options(future.globals.maxSize = 20000 * 1024^2)

cat("Base dir :", base_dir, "\n")
cat("Cores    :", n_cores, "\n\n")

# Output dirs
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################
# 0) Paths & parallel setup
###############################################################

base_dir <- "/home/khan/nlsc_scrna"
seurat_rds <- file.path(base_dir, "NSCLC_Seurat.rds")

if (!file.exists(seurat_rds)) {
  stop("NSCLC_Seurat.rds not found at: ", seurat_rds)
}

# parallel (Linux)
n_cores <- max(1, parallel::detectCores() - 2)
plan(multicore, workers = n_cores)
options(future.globals.maxSize = 20000 * 1024^2)

cat("Base dir :", base_dir, "\n")
cat("Cores    :", n_cores, "\n\n")

# Output dirs
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################
# 1) Load Seurat object
###############################################################

cat("[1] Loading NSCLC_Seurat.rds...\n")
sce <- readRDS(seurat_rds)

DefaultAssay(sce) <- "RNA"
cat("Cells:", ncol(sce), " | Genes:", nrow(sce), "\n\n")

cat("Metadata columns:\n")
print(colnames(sce@meta.data))

###############################################################
# 2) Quick QC overview (like Fig 1 extended)
###############################################################

library(ggplot2)

df <- sce@meta.data
df$Sample_Origin <- factor(df$Sample_Origin)

# Single metric violin: percent.mt
p_mt <- ggplot(df, aes(x = Sample_Origin, y = percent.mt, fill = Sample_Origin)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.1, color = "black") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  ylab("Percent mitochondrial reads") +
  ggtitle("Mitochondrial content by Sample_Origin")

ggsave("QC_percent_mt_by_origin.png", p_mt, width = 6, height = 4, dpi = 300)

metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

for (m in metrics) {
  p <- ggplot(df, aes(x = Sample_Origin, y = .data[[m]], fill = Sample_Origin)) +
    geom_violin(trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.1, color = "black") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    ylab(m) +
    ggtitle(paste(m, "by Sample_Origin"))
  
  ggsave(paste0("QC_", m, "_by_origin.png"),
         p, width = 6, height = 4, dpi = 300)
}
###############################################################
# 3) Normalization & variance modeling
#    Option A: SCTransform (recommended)
#    Option B: LogNormalize (fallback)
###############################################################

use_sctransform <- TRUE  # set FALSE if memory is an issue

if (use_sctransform) {
  cat("[3] Running SCTransform (this may take some time)...\n")
  # Regress out mitochondrial percentage, keep 3000 HVGs
  sce <- SCTransform(
    sce,
    vst.flavor = "v2",
    variable.features.n = 3000,
    vars.to.regress = "percent.mt",
    verbose = TRUE
  )
  DefaultAssay(sce) <- "SCT"
} else {
  cat("[3] Using standard log-normalization...\n")
  sce <- NormalizeData(sce, verbose = FALSE)
  sce <- FindVariableFeatures(sce, nfeatures = 3000, selection.method = "vst")
  sce <- ScaleData(sce, vars.to.regress = "percent.mt", verbose = FALSE)
  DefaultAssay(sce) <- "RNA"
}

###############################################################
# 4) Dimensionality reduction: PCA + UMAP
###############################################################

cat("[4] Running PCA...\n")

sce <- RunPCA(sce, npcs = 50, verbose = FALSE)

# inspect first PCs
p_pca <- DimPlot(sce, reduction = "pca", group.by = "Sample_Origin") +
  ggtitle("PCA colored by Sample_Origin")

ggsave(file.path(fig_dir, "PCA_by_origin.png"),
       p_pca, width = 6, height = 5, dpi = 300)

cat("[5] Running neighbors & UMAP (dims 1–30)...\n")

sce <- FindNeighbors(sce, dims = 1:30, verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:30, verbose = FALSE)

p_umap_origin <- DimPlot(sce, reduction = "umap", group.by = "Sample_Origin") +
  ggtitle("UMAP colored by Sample_Origin")

ggsave(file.path(fig_dir, "UMAP_by_Sample_Origin.png"),
       p_umap_origin, width = 7, height = 6, dpi = 300)
###############################################################
# 5) Clustering
#    Paper-level analysis often uses multiple resolutions
###############################################################

cat("[6] Graph-based clustering...\n")

sce <- FindClusters(sce, resolution = 0.6)   # broad
sce <- FindClusters(sce, resolution = 1.0, verbose = FALSE)  # finer

# store preferred cluster column
sce$cluster_res0.6 <- sce$SCT_snn_res.0.6 %||% sce$RNA_snn_res.0.6
sce$cluster_res1.0 <- sce$SCT_snn_res.1  %||% sce$RNA_snn_res.1

# UMAP colored by cluster (0.6)
p_umap_clust <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "cluster_res0.6",
  label = TRUE,
  repel = TRUE
) + ggtitle("UMAP – clusters (resolution 0.6)")

ggsave(file.path(fig_dir, "UMAP_by_cluster_res0.6.png"),
       p_umap_clust, width = 7, height = 6, dpi = 300)
###############################################################
# 6) Transfer / Use published cell type annotation
###############################################################

# You already have:
#   Cell_type           = broad lineage (Myeloid, T, B, Epithelial, Fibroblast...)
#   Cell_type.refined   = refined type
#   Cell_subtype        = subclusters like mo-Mac, Treg, etc.

cat("[7] Using provided annotation for cell types...\n")

# sanity check
table_CellType <- table(sce$Cell_type)
cat("Broad Cell_type distribution:\n")
print(table_CellType)

# UMAP by broad Cell_type
p_umap_celltype <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "Cell_type",
  label = TRUE,
  repel = TRUE
) + ggtitle("UMAP – Broad Cell Types")

ggsave(file.path(fig_dir, "UMAP_by_Cell_type.png"),
       p_umap_celltype, width = 7, height = 6, dpi = 300)

# UMAP by refined subtype
p_umap_subtype <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "Cell_type.refined",
  label = FALSE
) + ggtitle("UMAP – Refined Cell Types")

ggsave(file.path(fig_dir, "UMAP_by_Cell_type_refined.png"),
       p_umap_subtype, width = 8, height = 6, dpi = 300)
###############################################################
# 7) Composition plots (like atlas barplots)
###############################################################

cat("[8] Building composition plots...\n")

md <- sce@meta.data

# Example 1: fraction of Cell_type per Sample_Origin (nLung, Tumor, LN, etc.)
comp1 <- md %>%
  group_by(Sample_Origin, Cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample_Origin) %>%
  mutate(freq = n / sum(n))

p_comp1 <- ggplot(comp1, aes(x = Sample_Origin, y = freq, fill = Cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 12) +
  ggtitle("Cell-type composition across Sample_Origin") +
  xlab("Sample Origin") + ylab("Proportion of cells")

ggsave(file.path(fig_dir, "Composition_SampleOrigin_CellType.png"),
       p_comp1, width = 7, height = 5, dpi = 300)

# Example 2: composition per Sample (patient-level)
comp2 <- md %>%
  group_by(Sample, Cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(freq = n / sum(n))

p_comp2 <- ggplot(comp2, aes(x = Sample, y = freq, fill = Cell_type)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Proportion of cells") +
  ggtitle("Cell-type composition per sample")

ggsave(file.path(fig_dir, "Composition_Sample_CellType.png"),
       p_comp2, width = 10, height = 5, dpi = 300)
###############################################################
# 8) Marker gene visualization (Nature-style UMAP markers)
###############################################################

cat("[9] Marker expression plots...\n")

# canonical markers for major compartments
marker_list <- list(
  Epithelial = c("EPCAM", "KRT8", "KRT18"),
  Myeloid    = c("LYZ", "CD68", "LST1"),
  Tcells     = c("CD3D", "CD3E", "CD4", "CD8A"),
  Bcells     = c("MS4A1", "CD79A"),
  Fibroblast = c("COL1A1", "DCN"),
  Endothelial = c("PECAM1", "VWF")
)

# UMAP FeaturePlots for a few key markers
fp1 <- FeaturePlot(sce, reduction = "umap",
                   features = c("EPCAM", "KRT8", "LYZ", "CD3D"),
                   ncol = 2, order = TRUE)
ggsave(file.path(fig_dir, "UMAP_markers_epithelium_immune.png"),
       fp1, width = 8, height = 8, dpi = 300)

# DotPlot summarizing markers by Cell_type
all_markers <- unique(unlist(marker_list))
dot1 <- DotPlot(sce, features = all_markers, group.by = "Cell_type") +
  RotatedAxis() +
  ggtitle("Marker expression by broad Cell_type")

ggsave(file.path(fig_dir, "DotPlot_markers_by_Cell_type.png"),
       dot1, width = 9, height = 5, dpi = 300)

###############################################################
# 9) Example: Differential expression (Tumor vs Normal in epithelial)
###############################################################

cat("[10] Example DE: tumor vs normal epithelium...\n")

# Focus on epithelial cells only
epi <- subset(sce, subset = Cell_type == "Epithelial cells" | Cell_type == "Epithelial")

# Define tumor vs normal based on Sample_Origin (adjust if different labels)
table(epi$Sample_Origin)

# Example: Tumor = tLung / primary tumor; Normal = nLung
Idents(epi) <- factor(epi$Sample_Origin)

if (all(c("tLung", "nLung") %in% levels(Idents(epi)))) {
  de_res <- FindMarkers(epi, ident.1 = "tLung", ident.2 = "nLung",
                        logfc.threshold = 0.25, min.pct = 0.1)
  de_res <- de_res[order(de_res$p_val_adj), ]
  write.csv(de_res,
            file.path(base_dir, "DE_tumor_vs_normal_epithelium.csv"),
            quote = TRUE)
  cat("DE results saved: DE_tumor_vs_normal_epithelium.csv\n")
} else {
  cat("Warning: Could not find both 'tLung' and 'nLung' in Sample_Origin for epithelial subset.\n")
}

###############################################################
# 10) Save updated object
###############################################################

saveRDS(sce, file.path(base_dir, "NSCLC_Seurat_processed.rds"))
cat("\n✅ DONE: NSCLC_Seurat_processed.rds saved.\n")
cat("Figures written to: ", fig_dir, "\n")

#!/usr/bin/env Rscript

###############################################################
# 2_epithelial_subatlas.R
# NSCLC (GSE131907) – Epithelial sub-atlas in R / Seurat
###############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(future)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

###############################################################
# 0) Setup paths & parallelization
###############################################################

base_dir <- "/home/khan/nlsc_scrna"

# Use the processed object from your previous script
seurat_in  <- file.path(base_dir, "NSCLC_Seurat_processed.rds")
epi_out_rds <- file.path(base_dir, "NSCLC_Epithelial.rds")
fig_dir   <- file.path(base_dir, "figures_epithelial")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(seurat_in)) {
  stop("NSCLC_Seurat_processed.rds not found at: ", seurat_in,
       "\nRun your global script first.")
}

n_cores <- max(1, parallel::detectCores() - 2)
plan(multicore, workers = n_cores)
options(future.globals.maxSize = 20000 * 1024^2)

cat("Base dir :", base_dir, "\n")
cat("Cores    :", n_cores, "\n\n")

###############################################################
# 1) Load full object & subset epithelial cells
###############################################################

cat("[1] Loading full NSCLC object...\n")
sce <- readRDS(seurat_in)
DefaultAssay(sce) <- "RNA"

cat("Total cells:", ncol(sce), "| genes:", nrow(sce), "\n")
cat("Cell_type levels:\n")
print(table(sce$Cell_type))

# Epithelial labels differ sometimes; handle both
epi <- subset(
  sce,
  subset = Cell_type %in% c("Epithelial cells", "Epithelial", "Epithelial cell")
)

cat("Epithelial cells:", ncol(epi), " | genes:", nrow(epi), "\n\n")

# Ensure basic metadata
epi$Sample_Origin <- factor(epi$Sample_Origin)
epi$Sample        <- factor(epi$Sample)

###############################################################
# 2) QC summary for epithelial only (optional but useful)
###############################################################

cat("[2] QC summary for epithelial sub-atlas...\n")

df_epi <- epi@meta.data

qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

for (m in qc_metrics) {
  p <- ggplot(df_epi, aes(x = Sample_Origin, y = .data[[m]], fill = Sample_Origin)) +
    geom_violin(trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.1, outlier.size = 0.1, color = "black") +
    scale_y_continuous(labels = comma) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    ylab(m) +
    ggtitle(paste("Epithelial:", m, "by Sample_Origin"))
  
  ggsave(file.path(fig_dir, paste0("Epi_QC_", m, "_by_origin.png")),
         p, width = 6, height = 4, dpi = 300)
}

cat("Epithelial QC plots written to:", fig_dir, "\n\n")

###############################################################
# 3) Epithelial-specific normalization & clustering
###############################################################

use_sctransform <- TRUE

if (use_sctransform) {
  cat("[3] Running SCTransform on epithelial cells...\n")
  epi <- SCTransform(
    epi,
    vst.flavor = "v2",
    variable.features.n = 3000,
    vars.to.regress = "percent.mt",
    verbose = TRUE
  )
  DefaultAssay(epi) <- "SCT"
} else {
  cat("[3] Using log-normalization on epithelial cells...\n")
  epi <- NormalizeData(epi, verbose = FALSE)
  epi <- FindVariableFeatures(epi, nfeatures = 3000)
  epi <- ScaleData(epi, vars.to.regress = "percent.mt", verbose = FALSE)
  DefaultAssay(epi) <- "RNA"
}

cat("[4] PCA + UMAP on epithelial cells...\n")
epi <- RunPCA(epi, npcs = 50, verbose = FALSE)
epi <- FindNeighbors(epi, dims = 1:30, verbose = FALSE)
epi <- FindClusters(epi, resolution = 0.8, verbose = FALSE)
epi <- RunUMAP(epi, dims = 1:30, verbose = FALSE)

epi$epi_cluster <- epi$SCT_snn_res.0.8 %||% epi$RNA_snn_res.0.8

p_umap_epi <- DimPlot(
  epi,
  reduction = "umap",
  group.by = "epi_cluster",
  label = TRUE,
  repel = TRUE
) + ggtitle("Epithelial sub-atlas – clusters (res=0.8)")

ggsave(file.path(fig_dir, "Epi_UMAP_clusters.png"),
       p_umap_epi, width = 7, height = 6, dpi = 300)

# UMAP by sample origin (tumor vs normal)
p_umap_origin <- DimPlot(
  epi,
  reduction = "umap",
  group.by = "Sample_Origin"
) + ggtitle("Epithelial UMAP – Sample_Origin")

ggsave(file.path(fig_dir, "Epi_UMAP_by_Sample_Origin.png"),
       p_umap_origin, width = 7, height = 6, dpi = 300)
###############################################################
# 4) Define tumor vs normal epithelial states
###############################################################

cat("[5] Defining TumorStatus from Sample_Origin...\n")

# This depends on how Sample_Origin is coded in this dataset.
# Common values in GSE131907: nLung (normal), tLung (tumor), LN, Effusion, Broncho, etc.

table(epi$Sample_Origin)

epi$TumorStatus <- "Other"
epi$TumorStatus[epi$Sample_Origin == "nLung"] <- "Normal"
epi$TumorStatus[epi$Sample_Origin == "tLung"] <- "Tumor"

epi$TumorStatus <- factor(epi$TumorStatus, levels = c("Normal", "Tumor", "Other"))

cat("TumorStatus distribution:\n")
print(table(epi$TumorStatus))

p_umap_ts <- DimPlot(
  epi,
  reduction = "umap",
  group.by = "TumorStatus"
) + ggtitle("Epithelial UMAP – TumorStatus")

ggsave(file.path(fig_dir, "Epi_UMAP_by_TumorStatus.png"),
       p_umap_ts, width = 7, height = 6, dpi = 300)
###############################################################
# 5) Gene program scoring – EMT, proliferation, interferon, LUAD/LUSC-like
###############################################################

cat("[6] Scoring gene programs (EMT, Proliferation, IFN, LUAD/LUSC-like)...\n")

DefaultAssay(epi) <- "SCT"

emt_genes    <- c("VIM","ZEB1","ZEB2","SNAI1","SNAI2","FN1","ITGA5","ITGB1")
prolif_genes <- c("MKI67","TOP2A","PCNA","BIRC5","CCNB1","CDC20","AURKA","AURKB")
ifn_genes    <- c("IFIT1","IFIT3","ISG15","MX1","OAS1","STAT1","IRF7")

luad_genes <- c("EPCAM","KRT8","KRT18","KRT19","MUC1")
lusc_genes <- c("KRT5","KRT6A","KRT14","TP63","DSG3")

epi <- AddModuleScore(epi, features = list(emt_genes),    name = "EMT_Score")
epi <- AddModuleScore(epi, features = list(prolif_genes), name = "Prolif_Score")
epi <- AddModuleScore(epi, features = list(ifn_genes),    name = "IFN_Score")
epi <- AddModuleScore(epi, features = list(luad_genes),   name = "LUAD_Score")
epi <- AddModuleScore(epi, features = list(lusc_genes),   name = "LUSC_Score")

score_features <- c(
  "EMT_Score1",
  "Prolif_Score1",
  "IFN_Score1",
  "LUAD_Score1",
  "LUSC_Score1"
)
for (feat in score_features) {
  
  p <- FeaturePlot(
    epi,
    reduction = "umap",
    features = feat,
    cols = c("grey90", "darkred"),
    min.cutoff = "q5",
    max.cutoff = "q95"
  ) +
    ggtitle(paste("Epithelial", feat)) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    filename = file.path(fig_dir, paste0("Epi_UMAP_", feat, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}
###############################################################
# 6) Marker-based refinement of epithelial clusters
###############################################################

cat("[7] Marker expression per epithelial cluster...\n")

# Quickly find markers across epithelial clusters (res 0.8)
Idents(epi) <- epi$epi_cluster

markers_epi <- FindAllMarkers(
  epi,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

markers_epi <- markers_epi[order(markers_epi$avg_log2FC, decreasing = TRUE), ]

write.csv(markers_epi,
          file.path(base_dir, "Epi_cluster_markers_res0.8.csv"),
          row.names = FALSE)

# Heatmap of top markers per cluster
top_markers <- markers_epi %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

DefaultAssay(epi) <- "SCT"
epi <- ScaleData(epi, features = top_markers, verbose = FALSE)

p_heat <- DoHeatmap(
  epi,
  features = top_markers,
  group.by = "epi_cluster",
  size = 3
) + ggtitle("Top 10 markers per epithelial cluster")

ggsave(file.path(fig_dir, "Epi_heatmap_top_markers.png"),
       p_heat, width = 10, height = 8, dpi = 300)
###############################################################
# 7) Tumor vs Normal epithelial DE analysis (like paper)
###############################################################

cat("[8] Differential expression: Tumor vs Normal epithelium...\n")

epi_sub <- subset(epi, subset = TumorStatus %in% c("Normal", "Tumor"))

Idents(epi_sub) <- epi_sub$TumorStatus

if (all(c("Normal","Tumor") %in% levels(Idents(epi_sub)))) {
  de_epi <- FindMarkers(
    epi_sub,
    ident.1 = "Tumor",
    ident.2 = "Normal",
    logfc.threshold = 0.25,
    min.pct = 0.1
  )
  de_epi <- de_epi[order(de_epi$p_val_adj),]
  write.csv(de_epi,
            file.path(base_dir, "Epi_DE_Tumor_vs_Normal.csv"),
            row.names = TRUE)
  cat("DE table saved: Epi_DE_Tumor_vs_Normal.csv\n")
  
  # Quick volcano-style scatter (not real volcano but close)
  de_df <- de_epi
  de_df$gene <- rownames(de_df)
  
  p_de <- ggplot(de_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(alpha = 0.4) +
    theme_bw(base_size = 12) +
    xlab("log2 fold-change (Tumor vs Normal)") +
    ylab("-log10 adjusted p-value") +
    ggtitle("Tumor vs Normal epithelial – DE summary")
  
  ggsave(file.path(fig_dir, "Epi_DE_Tumor_vs_Normal_scatter.png"),
         p_de, width = 6, height = 5, dpi = 300)
  
} else {
  cat("WARNING: Could not find both 'Normal' and 'Tumor' levels in TumorStatus for epithelium.\n")
}

###############################################################
# 8) Save epithelial object
###############################################################

saveRDS(epi, epi_out_rds)
cat("\n✅ DONE: Epithelial sub-atlas saved as:\n  ", epi_out_rds, "\n")
cat("Figures written to:\n  ", fig_dir, "\n")

#!/usr/bin/env Rscript

###############################################################
# 3_epithelial_states_trajectory.R
# Refine epithelial states, compare EMT/Prolif, and run pseudotime
###############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(future)
})

base_dir <- "/home/khan/nlsc_scrna"
epi_rds  <- file.path(base_dir, "NSCLC_Epithelial.rds")
fig_dir  <- file.path(base_dir, "figures_epithelial_states")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(epi_rds)) {
  stop("NSCLC_Epithelial.rds not found at: ", epi_rds)
}

n_cores <- max(1, parallel::detectCores() - 2)
plan(multicore, workers = n_cores)
options(future.globals.maxSize = 20000 * 1024^2)

cat("Base dir:", base_dir, "\n")
cat("Cores   :", n_cores, "\n\n")

###############################################################
# 1) Load epithelial object
###############################################################

cat("[1] Loading epithelial sub-atlas...\n")
epi <- readRDS(epi_rds)
DefaultAssay(epi) <- "SCT"

cat("Cells:", ncol(epi), "| Genes:", nrow(epi), "\n")
cat("Available metadata:\n")
print(colnames(epi@meta.data))

if (is.null(epi$epi_cluster)) {
  stop("epi_cluster column not found. Make sure you ran 2_epithelial_subatlas.R.")
}
###############################################################
# 2) Cluster-level summaries of module scores
###############################################################

cat("[2] Computing cluster-level summaries (EMT / Prolif / IFN / LUAD / LUSC)...\n")

scores <- c("EMT_Score1","Prolif_Score1","IFN_Score1","LUAD_Score1","LUSC_Score1")

if (!all(scores %in% colnames(epi@meta.data))) {
  stop("Some module scores are missing. Ensure 2_epithelial_subatlas.R was run.")
}

meta_epi <- epi@meta.data %>%
  mutate(epi_cluster = as.factor(epi_cluster))

cluster_summary <- meta_epi %>%
  group_by(epi_cluster) %>%
  summarise(
    n_cells      = n(),
    EMT_mean     = mean(EMT_Score1, na.rm = TRUE),
    Prolif_mean  = mean(Prolif_Score1, na.rm = TRUE),
    IFN_mean     = mean(IFN_Score1, na.rm = TRUE),
    LUAD_mean    = mean(LUAD_Score1, na.rm = TRUE),
    LUSC_mean    = mean(LUSC_Score1, na.rm = TRUE)
  ) %>%
  arrange(epi_cluster)

write.csv(cluster_summary,
          file.path(base_dir, "Epi_cluster_program_summary.csv"),
          row.names = FALSE)

print(cluster_summary)

# Heatmap-like plot of program means per cluster (long format)
cluster_long <- cluster_summary %>%
  tidyr::pivot_longer(
    cols = ends_with("_mean"),
    names_to = "Program",
    values_to = "MeanScore"
  )

p_prog_heat <- ggplot(cluster_long,
                      aes(x = epi_cluster, y = Program, fill = MeanScore)) +
  geom_tile(color = "grey60") +
  scale_fill_gradient2(low = "navy", mid = "white", high = "darkred") +
  theme_bw(base_size = 11) +
  ggtitle("Average module scores per epithelial cluster") +
  xlab("Epithelial cluster") + ylab("Program")

ggsave(file.path(fig_dir, "Epi_cluster_program_heatmap.png"),
       p_prog_heat, width = 7, height = 4.5, dpi = 300)
###############################################################
# 3) Assign biologically meaningful labels (manual mapping)
###############################################################
# NOTE:
# You MUST look at `Epi_cluster_program_summary.csv` and marker genes
# to refine this mapping. Below is a TEMPLATE mapping you should adjust.

cat("[3] Creating example biological labels for epithelial clusters...\n")

# Get existing cluster IDs
clusters <- levels(meta_epi$epi_cluster)
print(clusters)

# ---- EDIT THIS VECTOR AFTER INSPECTING cluster_summary & markers ----
# Example placeholder labels with same length as number of clusters:
# e.g. c("Proliferative", "EMT-high", "Basal-like", "AT2-like", ...)
new_labels <- setNames(
  object = paste0("State_", clusters),
  nm     = clusters
)

# Example for 6 clusters:
# new_labels <- c(
#   `0` = "AT2-like (LUAD)",
#   `1` = "EMT-high",
#   `2` = "Basal-like (LUSC-like)",
#   `3` = "Proliferative",
#   `4` = "IFN-high",
#   `5` = "Normal-like"
# )

if (length(new_labels) != length(clusters)) {
  stop("new_labels length does not match number of clusters. Edit mapping.")
}

epi$epi_state <- plyr::mapvalues(
  x = epi$epi_cluster,
  from = names(new_labels),
  to   = new_labels
)

epi$epi_state <- factor(epi$epi_state)

p_umap_state <- DimPlot(epi, reduction = "umap",
                        group.by = "epi_state",
                        label = TRUE, repel = TRUE) +
  ggtitle("Epithelial states (manual labels)")

ggsave(file.path(fig_dir, "Epi_UMAP_epi_state.png"),
       p_umap_state, width = 7, height = 6, dpi = 300)

###############################################################
# 4) EMT / Proliferation across clusters & patients
###############################################################
library(ggplot2)
library(dplyr)

features_to_plot <- c("EMT_Score1", "Prolif_Score1", "IFN_Score1")

plot_df <- FetchData(
  epi,
  vars = c(features_to_plot, "epi_state")
)
for (feat in features_to_plot) {
  
  p <- ggplot(
    plot_df,
    aes(x = epi_state, y = .data[[feat]], fill = epi_state)
  ) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.6) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = paste(feat, "by epithelial state"),
      x = "Epithelial state",
      y = "Module score"
    )
  
  ggsave(
    filename = file.path(fig_dir, paste0("Epi_Vln_", feat, "_by_state.png")),
    plot = p,
    width = 6.5,
    height = 4,
    dpi = 300
  )
}
###############################################################
# 5) Trajectory analysis using monocle3
###############################################################

cat("[5] Running monocle3 trajectory (pseudotime)...\n")

suppressPackageStartupMessages({
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    stop("Package 'monocle3' is not installed. Install it before running this section.")
  }
  if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
    stop("Package 'SeuratWrappers' is not installed. Install it before running this section.")
  }
  library(monocle3)
  library(SeuratWrappers)
})

# Convert to CellDataSet
cds <- as.cell_data_set(epi)

# Use the same dimensional reduction (UMAP) if available
# Ensure cluster and partition fields exist
cds <- cluster_cells(cds, reduction = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# Define root cells = Normal, low EMT epithelial cells
meta_cds <- as.data.frame(colData(cds))

if (!"TumorStatus" %in% colnames(meta_cds)) {
  stop("TumorStatus not found on cds. Ensure epi$TumorStatus is set.")
}

emt_cut <- quantile(meta_cds$EMT_Score1, probs = 0.2, na.rm = TRUE)

root_candidates <- rownames(meta_cds)[
  meta_cds$TumorStatus == "Normal" &
    meta_cds$EMT_Score1 <= emt_cut
]

if (length(root_candidates) < 50) {
  warning("Few root candidates; using all Normal epithelium as root.")
  root_candidates <- rownames(meta_cds)[meta_cds$TumorStatus == "Normal"]
}

cds <- order_cells(cds, root_cells = root_candidates)

# Extract pseudotime and add back into Seurat
pseudotime_vals <- pseudotime(cds)
epi$Pseudotime <- pseudotime_vals[colnames(epi)]

# Plot pseudotime on UMAP (Seurat)
p_ps <- FeaturePlot(
  epi,
  reduction = "umap",
  features = "Pseudotime",
  cols = c("navy", "gold", "red"),
  min.cutoff = "q1",
  max.cutoff = "q99"
) + ggtitle("Epithelial pseudotime")

ggsave(file.path(fig_dir, "Epi_UMAP_pseudotime.png"),
       p_ps, width = 7, height = 6, dpi = 300)

# Relationship of pseudotime with EMT & Prolif
p_ps_emt <- ggplot(epi@meta.data,
                   aes(x = Pseudotime, y = EMT_Score1)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", color = "red") +
  theme_bw(base_size = 11) +
  xlab("Pseudotime") + ylab("EMT score") +
  ggtitle("EMT dynamics along pseudotime")

ggsave(file.path(fig_dir, "Epi_Pseudotime_vs_EMT.png"),
       p_ps_emt, width = 6, height = 4.5, dpi = 300)

p_ps_prolif <- ggplot(epi@meta.data,
                      aes(x = Pseudotime, y = Prolif_Score1)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess", color = "blue") +
  theme_bw(base_size = 11) +
  xlab("Pseudotime") + ylab("Proliferation score") +
  ggtitle("Proliferation dynamics along pseudotime")

ggsave(file.path(fig_dir, "Epi_Pseudotime_vs_Prolif.png"),
       p_ps_prolif, width = 6, height = 4.5, dpi = 300)
###############################################################
# 5) Trajectory analysis using monocle3
###############################################################

###############################################################
# 5) Trajectory analysis using monocle3
###############################################################

cat("[5] Running monocle3 trajectory (pseudotime)...\n")

suppressPackageStartupMessages({
  library(monocle3)
  library(SeuratWrappers)
})

# 1. Convert to CellDataSet
cds <- as.cell_data_set(epi)

# 2. Re-calculate partitions and learn graph
# Note: use_partition = FALSE can help if your tumor cells are too 'disconnected'
cds <- cluster_cells(cds, reduction = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# 3. Define root cells (Normal + Low EMT)
meta_cds <- as.data.frame(colData(cds))
emt_cut <- quantile(meta_cds$EMT_Score1, probs = 0.2, na.rm = TRUE)

root_candidates <- rownames(meta_cds)[
  meta_cds$TumorStatus == "Normal" & 
    meta_cds$EMT_Score1 <= emt_cut
]

# Fallback if EMT filter is too strict
if (length(root_candidates) < 20) {
  root_candidates <- rownames(meta_cds)[meta_cds$TumorStatus == "Normal"]
}

# 4. Order cells
cds <- order_cells(cds, root_cells = root_candidates)

# 5. FIX: Handle Infinite Pseudotime values
pseudotime_vals <- pseudotime(cds)
max_pseudo <- max(pseudotime_vals[is.finite(pseudotime_vals)], na.rm = TRUE)

# Replace Inf with the maximum finite value
pseudotime_vals[is.infinite(pseudotime_vals)] <- max_pseudo

# 6. Add back to Seurat object
epi$Pseudotime <- pseudotime_vals[colnames(epi)]

# 7. Final Plotting
p_ps <- FeaturePlot(
  epi,
  reduction = "umap",
  features = "Pseudotime",
  cols = c("navy", "gold", "red")
) + ggtitle("Epithelial Pseudotime (Corrected)")

ggsave(file.path(fig_dir, "Epi_UMAP_Pseudotime.png"), p_ps, width = 7, height = 6, dpi = 300)
###############################################################
# 6) Save updated epithelial object
###############################################################

saveRDS(epi, file.path(base_dir, "NSCLC_Epithelial_states_trajectory.rds"))
cat("\n✅ DONE: Epithelial states + pseudotime saved as:\n")
cat("   NSCLC_Epithelial_states_trajectory.rds\n")
cat("Figures in:\n   ", fig_dir, "\n")

#!/usr/bin/env Rscript

###############################################################
# 4_myeloid_subatlas.R
# Myeloid sub-atlas (monocytes, macrophages, DCs)
###############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(future)
  library(dplyr)
  library(ggplot2)
})

base_dir <- "/home/khan/nlsc_scrna"
seurat_in <- file.path(base_dir, "NSCLC_Seurat_processed.rds")
fig_dir   <- file.path(base_dir, "figures_myeloid")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(seurat_in)) {
  stop("NSCLC_Seurat_processed.rds not found.")
}

n_cores <- max(1, parallel::detectCores() - 2)
plan(multicore, workers = n_cores)
options(future.globals.maxSize = 20000 * 1024^2)

###############################################################
# 1) Load global object & subset myeloid cells
###############################################################

sce <- readRDS(seurat_in)
DefaultAssay(sce) <- "RNA"

cat("Global cells:", ncol(sce), "\n")
cat("Cell_type levels:\n")
print(table(sce$Cell_type))

my <- subset(
  sce,
  subset = Cell_type %in% c("Myeloid cells", "Myeloid", "myeloid")
)

cat("Myeloid cells:", ncol(my), "\n")

my$Sample_Origin <- factor(my$Sample_Origin)
my$Sample        <- factor(my$Sample)

###############################################################
# 2) SCTransform + PCA + UMAP + clustering
###############################################################

my <- SCTransform(my, vst.flavor = "v2", variable.features.n = 3000,
                  vars.to.regress = "percent.mt", verbose = TRUE)
DefaultAssay(my) <- "SCT"

my <- RunPCA(my, npcs = 50, verbose = FALSE)
my <- FindNeighbors(my, dims = 1:30, verbose = FALSE)
my <- FindClusters(my, resolution = 0.8, verbose = FALSE)
my <- RunUMAP(my, dims = 1:30, verbose = FALSE)

my$myeloid_cluster <- my$SCT_snn_res.0.8 %||% my$RNA_snn_res.0.8

p_umap_my <- DimPlot(my, reduction = "umap",
                     group.by = "myeloid_cluster",
                     label = TRUE, repel = TRUE) +
  ggtitle("Myeloid sub-atlas – clusters")

ggsave(file.path(fig_dir, "Myeloid_UMAP_clusters.png"),
       p_umap_my, width = 7, height = 6, dpi = 300)

p_umap_origin <- DimPlot(my, reduction = "umap",
                         group.by = "Sample_Origin") +
  ggtitle("Myeloid UMAP – Sample_Origin")

ggsave(file.path(fig_dir, "Myeloid_UMAP_by_Sample_Origin.png"),
       p_umap_origin, width = 7, height = 6, dpi = 300)

###############################################################
# 3) Myeloid marker visualization
###############################################################

marker_list <- list(
  Mono       = c("LYZ","S100A8","S100A9","VCAN"),
  Macrophage = c("CD68","LST1","MSR1","C1QA","C1QB"),
  DC         = c("LILRA4","GZMB","IRF8","CLEC9A","BATF3"),
  M2_like    = c("MRC1","CD163","MSR1"),
  M1_like    = c("IL1B","TNF","CXCL10")
)

all_markers <- unique(unlist(marker_list))

fp_my <- FeaturePlot(my, reduction = "umap",
                     features = c("LYZ","CD68","MRC1","IL1B"),
                     ncol = 2, order = TRUE)
ggsave(file.path(fig_dir, "Myeloid_UMAP_key_markers.png"),
       fp_my, width = 8, height = 8, dpi = 300)

dot_my <- DotPlot(my, features = all_markers, group.by = "myeloid_cluster") +
  RotatedAxis() +
  ggtitle("Myeloid markers by cluster")

ggsave(file.path(fig_dir, "Myeloid_DotPlot_markers_by_cluster.png"),
       dot_my, width = 9, height = 5, dpi = 300)

###############################################################
# 4) Cluster markers & basic DE (e.g. tLung vs nLung monocytes/macrophages)
###############################################################

Idents(my) <- my$myeloid_cluster

my_markers <- FindAllMarkers(my, only.pos = TRUE,
                             min.pct = 0.25, logfc.threshold = 0.25)
my_markers <- my_markers[order(my_markers$avg_log2FC, decreasing = TRUE),]

write.csv(my_markers,
          file.path(base_dir, "Myeloid_cluster_markers_res0.8.csv"),
          row.names = FALSE)

# TumorStatus (reuse same definition)
my$TumorStatus <- "Other"
my$TumorStatus[my$Sample_Origin == "nLung"] <- "Normal"
my$TumorStatus[my$Sample_Origin == "tLung"] <- "Tumor"
my$TumorStatus <- factor(my$TumorStatus, levels = c("Normal","Tumor","Other"))

# Example: DE in myeloid cells Tumor vs Normal
my_sub <- subset(my, subset = TumorStatus %in% c("Normal","Tumor"))
Idents(my_sub) <- my_sub$TumorStatus

if (all(c("Normal","Tumor") %in% levels(Idents(my_sub)))) {
  de_my <- FindMarkers(my_sub,
                       ident.1 = "Tumor",
                       ident.2 = "Normal",
                       logfc.threshold = 0.25,
                       min.pct = 0.1)
  de_my <- de_my[order(de_my$p_val_adj),]
  write.csv(de_my,
            file.path(base_dir, "Myeloid_DE_Tumor_vs_Normal.csv"),
            row.names = TRUE)
}

saveRDS(my, file.path(base_dir, "NSCLC_Myeloid.rds"))
cat("✅ Myeloid sub-atlas saved: NSCLC_Myeloid.rds\n")
cat("Figures in:", fig_dir, "\n")

#!/usr/bin/env Rscript

###############################################################
# 5_TNK_subatlas.R
# T and NK cell sub-atlas
###############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(future)
  library(dplyr)
  library(ggplot2)
})

base_dir <- "/home/khan/nlsc_scrna"
seurat_in <- file.path(base_dir, "NSCLC_Seurat_processed.rds")
fig_dir   <- file.path(base_dir, "figures_TNK")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(seurat_in)) {
  stop("NSCLC_Seurat_processed.rds not found.")
}

n_cores <- max(1, parallel::detectCores() - 2)
plan(multicore, workers = n_cores)
options(future.globals.maxSize = 20000 * 1024^2)

###############################################################
# 1) Load global object & subset T/NK cells
###############################################################

sce <- readRDS(seurat_in)
DefaultAssay(sce) <- "RNA"

cat("Global cells:", ncol(sce), "\n")
cat("Cell_type levels:\n")
print(table(sce$Cell_type))

tnk <- subset(
  sce,
  subset = Cell_type %in% c("T cells", "NK cells", "T/NK cells", "T and NK cells")
)

cat("T/NK cells:", ncol(tnk), "\n")

tnk$Sample_Origin <- factor(tnk$Sample_Origin)
tnk$Sample        <- factor(tnk$Sample)

###############################################################
# 2) SCTransform + PCA + UMAP + clustering
###############################################################

tnk <- SCTransform(tnk, vst.flavor = "v2", variable.features.n = 3000,
                   vars.to.regress = "percent.mt", verbose = TRUE)
DefaultAssay(tnk) <- "SCT"

tnk <- RunPCA(tnk, npcs = 50, verbose = FALSE)
tnk <- FindNeighbors(tnk, dims = 1:30, verbose = FALSE)
tnk <- FindClusters(tnk, resolution = 0.8, verbose = FALSE)
tnk <- RunUMAP(tnk, dims = 1:30, verbose = FALSE)

tnk$tnk_cluster <- tnk$SCT_snn_res.0.8 %||% tnk$RNA_snn_res.0.8

p_umap_tnk <- DimPlot(tnk, reduction = "umap",
                      group.by = "tnk_cluster",
                      label = TRUE, repel = TRUE) +
  ggtitle("T/NK sub-atlas – clusters")

ggsave(file.path(fig_dir, "TNK_UMAP_clusters.png"),
       p_umap_tnk, width = 7, height = 6, dpi = 300)

p_umap_origin <- DimPlot(tnk, reduction = "umap",
                         group.by = "Sample_Origin") +
  ggtitle("T/NK UMAP – Sample_Origin")

ggsave(file.path(fig_dir, "TNK_UMAP_by_Sample_Origin.png"),
       p_umap_origin, width = 7, height = 6, dpi = 300)

###############################################################
# 3) Fixed T/NK marker visualization 
###############################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Ensure the correct assay is active (SCT was used in the study for states)
DefaultAssay(tnk) <- "SCT"

# Define the marker list used to identify lineages [cite: 46, 175, 839]
marker_list <- list(
  T_core      = c("CD3D","CD3E","CD4","CD8A"),
  T_exhausted = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4"), # High in metastatic stages 
  T_act       = c("IL2RA","ICOS","TNFRSF9","CD69"),
  Cytotoxic   = c("GZMB","PRF1","NKG7","GNLY","FGFBP2"), # High in normal lung 
  NK          = c("KLRD1","KLRB1","NCR1")
)

all_markers <- unique(unlist(marker_list))

# 1. FIXED FEATURE PLOT
# Using the key markers that define the transition from normal to metastatic [cite: 842, 843]
fp_tnk <- FeaturePlot(
  tnk, 
  reduction = "umap",
  features = c("CD3D","CD8A","PDCD1","GZMB"),
  ncol = 2, 
  order = TRUE
) & 
  theme(plot.title = element_text(size = 10))

ggsave(file.path(fig_dir, "TNK_UMAP_key_markers.png"), 
       plot = fp_tnk, width = 8, height = 8, dpi = 300)

# 2. FIXED DOT PLOT
# Check if 'tnk_cluster' exists; if not, use 'seurat_clusters'
group_var <- if("tnk_cluster" %in% colnames(tnk@meta.data)) "tnk_cluster" else "seurat_clusters"

dot_tnk <- DotPlot(
  tnk, 
  features = marker_list, # Passing the list creates organized gene headers
  group.by = group_var
) + 
  RotatedAxis() + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) +
  ggtitle(paste("T/NK markers by", group_var))

# Apply Nature-style aesthetic adjustments 
dot_tnk <- dot_tnk + theme(
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8)
)

ggsave(file.path(fig_dir, "TNK_DotPlot_Markers.png"), 
       plot = dot_tnk, width = 12, height = 6, dpi = 300)
###############################################################
# 4) Fixed functional state scoring & Nature-style plotting
###############################################################
library(Seurat)
library(ggplot2)
library(patchwork)

# Set Default Assay to SCT as per the study's trajectory analysis
DefaultAssay(tnk) <- "SCT"

# Define signatures from the LUAD atlas
# Exhaustion: Markers primarily found in metastatic tissues (mLN, mBrain)
exhaustion <- list(c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4"))
# Activation: Markers for early/active immune responses
activation <- list(c("CD69","TNFRSF9","ICOS","IL2RA"))
# Cytotoxicity: Markers primarily in normal lung (nLung) effector cells
cytotoxic  <- list(c("GZMB","PRF1","NKG7","GNLY","FGFBP2"))

# Add Module Scores (Seurat appends '1' to the name)
tnk <- AddModuleScore(tnk, features = exhaustion, name = "Exhaust_Score")
tnk <- AddModuleScore(tnk, features = activation, name = "Activate_Score")
tnk <- AddModuleScore(tnk, features = cytotoxic, name = "Cytotox_Score")

score_feats <- c("Exhaust_Score1", "Activate_Score1", "Cytotox_Score1")

# Determine the grouping variable: checks for 'tnk_cluster', defaults to 'seurat_clusters'
group_var <- if("tnk_cluster" %in% colnames(tnk@meta.data)) "tnk_cluster" else "seurat_clusters"

# 1. FIXED FEATUREPLOT LOOP (Nature Figure 5c Style)
for (feat in score_feats) {
  # Use & instead of + to avoid S4SXP errors inside loops
  p <- FeaturePlot(tnk, reduction = "umap", features = feat,
                   cols = c("grey90", "darkred"), 
                   min.cutoff = "q5", max.cutoff = "q95") & 
    ggtitle(paste("T/NK", feat)) &
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
  
  ggsave(file.path(fig_dir, paste0("TNK_UMAP_", feat, ".png")),
         plot = p, width = 6, height = 5, dpi = 300)
}

# 2. FIXED VLNPLOT LOOP (Nature Figure 6b Style)
for (feat in score_feats) {
  # pt.size = 0 removes outliers for publication clarity
  p <- VlnPlot(tnk, features = feat, group.by = group_var, pt.size = 0) & 
    theme_bw(base_size = 11) & 
    ggtitle(paste(feat, "by", group_var)) &
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(file.path(fig_dir, paste0("TNK_Vln_", feat, "_by_cluster.png")),
         plot = p, width = 6.5, height = 4, dpi = 300)
}
###############################################################
# 5) Save T/NK object
###############################################################

saveRDS(tnk, file.path(base_dir, "NSCLC_TNK.rds"))
cat("✅ T/NK sub-atlas saved: NSCLC_TNK.rds\n")
cat("Figures in:", fig_dir, "\n")

###############################################################
# 5) Figures composite
###############################################################
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

base_dir <- "/home/khan/nlsc_scrna"
fig_dir <- file.path(base_dir, "final_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Replace with your real figure paths
pA <- ggdraw() + draw_image(file.path(base_dir, "figures/UMAP_by_Sample_Origin.png"))
pB <- ggdraw() + draw_image(file.path(base_dir, "figures/Composition_Sample_CellType.png"))
pC <- ggdraw() + draw_image(file.path(base_dir, "figures/Epithelial_UMAP_clusters.png"))  # if you have it
pD <- ggdraw() + draw_image(file.path(base_dir, "trajectory_myeloid/Myeloid_UMAP_clusters.png"))
pE <- ggdraw() + draw_image(file.path(base_dir, "cellchat/CellChat_pathway_heatmap.pdf"))

final <- (pA | pB) / (pC | pD) / pE +
  plot_annotation(tag_levels = "A")

ggsave(file.path(fig_dir, "Figure1_NatureStyle.png"),
       final, width = 14, height = 16, dpi = 500)

cat("✅ Saved: ", file.path(fig_dir, "Figure1_NatureStyle.png"), "\n")
