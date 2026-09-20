#Bachelor Thesis - Liv Amin 
#Single-cell RNA-seq analysis of ETV6::RUNX1 

#---------libararies-------------

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(celldex)
library(SingleCellExperiment)
library(SingleR)
library(ggrepel)
library(writexl)
library(decoupleR)

#human and mouse reference cellranger multi output 
#sample 1: preleukemic with splenomegaly 
#sample2 & 3 : preleukemic without splenomegaly 
#sample 4: control 

#-----------paths---------

base <- "/omics/odcf/analysis/OE0650_projects/rach_seq/src/l172n/"

sample1_path <- "/omics/odcf/analysis/OE0650_projects/rach_seq/src/l172n/cellranger/1295_1_2_human_mouse_multi/outs/per_sample_outs/1295-1/sample_filtered_feature_bc_matrix"
sample2_path <- "/omics/odcf/analysis/OE0650_projects/rach_seq/src/l172n/cellranger/1295_1_2_human_mouse_multi/outs/per_sample_outs/1295-2/sample_filtered_feature_bc_matrix"
sample3_path <- "/omics/odcf/analysis/OE0650_projects/rach_seq/src/l172n/cellranger/1295_3_4_human_mouse_multi/outs/per_sample_outs/1295-3/sample_filtered_feature_bc_matrix"
sample4_path <- "/omics/odcf/analysis/OE0650_projects/rach_seq/src/l172n/cellranger/1295_3_4_human_mouse_multi/outs/per_sample_outs/1295-4/sample_filtered_feature_bc_matrix"

#-----------load data-----------

# Read10X loads the sparse count matrix from cellranger output 
sample1 <- Read10X(sample1_path)
sample2 <- Read10X(sample2_path)
sample3 <- Read10X(sample3_path)
sample4 <- Read10X(sample4_path)


seurat1 <- CreateSeuratObject(counts = sample1,  project = "Sample1")
seurat2 <- CreateSeuratObject(counts = sample2,  project = "Sample2")
seurat3 <- CreateSeuratObject(counts = sample3,  project = "Sample3")
control <- CreateSeuratObject(counts = sample4,  project = "Control")

seurat1$sample_group <- "Sample1"
seurat2$sample_group <- "Sample2"
seurat3$sample_group <- "Sample3"
control$sample_group <- "Control"


#----------QC and filtering---------------

#mitochondrial fraction, restricted to the murine part of the reference
seurat1[["percent.mt"]] <- PercentageFeatureSet(seurat1, pattern = "^GRCm39-mt-")
seurat2[["percent.mt"]] <- PercentageFeatureSet(seurat2, pattern = "^GRCm39-mt-")
seurat3[["percent.mt"]] <- PercentageFeatureSet(seurat3, pattern = "^GRCm39-mt-")
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^GRCm39-mt-")

#QC metrics per sample before filtering
VlnPlot(seurat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(seurat2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(seurat3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seurat1, "nCount_RNA","percent.mt") + FeatureScatter(seurat1,"nCount_RNA","nFeature_RNA")
FeatureScatter(seurat2, "nCount_RNA","percent.mt") + FeatureScatter(seurat2,"nCount_RNA","nFeature_RNA")
FeatureScatter(seurat3, "nCount_RNA","percent.mt") + FeatureScatter(seurat3,"nCount_RNA","nFeature_RNA")
FeatureScatter(control, "nCount_RNA","percent.mt") + FeatureScatter(control,"nCount_RNA","nFeature_RNA")

#filter cells 

seurat1 <- subset(seurat1, subset = nFeature_RNA > 1500 & percent.mt < 12.5)
seurat2 <- subset(seurat2, subset = nFeature_RNA > 1500 & percent.mt < 12.5)
seurat3 <- subset(seurat3, subset = nFeature_RNA > 1500 & percent.mt < 12.5)
control <- subset(control, subset = nFeature_RNA > 1500 & percent.mt < 12.5)


#----------preprocessing and integration-------- 

#merging all samples and normalization 

merged_obj <- merge(seurat1, y = list(seurat2, seurat3, control), add.cell.ids = c("S1", "S2", "S3", "C"), project = "AllSamples")
merged_obj <- NormalizeData(merged_obj)

#variable features per layer
#FindVariableFeatures identifies genes that show high variability across cells
merged_obj <- FindVariableFeatures(merged_obj)
VariableFeaturePlot(merged_obj)

#scaling data -> Put all genes on the same scale so they are comparable
#ScaleData prepares for PCA and clustering 
merged_obj <- ScaleData(merged_obj)

#----------dimensionality reduction and unintegrated clustering------

#pca 
merged_obj <- RunPCA(merged_obj, reduction.name = "pca")
DimPlot(merged_obj, reduction = "pca")
ElbowPlot(merged_obj)

merged_obj <- FindNeighbors(merged_obj, dims = 1:20, reduction = "pca")
merged_obj <- FindClusters(merged_obj, resolution = 0.5, cluster.name = "unintegrated_clusters")

#UMAP (unintegrated)
merged_obj <- RunUMAP(merged_obj, dims = 1:20, reduction = "pca", reduction.name = "umap.unintegrated")

merged_obj_umap_plot <- DimPlot(
  merged_obj,
  reduction = "umap.unintegrated",
  group.by = c("orig.ident", "unintegrated_clusters"),
  combine = FALSE, label.size = 2
)
merged_obj_umap_plot

# CCA integration 
obj <- IntegrateLayers(
  object = merged_obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5, cluster.name = "cca_clusters")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:20, reduction.name = "umap.cca")

integerated_obj_umap_plot <- DimPlot(
  obj,
  reduction = "umap.cca",
  group.by = c("orig.ident", "cca_clusters"),
  combine = FALSE, label.size = 2
)
integerated_obj_umap_plot

# compare before and after integration 
merged_obj_umap_plot[[1]] | integerated_obj_umap_plot[[1]]


#------cell type annotation (SingleR)---------

#ImmGenData -> a mouse immune cell reference 
ref_imm <- celldex::ImmGenData()
colnames(colData(ref_imm))
unique(ref_imm$label.main)

#rejoin layers
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])

#keep only mouse genes
mouse_genes <- grep("^GRCm39-", rownames(obj), value = TRUE)
obj_mouse <- obj[mouse_genes, ]

# emove the GRCm39- prefix so names match ImmGen symbols
rownames(obj_mouse) <- sub("^GRCm39-", "", rownames(obj_mouse))

#check
head(rownames(obj_mouse), 20)
obj_mouse[["RNA"]] <- JoinLayers(obj_mouse[["RNA"]])

#run SingleR for cell type annotation 
sce_mouse <- as.SingleCellExperiment(obj_mouse)

singleR_results <- SingleR(
  test = sce_mouse,
  ref = ref_imm,
  labels = ref_imm$label.main
)

#store the predicted labels from SingleR and stores them as metadata in the Seurat objects
obj_mouse$SingleR_label <- singleR_results$labels
# Transfer SingleR labels to obj
obj$SingleR_label <- singleR_results$labels

DimPlot(
  obj_mouse,
  reduction = "umap.cca",
  group.by = "SingleR_label",
  label = TRUE,
  repel = TRUE
) + NoLegend()

DimPlot(obj,
        reduction = "umap.cca",
        group.by  = "SingleR_label",
        label     = TRUE,
        repel     = TRUE,
        label.size = 5,        # make labels bigger 
        pt.size   = 0.5) +
  NoLegend() +
  ggtitle("UMAP with SingleR labels") +
  theme(plot.title = element_text(face = "bold"),text = element_text(face = "bold"))

#------fusion gene detection---------- 

#barcodes of cells with atleast two reads over RUNX1 exons 5-9
barcode_base <- paste0(base, "txt_files/")

runx1_barcodes <- list(
  Sample1 = readLines(paste0(barcode_base, "RUNX1_exons5to9_Sample1_barcodes_min2.txt")),
  Sample2 = readLines(paste0(barcode_base, "RUNX1_exons5to9_Sample2_barcodes_min2.txt")),
  Sample3 = readLines(paste0(barcode_base, "RUNX1_exons5to9_Sample3_barcodes_min2.txt")),
  Control = readLines(paste0(barcode_base, "RUNX1_exons5to9_Control_barcodes_min2.txt"))
)

# match barcodes to obj's actual cell names 
prefix_map <- c(Sample1 = "S1_", Sample2 = "S2_", Sample3 = "S3_", Control = "C_")

runx1_cells <- unlist(lapply(names(runx1_barcodes), function(s) {
  paste0(prefix_map[s], runx1_barcodes[[s]], "-1")
}))

sum(runx1_cells %in% colnames(obj))
length(runx1_cells)

#-----RUNX1 positive cells per cell type and sample-----

obj_meta <- data.frame(
  cell = colnames(obj),
  sample_group = obj$sample_group,
  SingleR_label = obj$SingleR_label
)
obj_meta$RUNX1_pos <- obj_meta$cell %in% runx1_cells

# ---- barplot ----
obj_meta %>%
  filter(sample_group %in% c("Sample1", "Sample2", "Sample3", "Control")) %>%
  group_by(sample_group, SingleR_label) %>%
  summarise(
    n_RUNX1_pos = sum(RUNX1_pos),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = SingleR_label, y = n_RUNX1_pos, fill = sample_group)) +
    geom_col(position = "dodge") +
    RotatedAxis() +
    scale_fill_manual(values = c("Sample1" = "#0d59dd", 
                                  "Sample2" = "#42ad7d", 
                                  "Sample3" = "#b15cc0",
                                  "Control" = "#f40505")) +
    labs(title = "RUNX1 exon5-9+ cells: Sample1, Sample2, Sample3, Control",
         x = "Cell type", y = "Number of RUNX1 exon5-9+ cells") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
      axis.text.y = element_text(size = 16),
      axis.title  = element_text(size = 18),
      plot.title  = element_text(size = 19),
      legend.title = element_text(size = 17),
      legend.text  = element_text(size = 16))


#number of RUNX1 positive cells within each cell type per sample
runx1_counts <- obj_meta %>%
  group_by(SingleR_label, sample_group) %>%
  summarise(
    cells    = n(),
    positive = sum(RUNX1_pos),
    .groups  = "drop"
  ) %>%
  filter(positive > 0) %>%
  arrange(SingleR_label, sample_group)

print(runx1_counts, n = Inf)

##----------------DE analysis----------------

# Samples 2 and 3 are the samples without splenomegaly and are pooled into one group
obj$de_group <- as.character(obj$sample_group)
obj$de_group[obj$de_group %in% c("Sample2", "Sample3")] <- "S23"

# test murine features only, so that human transgene features cannot
# come out as differential between transgenic and control
mouse_genes <- grep("^GRCm39-", rownames(obj), value = TRUE)


#----------DE: Samples 2 and 3 vs control----------

sub <- subset(obj, subset = SingleR_label == "Stem cells" & de_group %in% c("S23", "Control"))
Idents(sub) <- "de_group"
de_stem_S23vsCtrl <- FindMarkers(sub, ident.1 = "S23", ident.2 = "Control",
                                 features = mouse_genes, test.use = "wilcox",
                                 min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj, subset = SingleR_label == "Neutrophils" & de_group %in% c("S23", "Control"))
Idents(sub) <- "de_group"
de_neut_S23vsCtrl <- FindMarkers(sub, ident.1 = "S23", ident.2 = "Control",
                                 features = mouse_genes, test.use = "wilcox",
                                 min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj, subset = SingleR_label == "Monocytes" & de_group %in% c("S23", "Control"))
Idents(sub) <- "de_group"
de_mono_S23vsCtrl <- FindMarkers(sub, ident.1 = "S23", ident.2 = "Control",
                                 features = mouse_genes, test.use = "wilcox",
                                 min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj, subset = SingleR_label == "B cells, pro" & de_group %in% c("S23", "Control"))
Idents(sub) <- "de_group"
de_bpro_S23vsCtrl <- FindMarkers(sub, ident.1 = "S23", ident.2 = "Control",
                                 features = mouse_genes, test.use = "wilcox",
                                 min.pct = 0.1, logfc.threshold = 0)


#----------DE: Sample 1 vs Samples 2 and 3----------

sub <- subset(obj, subset = SingleR_label == "Stem cells" & de_group %in% c("Sample1", "S23"))
Idents(sub) <- "de_group"
de_stem_S1vsS23 <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "S23",
                               features = mouse_genes, test.use = "wilcox",
                               min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj, subset = SingleR_label == "Neutrophils" & de_group %in% c("Sample1", "S23"))
Idents(sub) <- "de_group"
de_neut_S1vsS23 <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "S23",
                               features = mouse_genes, test.use = "wilcox",
                               min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj, subset = SingleR_label == "Monocytes" & de_group %in% c("Sample1", "S23"))
Idents(sub) <- "de_group"
de_mono_S1vsS23 <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "S23",
                               features = mouse_genes, test.use = "wilcox",
                               min.pct = 0.1, logfc.threshold = 0)


#----------DE: Sample 1 vs control----------

sub <- subset(obj, subset = SingleR_label == "Stem cells" & de_group %in% c("Sample1", "Control"))
Idents(sub) <- "de_group"
de_stem_S1vsCtrl <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "Control",
                                features = mouse_genes, test.use = "wilcox",
                                min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj, subset = SingleR_label == "Neutrophils" & de_group %in% c("Sample1", "Control"))
Idents(sub) <- "de_group"
de_neut_S1vsCtrl <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "Control",
                                features = mouse_genes, test.use = "wilcox",
                                min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj, subset = SingleR_label == "Monocytes" & de_group %in% c("Sample1", "Control"))
Idents(sub) <- "de_group"
de_mono_S1vsCtrl <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "Control",
                                features = mouse_genes, test.use = "wilcox",
                                min.pct = 0.1, logfc.threshold = 0)


#----------volcano plot----------

figures_dir <- paste0(base, "thesis_figures/")
dir.create(figures_dir, showWarnings = FALSE)

make_volcano <- function(df, title,
                         label_n    = 15,
                         label_size = 3.2,
                         base_size  = 10) {

  df$names          <- sub("^GRCm39-", "", rownames(df))
  df$logfoldchanges <- df$avg_log2FC
  df$pvals_adj      <- df$p_val_adj

  df$sig_type <- "Not significant"
  df$sig_type[df$pvals_adj < 0.05 & df$logfoldchanges >  0.58] <- "Upregulated"
  df$sig_type[df$pvals_adj < 0.05 & df$logfoldchanges < -0.58] <- "Downregulated"

  # top genes per direction, labelled by adjusted p-value
  top <- df[df$sig_type == "Upregulated", ]
  top <- head(top[order(top$pvals_adj), ], label_n)
  bot <- df[df$sig_type == "Downregulated", ]
  bot <- head(bot[order(bot$pvals_adj), ], label_n)
  to_label <- rbind(top, bot)

  # Ly6a labelled as well, if significant
  ly6a <- df[df$names == "Ly6a" & df$sig_type != "Not significant", ]
  if (nrow(ly6a) > 0) {
    to_label <- rbind(to_label, ly6a)
    to_label <- to_label[!duplicated(to_label$names), ]
  }

  to_label$label_color <- ifelse(to_label$sig_type == "Downregulated",
                                 "#1a0294", "#a0047c")


  ggplot(df, aes(logfoldchanges, -log10(pvals_adj + 1e-300), color = sig_type)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = c("Downregulated"   = "#7156f6",
                                  "Not significant" = "#CCCCCC",
                                  "Upregulated"     = "#f556d0")) +
    geom_text_repel(data = to_label, aes(label = names),
                    color = to_label$label_color,
                    size = label_size, fontface = "bold",
                    max.overlaps = 30, box.padding = 0.3,
                    segment.color = NA, show.legend = FALSE) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
    coord_cartesian(xlim = c(-5, 5)) +
    labs(title = title, x = "Log2 Fold Change",
         y = "-log10 adj. p-value", color = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "right")
}


#volcano plots and significant gene counts

de_list <- list(
  "Stem cells_S23vsControl"  = de_stem_S23vsCtrl,
  "Neutrophils_S23vsControl" = de_neut_S23vsCtrl,
  "Monocytes_S23vsControl"   = de_mono_S23vsCtrl,
  "B cells pro_S23vsControl" = de_bpro_S23vsCtrl,
  "Stem cells_S1vsS23"       = de_stem_S1vsS23,
  "Neutrophils_S1vsS23"      = de_neut_S1vsS23,
  "Monocytes_S1vsS23"        = de_mono_S1vsS23,
  "Stem cells_S1vsControl"   = de_stem_S1vsCtrl,
  "Neutrophils_S1vsControl"  = de_neut_S1vsCtrl,
  "Monocytes_S1vsControl"    = de_mono_S1vsCtrl
)

de_counts <- list()

for (title in names(de_list)) {

  df <- de_list[[title]]
  tag <- gsub(" ", "_", title)

  p <- make_volcano(df, title)
  ggsave(paste0(figures_dir, "DE_", tag, ".png"),
         plot = p, dpi = 300, width = 8.5, height = 5.5, units = "in")
  print(p)

  df$names          <- sub("^GRCm39-", "", rownames(df))
  df$logfoldchanges <- df$avg_log2FC
  df$pvals_adj      <- df$p_val_adj
  df$sig_type <- "Not significant"
  df$sig_type[df$pvals_adj < 0.05 & df$logfoldchanges >  0.58] <- "Upregulated"
  df$sig_type[df$pvals_adj < 0.05 & df$logfoldchanges < -0.58] <- "Downregulated"

  write.csv(df, paste0(figures_dir, "DE_", tag, ".csv"))

  de_counts[[title]] <- data.frame(
    comparison   = title,
    genes_tested = nrow(df),
    up           = sum(df$sig_type == "Upregulated"),
    down         = sum(df$sig_type == "Downregulated")
  )
}

de_counts <- do.call(rbind, de_counts)
de_counts$total_sig <- de_counts$up + de_counts$down
print(de_counts)
write.csv(de_counts, paste0(figures_dir, "DE_summary_counts.csv"), row.names = FALSE)

##save and reload the processed objects
saveRDS(obj,       paste0(base, "obj_processed.rds"))
saveRDS(obj_mouse, paste0(base, "obj_mouse_processed.rds"))

#---------------COLLECTRI--------------------

# obj       <- readRDS(paste0(base, "obj_processed.rds"))
# obj_mouse <- readRDS(paste0(base, "obj_mouse_processed.rds"))

#paths and regulon 

output_dir  <- paste0(base, "collectri_output/")
scores_path <- paste0(output_dir, "tf_activity_all.rds")
dir.create(output_dir, showWarnings = FALSE)

cell_types <- c("Stem cells", "Neutrophils", "Monocytes")

# collecTRI mouse regulon: for each TF, its target genes and the mode of regulation
collecTRI <- get_collectri(organism = "mouse", split_complexes = FALSE)


#ULM 
# he univariate linear model scores, per cell and per TF, whether the targets of that TF are expressed above or below the background of that cell.
#This step is computationally heavy so results are saved and reused.

if (!file.exists(scores_path)) {

  # normalized expression matrix (RNA assay, log-normalized layer)
  mat <- GetAssayData(obj_mouse, assay = "RNA", layer = "data")

  activity_scores <- run_ulm(mat = mat, network = collecTRI)

  # reshape from long to wide format: rows = TFs, columns = cells
  activity_scores_wide <- activity_scores %>%
    filter(statistic == "ulm") %>%
    select(source, condition, score) %>%
    pivot_wider(names_from = condition, values_from = score) %>%
    column_to_rownames("source") %>%
    as.matrix()

  saveRDS(activity_scores_wide, scores_path)
  rm(mat, activity_scores)
  gc()
}

activity_scores_wide <- readRDS(scores_path)


#attach as assays 

#TFact holds the activity scores, which lets Seurat functions (FindMarkers, FeaturePlot, etc.) be applied to them exactly as to gene expression
obj_mouse[["TFact"]] <- CreateAssayObject(counts = activity_scores_wide)
obj_mouse            <- ScaleData(obj_mouse, assay = "TFact")

#TFexpr holds the expression of the TFs themselves, so that inferred activity and measured expression can be compared directly
tf_genes_present      <- rownames(activity_scores_wide)[rownames(activity_scores_wide) %in% rownames(obj_mouse)]
mat_TF                <- GetAssayData(obj_mouse, assay = "RNA", layer = "data")[tf_genes_present, ]
obj_mouse[["TFexpr"]] <- CreateAssayObject(counts = mat_TF)
obj_mouse             <- ScaleData(obj_mouse, assay = "TFexpr")
rm(mat_TF)

#same grouping as the DE analysis
obj_mouse$de_group <- as.character(obj_mouse$sample_group)
obj_mouse$de_group[obj_mouse$de_group %in% c("Sample2", "Sample3")] <- "S23"


#TF activity: Samples 2 and 3 vs control 

sub <- subset(obj_mouse, subset = SingleR_label == "Stem cells" & de_group %in% c("S23", "Control"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_stem_S23vsCtrl <- FindMarkers(sub, ident.1 = "S23", ident.2 = "Control",
                                 test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj_mouse, subset = SingleR_label == "Neutrophils" & de_group %in% c("S23", "Control"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_neut_S23vsCtrl <- FindMarkers(sub, ident.1 = "S23", ident.2 = "Control",
                                 test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj_mouse, subset = SingleR_label == "Monocytes" & de_group %in% c("S23", "Control"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_mono_S23vsCtrl <- FindMarkers(sub, ident.1 = "S23", ident.2 = "Control",
                                 test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj_mouse, subset = SingleR_label == "B cells, pro" & de_group %in% c("S23", "Control"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_bpro_S23vsCtrl <- FindMarkers(sub, ident.1 = "S23", ident.2 = "Control",
                                 test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)


#TF activity: Sample 1 vs Samples 2 and 3 

sub <- subset(obj_mouse, subset = SingleR_label == "Stem cells" & de_group %in% c("Sample1", "S23"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_stem_S1vsS23 <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "S23",
                               test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj_mouse, subset = SingleR_label == "Neutrophils" & de_group %in% c("Sample1", "S23"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_neut_S1vsS23 <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "S23",
                               test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj_mouse, subset = SingleR_label == "Monocytes" & de_group %in% c("Sample1", "S23"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_mono_S1vsS23 <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "S23",
                               test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)


#TF activity: Sample 1 vs control 

sub <- subset(obj_mouse, subset = SingleR_label == "Stem cells" & de_group %in% c("Sample1", "Control"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_stem_S1vsCtrl <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "Control",
                                test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj_mouse, subset = SingleR_label == "Neutrophils" & de_group %in% c("Sample1", "Control"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_neut_S1vsCtrl <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "Control",
                                test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)

sub <- subset(obj_mouse, subset = SingleR_label == "Monocytes" & de_group %in% c("Sample1", "Control"))
DefaultAssay(sub) <- "TFact"; Idents(sub) <- "de_group"
tf_mono_S1vsCtrl <- FindMarkers(sub, ident.1 = "Sample1", ident.2 = "Control",
                                test.use = "wilcox", min.pct = 0.1, logfc.threshold = 0)

rm(sub)


#TF volcano plot function 

make_tf_volcano <- function(df, title,
                            label_n    = 20,
                            label_size = 3.2,
                            base_size  = 11,
                            highlight  = character(0)) {

  df$TF <- rownames(df)

  df$sig_type <- "Not significant"
  df$sig_type[df$p_val_adj < 0.05 & df$avg_log2FC >  0.58] <- "Upregulated"
  df$sig_type[df$p_val_adj < 0.05 & df$avg_log2FC < -0.58] <- "Downregulated"

  # top TFs per direction, labelled by adjusted p-value
  top <- df[df$sig_type == "Upregulated", ]
  top <- head(top[order(top$p_val_adj), ], label_n)
  bot <- df[df$sig_type == "Downregulated", ]
  bot <- head(bot[order(bot$p_val_adj), ], label_n)
  to_label <- rbind(top, bot)

  # TFs labelled regardless of rank, e.g. Pax5
  if (length(highlight) > 0) {
    to_label <- rbind(to_label, df[df$TF %in% highlight, ])
    to_label <- to_label[!duplicated(to_label$TF), ]
  }

  to_label$label_color <- ifelse(to_label$sig_type == "Downregulated", "#057c0b",
                          ifelse(to_label$sig_type == "Upregulated",   "#9004a0",
                                 "#666666"))

  # 1e-300 keeps p-values that underflow to zero from becoming Inf
  ggplot(df, aes(avg_log2FC, -log10(p_val_adj + 1e-300), color = sig_type)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = c("Downregulated"   = "#0ad132",
                                  "Not significant" = "#CCCCCC",
                                  "Upregulated"     = "#b25fdc")) +
    geom_text_repel(data = to_label, aes(label = TF),
                    color = to_label$label_color,
                    size = label_size, fontface = "bold",
                    max.overlaps = 30, box.padding = 0.3,
                    segment.color = NA, show.legend = FALSE) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
    coord_cartesian(xlim = c(-5, 5)) +
    labs(title = title, x = "Log2 Fold Change",
         y = "-log10 adj. p-value", color = "") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "right")
}


#TF volcano plots and significant TF counts 

tf_list <- list(
  "Stem cells_S23vsControl"  = tf_stem_S23vsCtrl,
  "Neutrophils_S23vsControl" = tf_neut_S23vsCtrl,
  "Monocytes_S23vsControl"   = tf_mono_S23vsCtrl,
  "B cells pro_S23vsControl" = tf_bpro_S23vsCtrl,
  "Stem cells_S1vsS23"       = tf_stem_S1vsS23,
  "Neutrophils_S1vsS23"      = tf_neut_S1vsS23,
  "Monocytes_S1vsS23"        = tf_mono_S1vsS23,
  "Stem cells_S1vsControl"   = tf_stem_S1vsCtrl,
  "Neutrophils_S1vsControl"  = tf_neut_S1vsCtrl,
  "Monocytes_S1vsControl"    = tf_mono_S1vsCtrl
)

tf_counts <- list()

for (title in names(tf_list)) {

  df  <- tf_list[[title]]
  tag <- gsub(" ", "_", title)

  p <- make_tf_volcano(df, paste("TF activity -", title), highlight = "Pax5")
  ggsave(paste0(figures_dir, "TF_", tag, ".png"),
         plot = p, dpi = 300, width = 8, height = 6, units = "in")
  print(p)

  df$TF <- rownames(df)
  df$sig_type <- "Not significant"
  df$sig_type[df$p_val_adj < 0.05 & df$avg_log2FC >  0.58] <- "Upregulated"
  df$sig_type[df$p_val_adj < 0.05 & df$avg_log2FC < -0.58] <- "Downregulated"

  write.csv(df, paste0(figures_dir, "TF_", tag, ".csv"), row.names = FALSE)

  tf_counts[[title]] <- data.frame(
    comparison = title,
    tfs_tested = nrow(df),
    up         = sum(df$sig_type == "Upregulated"),
    down       = sum(df$sig_type == "Downregulated")
  )
}

tf_counts <- do.call(rbind, tf_counts)
tf_counts$total_sig <- tf_counts$up + tf_counts$down
print(tf_counts)
write.csv(tf_counts, paste0(figures_dir, "TF_summary_counts.csv"), row.names = FALSE)

