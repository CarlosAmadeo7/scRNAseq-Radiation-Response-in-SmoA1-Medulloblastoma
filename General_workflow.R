#' ---
#' title: "sc-RNA seq analysis on medulloblastoma samples exposed to radiation (5h, 24h) and Control"
#' author: "Carlos Alfaro"
#' date: "2025-07-02"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ----Libraries-----------------------------------------------------------------------------------------------------
library(SingleCellExperiment);library(Seurat);library(tidyverse);library(Matrix);library(scales);
library(cowplot);library(RCurl);library(openxlsx);library(knitr);library(monocle3);library(SeuratWrappers);library(hdf5r);library(readxl);library(scToppR);library(clustree);library(dittoSeq);library(presto);library(SeuratExtend); library(harmony)

#' 
## ----Loading-------------------------------------------------------------------------------------------------------
#### Loading dataset for the paper
file_paths <- list(
    "Sham_rep1" = "data/Sham/ReplicateA1/",
    "Sham_rep2" = "data/Sham/ReplicateA5/",
    "Sham_rep3" = "data/Sham/ReplicateB2/",
    "5h_rep1" = "data/5h/ReplicateA4/",
    "5h_rep2" = "data/5h/ReplicateB2/",
    "5h_rep3" = "data/5h/ReplicateC1/",
    "5h_rep4" = "data/5h/ReplicateC2/",
    "24h_rep1" = "data/24h/ReplicateA3/",
    "24h_rep2" = "data/24h/ReplicateB1/",
    "24h_rep3" = "data/24h/ReplicateB3/"
    )

raw_list<-list()
cat("Reading raw matrices...\n")
for (name in names(file_paths)) {
  raw_list[[name]] <- Read10X(data.dir = file_paths[[name]])}
cat("done.\n")
## Create a Seurat object 
min.cells <- 3
min.features <- 200
seurat_list<-list()
cat("Creating a Seurat object...")
for (name in names(file_paths)) {
  df <- raw_list[[name]]   
  seu <- CreateSeuratObject(counts = df, project = name,min.cells = min.cells,min.features = min.features)
  seurat_list[[name]] <- seu}
cat("done...!!!\n")

cat("Merging all samples into a single one Seurat...")
GT<-merge(seurat_list[[1]], y = seurat_list[2:10], add.cell.ids = names(seurat_list),project = "Mice_sC")
cat("done!...")

#' 
## ----Adding metadata to the merged object--------------------------------------------------------------------------
metadata_mapping <- data.frame(
  SampleID = c("Sham_rep1", "Sham_rep2", "Sham_rep3","5h_rep1","5h_rep2","5h_rep3","5h_rep4",
               "24h_rep1","24h_rep2","24h_rep3"),
  Treatment = c("Sham", "Sham", "Sham","1Gy","1Gy","1Gy","1Gy","1Gy","1Gy","1Gy"),
  Time = c("5h", "5h", "5h","5h","5h","5h","5h","24h","24h","24h"),
  Status = c("Control", "Control", "Control","5h", "5h", "5h", "5h","24h", "24h", "24h"))

# Sample name modification 
cell_sample_ids <- sapply(strsplit(colnames(GT), "_"), function(x) paste(x[-length(x)], collapse = "_"))
#Data frame for all the info
cell_metadata <- data.frame(
  SampleID = metadata_mapping$SampleID[match(cell_sample_ids, metadata_mapping$SampleID)],
  Treatment = metadata_mapping$Treatment[match(cell_sample_ids, metadata_mapping$SampleID)],
  Time = metadata_mapping$Time[match(cell_sample_ids, metadata_mapping$SampleID)],
  Status = metadata_mapping$Status[match(cell_sample_ids, metadata_mapping$SampleID)]
)
# Modify rownames
rownames(cell_metadata) <- colnames(GT)
# Merge metadata to original one 
GT <- AddMetaData(GT, metadata = cell_metadata)

#----------------------
# p.Mito, p.Ribo, etc.
ribosomal_genes <- grep("^Rp[sl]", rownames(GT), value = TRUE, ignore.case = TRUE)
#percent.mt
GT[['pMito']] <- PercentageFeatureSet(GT, pattern = '^mt-')
if(sum(GT[['pMito']], na.rm = TRUE) == 0) {GT[['pMito']] <- PercentageFeatureSet(GT, pattern = '^MT-')}
if(sum(GT[['pMito']], na.rm = TRUE) == 0) {GT[['pMito']] <- PercentageFeatureSet(GT, pattern = '^Mt-')}
#p.ribo
GT[["pRibo"]] <- PercentageFeatureSet(GT, pattern = "^Rp[sl]")

# Novelty score
GT@meta.data<-GT@meta.data |>dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA) |>dplyr::mutate(log10GenesPerUMI = log10(nUMI) / log10(nGene))
### Adding one column for sample
#GT@meta.data$Sample<- sapply(strsplit(colnames(GT), "_"), function(x) paste(x[-length(x)], collapse = "_"))
View(GT@meta.data)

# Save unfiltered seurat object 
cat("Saving unfiltered object...")
saveRDS(GT, file = "output/Mouse_combined_unfiltered_GT.rds")
cat('done!\n')
# Violin plot
dir.create("output/QC")
pdf("output/QC/QC_control_features_before_SampleID.pdf", width = 15, height = 8)
VlnPlot(GT, features = c("nUMI", "nGene", "pMito", "pRibo", "log10GenesPerUMI"), pt.size = 0,group.by = 'SampleID', ncol = 4) + theme(legend.position = "none")
invisible(dev.off())

#' 
## ----QC Control individual plots-----------------------------------------------------------------------------------
pdf("output/QC/Number_cells.pdf", width = 8, height = 7)
GT@meta.data %>% ggplot(aes(x=Status, fill=Status)) + geom_bar() +theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +theme(axis.title.x = element_blank())+
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  theme(plot.caption = element_text(size = 14, hjust = 0.5, face = "bold")) +ggtitle("Number of Cells") +
  labs(caption = "Number of cells/sample")
invisible(dev.off())

# Viz of the UMIs per seurat (nCount_RNA)
pdf("output/QC/UMI_per_sample.pdf", width = 8, height = 7)
GT@meta.data %>% ggplot(aes(color=SampleID, x=nUMI, fill= SampleID)) + geom_density(alpha = 0.2) + 
  scale_x_log10() + theme_classic() +theme(plot.caption = element_text(size = 14, hjust = 0.5, face = "bold")) +
  ylab("Cell Density") +geom_vline(xintercept = c(300,50000))+labs(caption = "nUMIs/sample")
invisible(dev.off())

# Viz the distribution of genes detected per cell via histogram (n_FeatureRNA)
pdf("output/QC/Genes_per_sample.pdf", width = 8, height = 7)
GT@meta.data %>% ggplot(aes(color=Status, x=nGene, fill= Status)) + geom_density(alpha = 0.2) + theme_classic() +
  theme(plot.caption = element_text(size = 14, hjust = 0.5, face = "bold")) +scale_x_log10() + 
  geom_vline(xintercept = 250)+labs(caption = "nGene/sample")
invisible(dev.off())

# Viz the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
pdf("output/QC/Novelty_score.pdf", width = 8, height = 7)
GT@meta.data %>%ggplot(aes(x=log10GenesPerUMI, color = Status, fill=Status)) +
  geom_density(alpha = 0.2) +theme_classic() +theme(plot.caption = element_text(size = 14, hjust = 0.5, face = "bold")) +
  geom_vline(xintercept = 0.70)+labs(caption = "log10GenesPerUMI/sample")
invisible(dev.off())

# Viz the percent of mito ratio
pdf("output/QC/p.Mito.pdf", width = 8, height = 7)
GT@meta.data %>% ggplot(aes(color=SampleID, x=pMito, fill=Status)) + geom_density(alpha = 0.2) + 
  scale_x_log10() + theme_classic() +theme(plot.caption = element_text(size = 14, hjust = 0.5, face = "bold"))+
  geom_vline(xintercept = 25)+labs(caption = "p.Mito/sample")
invisible(dev.off())

# Viz of p Ribo
pdf("output/QC/p.Ribo.pdf", width = 8, height = 7)
GT@meta.data %>% ggplot(aes(color=Status, x=pRibo, fill=Status)) + geom_density(alpha = 0.2) + 
  scale_x_log10() + theme_classic() +theme(plot.caption = element_text(size = 14, hjust = 0.5, face = "bold"))+labs(caption = "p.Ribo/sample")
invisible(dev.off())

# Viz the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf("output/QC/Umi_per_gene_log10.pdf", width = 8, height = 7)
GT@meta.data %>% ggplot(aes(x=nUMI, y=nGene, color=pMito)) + geom_point() + 
  scale_colour_gradient(low = "green", high = "black", name = "p.Mito") +stat_smooth(method = lm, se = FALSE, color = "blue") + scale_x_log10() + scale_y_log10() + theme_classic() + 
  theme(plot.caption = element_text(size = 14, hjust = 0.5, face = "bold"))+geom_vline(xintercept = 0.1) +
  geom_vline(xintercept = c(250,50000), linetype = "dashed", color = "red") + 
  geom_hline(yintercept = c(250,10000), linetype = "dashed", color = "red") + 
  facet_wrap(~Status) + labs(title = "UMI vs Gene Count", x = "Log10 of nUMI", y = "Log10 of nGene")+
  labs(caption = "Quality control metrics/sample")
invisible(dev.off())

#' 
#' 
## ----Filtering based on stats--------------------------------------------------------------------------------------
quantile(GT$nUMI,  probs = c(0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999), na.rm = TRUE)
quantile(GT$nGene, probs = c(0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999), na.rm = TRUE)
quantile(GT$pMito, probs = c(0.50, 0.80, 0.90, 0.95, 0.99), na.rm = TRUE)
quantile(GT$log10GenesPerUMI, probs = c(0.01, 0.05, 0.1, 0.5), na.rm = TRUE)

# Subset
GT_filtered<- subset(GT,nUMI > 300 & nUMI < 50000 &nGene > 250 &pMito <= 25 &log10GenesPerUMI >= 0.85)

cat("Saving the new GT_filtered object...")
saveRDS(GT_filtered, file = "output/Mice_combined_filtered_GT.rds")
cat("done...")

pdf("output/QC/QC_after_control_features_SampleID.pdf", width = 15, height = 8)
VlnPlot(GT_filtered, features = c("nUMI", "nGene", "pMito", "pRibo"), pt.size = 0,group.by = 'SampleID', ncol = 4) + theme(legend.position = "none")
invisible(dev.off())

GT_filtered@meta.data[GT_filtered@meta.data$SampleID == "5h_rep1",]

## ----Removing mito genes from downstream analysis------------------------------------------------------------------
GT_filtered<- readRDS("output/Mice_combined_filtered_GT.rds")

gene.list <- rownames(GT_filtered)
new.gene.list <- gsub('^mt-(.*)', NA, gene.list)
final.gene.list <- new.gene.list[!is.na(new.gene.list)]
# Removing the mitochondrial genes 
GT_filtered_1 <- subset(GT_filtered, features = final.gene.list)
# Checking the mitocondrial genes
grep("^mt-", rownames(GT_filtered_1), value = T, ignore.case = T)
saveRDS(GT_filtered_1, file = "output/Mouse_combined_filtered_no_mito_genes_GT.rds")

#' 
## ----SCTransform workflow------------------------------------------------------------------------------------------
# Join layers 
GT_filtered_1<-JoinLayers(GT_filtered_1)
## Split by sample 
cat("Splitting the object...")
seuObject_split <- SplitObject(GT_filtered_1, split.by = "SampleID")
cat("done happily...")

# Running SC Transform and regressing out the pmito, numi and pRibo
options(future.globals.maxSize = 8000 * 1024^2)  # 2GB
cat("SC_transfomr is processing now...")
for (i in 1:length(seuObject_split)) {
  message("Running SC Transform on : ", names(seuObject_split)[i])
  seuObject_split[[i]] <- SCTransform(seuObject_split[[i]],vars.to.regress = c("nUMI", "pMito", "pRibo"),verbose = FALSE) 
  gc()}
cat("done SC transform! ready for integration now!...")

#' 
## ----Calculating Phases score--------------------------------------------------------------------------------------
load("MouseCellCycleGenes.rda") 

seuObject_split <- lapply(seuObject_split, function(x) {
  x<-CellCycleScoring(x,s.features = s_genes,g2m.features = g2m_genes,set.ident = TRUE)
  return(x)})

# Create the Difference on Cell cycling stuff for each object in the dataset
seuObject_split <- lapply(seuObject_split, function(x) {x$CC.Difference <- x$S.Score - x$G2M.Score
  return(x)})
#save(seuObject_split, file = "Pancreas/output/full_matrix_before_integration.RData")

#' 
## ----SCTransform workflow regressing Phase & Integration-----------------------------------------------------------
options(future.globals.maxSize = 8000 * 1024^2)  
cat("SC_transfomr is processing now...")
for (i in 1:length(seuObject_split)) {
  message("Running SC Transform on : ", names(seuObject_split)[i])
  seuObject_split[[i]] <- SCTransform(seuObject_split[[i]],
                                      vars.to.regress = c("nUMI", "pMito", "pRibo", "S.Score","G2M.Score"), verbose = FALSE) 
  gc()}

# Select integration features across the groups
integ_features <- SelectIntegrationFeatures(object.list = seuObject_split, nfeatures = 3000)
head(integ_features,10)

# Preparing the SCT objects for integration
cat("Preparing the SCT object...")
seuObject_split <- PrepSCTIntegration(object.list = seuObject_split, anchor.features = integ_features)
cat("Done")

# Finding the anchors
cat("Finding integration anchors....")
start_time <- Sys.time()
integ_anchors <- FindIntegrationAnchors(object.list = seuObject_split,normalization.method = "SCT",anchor.features = integ_features)
end_time <- Sys.time()
cat("It took a while but finally done...it took", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")

length(integ_features)
table(integ_anchors@anchors[, "dataset1"])
table(integ_anchors@anchors[, "dataset2"])

# Integrate the data sets into a single Seurat object
library(future)
cat("Integrating dataset...")
seuObject_integrated <- IntegrateData(anchorset = integ_anchors,new.assay.name = "integrated",
  normalization.method = "SCT",dims = 1:50,k.weight = 100,sd.weight = 1,eps = 0.5,verbose = TRUE)
cat('done...')

#plan("multisession") 
# Set the integrated assay as default
DefaultAssay(seuObject_integrated) <- "integrated"
# Dimensionality reduction
cat("Starting dimensionality reduction...")
seuObject_integrated <- RunPCA(seuObject_integrated,features = NULL,weight.by.var = TRUE,ndims.print = 1:5,
                               nfeatures.print = 30,npcs = 50,reduction.name = "pca")
#pdf("PCA_dimensions.pdf", width = 9, height = 7)
ElbowPlot(seuObject_integrated, ndims = 50)
VizDimLoadings(seuObject_integrated, dims = 1:3, reduction = "pca")
DimHeatmap(seuObject_integrated, dims = 1:15, cells = 500, balanced = TRUE)
#invisible(dev.off())

# Using the first X component analysis in here
seuObject_integrated <- FindNeighbors(object=seuObject_integrated,reduction = "pca",dims = 1:40,  nn.eps = 0.5)
seuObject_integrated <- FindClusters(seuObject_integrated,resolution = seq(0.1, 1.2, by = 0.1),algorithm = 1, n.iter = 1000)  
cat("Done")

## ------------------------------------------------------------------------------------------------------------------
pdf("output/Cluster_tree_reductions.pdf", width = 12, height = 10)
clustree(seuObject_integrated@meta.data, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
invisible(dev.off())

#' 
## ----UMAP resolutions----------------------------------------------------------------------------------------------
set.seed(0708)
seuObject_integrated <- RunUMAP(seuObject_integrated, dims = 1:40, reduction = "pca")

resolutions <- seq(0.1, 1.2, by = 0.1)
dir.create("output/UMAP_by_resolution", showWarnings = FALSE)
# For loop to check resolutions
for (res in resolutions) {
  res_col <- paste0("integrated_snn_res.", res)
  Idents(seuObject_integrated) <- seuObject_integrated[[res_col]][,1]
  p <- DimPlot(seuObject_integrated, reduction = "umap", label = TRUE) + ggtitle(paste("Resolution", res)) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(filename = paste0("output/UMAP_by_resolution/UMAP_res_", res, ".png"),plot = p,width = 8, height = 6, dpi = 300)
  print(p)}

#' 
## ------------------------------------------------------------------------------------------------------------------
FeaturePlot(seuObject_integrated, features = "Aqp4", min.cutoff = "q10", max.cutoff = "q95")
DimPlot(seuObject_integrated, reduction = "umap", group.by = "integrated_snn_res.0.2", label =T)
DimPlot(seuObject_integrated, reduction = "umap", group.by = "integrated_snn_res.0.2", label =T, split.by = "SampleID")
DimPlot(seuObject_integrated, reduction = "umap", group.by = "integrated_snn_res.0.2", label =T, split.by = "Status")
DimPlot(seuObject_integrated, reduction = "umap", group.by = "integrated_snn_res.0.2", label =T, split.by = "Phase")

# Normalizing RNA counts for doublets detection 
DefaultAssay(seuObject_integrated) <- "RNA"
seuObject_integrated <- NormalizeData(object = seuObject_integrated,normalization.method = "LogNormalize",scale.factor = 10000)

#' 
#' #---------------------
#' # Doublets detection
## ------------------------------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
library(tidyverse);library(ggrepel);library(emmeans);library(SingleCellExperiment);library(scater)
;library(BiocParallel);library(ggpubr);library(speckle);library(magrittr);library(broom)
;library(muscat);library(Seurat);library(clustree);library(leiden)
;library(data.table);library(cowplot);library(scDblFinder);library(BiocSingular);library(scds)
})

#' 
## ------------------------------------------------------------------------------------------------------------------
# Loading Seurat Object
load("E:/Radiation_sc_paper_output/Seurat_integrated.RData")
# DietSeurat to reduce memmory burden
seuObject_integrated_slim <- DietSeurat(seuObject_integrated,counts = TRUE,data = TRUE,scale.data = FALSE,assays="RNA",
                                        dimreducs = c("pca","umap"))
DimPlot(seuObject_integrated_slim, reduction = "umap")

# Join the layers 
seuObject_integrated_slim_Joined<-JoinLayers(seuObject_integrated_slim)
# Transforming into a SingleCellExperiment 
sce <- as.SingleCellExperiment(seuObject_integrated_slim_Joined)

# Doubleting by SampleID
cat("Start running the doubleting detection.....")
set.seed(0708)
sce <- scDblFinder(sce,samples="SampleID", BPPARAM=BiocParallel::SnowParam(workers = 8), nfeatures = 3000,dims = 50,dbr.sd = 1)
cat("Finished..!")

# Adding doublet identification on the metadata
seuObject_integrated_slim_Joined@meta.data$Doublets <- sce$scDblFinder.class
# Subseting doublets from singlets
seuObject_slim_nodoub <- subset(seuObject_integrated_slim_Joined, subset = Doublets == "singlet")
table(sce$scDblFinder.class) ## number of singlets and doublets 
table(sce$scDblFinder.class, sce$SampleID) ## number of singlets and doublets in each condition 
round(100 * prop.table(table(sce$scDblFinder.class, sce$SampleID), margin = 2), 2)

#' 
## ------------------------------------------------------------------------------------------------------------------
DimPlot(seuObject_integrated, reduction = "umap", label = T, group.by = "integrated_snn_res.0.2")
DimPlot(seuObject_slim_nodoub, reduction = "umap", label =T, group.by = "integrated_snn_res.0.2")
DimPlot(seuObject_slim_nodoub, reduction = "umap", label =T, group.by = "integrated_snn_res.0.2", split.by = "SampleID")


DefaultAssay(seuObject_slim_nodoub)<-"SCT"
# Viz
dittoDimPlot(seuObject_slim_nodoub, var = "integrated_snn_res.0.2", reduction.use = "umap",  size = 0.5,opacity = 1,show.others = T,do.label = TRUE)

#UMAP plots of singlets
resolutions <- seq(0.1, 1.2, by = 0.1)
dir.create("output/nodoub_UMAP_by_resolution", showWarnings = FALSE)
# For loop to check resolutions
for (res in resolutions) {
  res_col <- paste0("integrated_snn_res.", res)
  Idents(seuObject_slim_nodoub) <- seuObject_slim_nodoub[[res_col]][,1]
  p <- DimPlot(seuObject_slim_nodoub, reduction = "umap", label = TRUE) +ggtitle(paste("Resolution", res)) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(filename = paste0("output/nodoub_UMAP_by_resolution/UMAP_res_", res, ".png"), plot = p,width = 8, height = 6, dpi = 300)
  print(p)}

#' 
#' ### Finding markers of each cluster ###
## ------------------------------------------------------------------------------------------------------------------
### Running presto 
DefaultAssay(seuObject_slim_nodoub) <- "RNA"

seuObject_slim_nodoub <- NormalizeData(object = seuObject_slim_nodoub,normalization.method = "LogNormalize", scale.factor = 10000)
seuObject_slim_nodoub<- JoinLayers(seuObject_slim_nodoub)
Idents(seuObject_slim_nodoub) <- seuObject_slim_nodoub$integrated_snn_res.0.4

wilcoxauc.Seurat <- function(X,group_by = NULL,assay = "data",groups_use = NULL, seurat_assay = "RNA",
    ...
) {
    requireNamespace("Seurat")
    X_matrix <- Seurat::GetAssayData(X, assay = seurat_assay, layer = assay)
    if (is.null(group_by)) {
        y <- Seurat::Idents(X)
    } else {
        y <- Seurat::FetchData(X, group_by) %>% unlist %>% as.character()
    }
    wilcoxauc(X_matrix, y, groups_use)
}

all_markers_clustID <- wilcoxauc.Seurat(seuObject_slim_nodoub, group_by ='integrated_snn_res.0.4')
unique(all_markers_clustID$group)
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>% dplyr::filter(padj < 0.05, logFC > 0.3)
top20 <- presto::top_markers(all_markers.Sign, n = 20, auc_min = 0.5, pval_max = 0.05)

## Running Wilcoxau analysis 
#all_markers_clustID <- presto::wilcoxauc(seuObject_slim_nodoub, 'integrated_snn_res.0.3', assay = 'data')
#unique(all_markers_clustID$group)
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
#all_markers.Sign <- all_markers_clustID %>%dplyr::filter(padj < 0.05, logFC > 0)
#top20 <- presto::top_markers(all_markers.Sign,n = 20,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(all_markers.Sign,
                     file = "Refinment/0.4_PrestoByCluster_Filteredmarkers_padjLT05_logfcGT0.xlsx",
                     colNames = TRUE,rowNames = FALSE,borders = "columns", sheetName="Markers")
openxlsx::write.xlsx(top20,
                     file = "Refinment/0.4_Mice_dataset_top20.xlsx",colNames = TRUE,
                     rowNames = FALSE,borders = "columns",sheetName="Markers")

#' 
## ----General annotation and Viz------------------------------------------------------------------------------------
#Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$integrated_snn_res.0.2
DimPlot(seuObject_slim_nodoub, group.by = "celltype")
DimPlot(seuObject_slim_nodoub, group.by = "Phase")
FeaturePlot(seuObject_slim_nodoub, features = "C1ql1", min.cutoff = "q10", max.cutoff = "q95")

Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$integrated_snn_res.0.3
## Initial celltype
labels <- c("0" = "Post mitotic neurons",
            "1" = "Proliferating Tumor",
            "2" = "Quiescent tumor",
            "3" = "Proliferating Tumor",
            "4" = "Post mitotic neurons",
            "5" = "Proliferating Tumor",
            "6" = "Microglia", 
            "7" = "Post mitotic neurons",
            "8" = "Astrocytes", 
            "9" = "Inmune cells",
            "10" = "Oligodendrocytes",
            "11" = "Neutrophils",
            "12" = "Endothelial cells", 
            "13" = "T cells",
            "14" = "Choroid plexus",
            "14" = "Choroid plexus")
Idents(seuObject_slim_nodoub) <- "integrated_snn_res.0.2"
seuObject_slim_nodoub@meta.data$celltype<- labels[as.character(Idents(seuObject_slim_nodoub))]


## More specific celltypes based on integrated_snn_res.0.4
labels <- c("0" = "Differentiating granule-neuron-like tumor cells",
            "1" = "GNP-like progenitor translation-high",
            "2" = "maturing granule-neuron-like",
            "3" = "Differentiating granule-neuron-like",
            "4" = "Cycling tumor S/G2M phase",
            "5" = "Cycling tumor G2M phase",
            "6" = "Tumor S-phase", 
            "7" = "GN transitional tumor cells",
            "8" = "Tumor-associated macrophages", 
            "9" = "Cycling tumor G2M phase",
            "10" = "Astrocytes",
            "11" = "BAM perivascular macrophages",
            "12" = "GNP-like neuronal progenitors", 
            "13" = "Immature granule neurons",
            "14" = "oligodendrocyte lineage",
            "15" = "Neutrophils",
            "16" = "Endothelial cells",
            "17" = "Mature cerebellar neurons",
            "18" = "Choroid plexus epithelium",
            "19" = "DAM microglia",
            "20" = "T lymphocytes")
Idents(seuObject_slim_nodoub) <- "integrated_snn_res.0.4"
seuObject_slim_nodoub@meta.data$celltype1<- labels[as.character(Idents(seuObject_slim_nodoub))]
#colnames(seuObject_slim_nodoub@meta.data)

## Dotplot based on celltype1
pdf("output/Markers/dotplot.pdf", width = 7, height = 10)
DotPlot(seuObject_slim_nodoub, features = c("Pf4",'Ms4a7','F13a1','Mrc1', # BAMs
                                            'Apoe','Trem2','Cst7','Lpl','Gpnmb') #DAM-like
                                           , cols = "RdYlBu") + RotatedAxis() + coord_flip()
invisible(dev.off())
 
### General cellular types
labels <- c("0" = "Tumor differentiating",
            "1" = "Tumor progenitor",
            "2" = "Tumor differentiating",
            "3" = "Tumor differentiating",
            "4" = "Tumor proliferating",
            "5" = "Tumor proliferating",
            "6" = "Tumor proliferating", 
            "7" = "Tumor differentiating",
            "8" = "TAM macrophages/ferritin",  #
            "9" = "Tumor proliferating",
            "10" = "Astrocytes",
            "11" = "BAM macrophages/monocytic", #
            "12" = "Tumor progenitor", 
            "13" = "Tumor differentiating",
            "14" = "Oligodendrocyte lineage",
            "15" = "Neutrophils",
            "16" = "Endothelial cells",
            "17" = "Mature cerebellar neurons",
            "18" = "Choroid plexus epithelium",
            "19" = "DAM Microglia",  #
            "20" = "T lymphocytes")
Idents(seuObject_slim_nodoub) <- "integrated_snn_res.0.4"
seuObject_slim_nodoub@meta.data$celltype2<- labels[as.character(Idents(seuObject_slim_nodoub))]

# Viz
Idents(seuObject_slim_nodoub)<-"celltype2"
pdf("Refinment/UMAP.pdf", width = 10, height = 8)
DimPlot2(seuObject_slim_nodoub, label = T, box = TRUE, label.color = "black", repel = TRUE, theme =  NoAxes()) + theme_umap_arrows() 
invisible(dev.off())

pdf("Refinment/Status_UMAP.pdf", width = 13, height = 7)
DimPlot2(seuObject_slim_nodoub, label = T, box = TRUE, 
         split.by = "Status",label.color = "black", repel = TRUE, theme =  NoAxes()) + theme_umap_arrows() 
invisible(dev.off())

# Proportions
ClusterDistrBar(origin = seuObject_slim_nodoub$Status, cluster = seuObject_slim_nodoub$celltype2)
ClusterDistrBar(origin = seuObject_slim_nodoub$Status, cluster = seuObject_slim_nodoub$celltype2, percent = FALSE)
ClusterDistrBar(origin = seuObject_slim_nodoub$Status, cluster = seuObject_slim_nodoub$celltype2, rev = TRUE, normalize = TRUE)

pdf("Refinment/Percentage_clusters.pdf", width = 11, height = 7)
ClusterDistrPlot(origin = seuObject_slim_nodoub$SampleID,cluster = seuObject_slim_nodoub$celltype2,condition = seuObject_slim_nodoub$Status, hide.ns = F, stat.method = "t.test", cols = c("red","blue","black"))
invisible(dev.off())

#' 
#' 
#' # Running Harmony workflow #
#' #---- Not used in output results 
## ------------------------------------------------------------------------------------------------------------------
DefaultAssay(seuObject_slim_nodoub) <- "SCT"
seuOject_nodoub_withHarmony <- RunHarmony(seuObject_slim_nodoub,assay.use= "SCT",group.by.vars = "Sample")

Reductions(seuOject_nodoub_withHarmony)
seuOject_nodoub_withHarmony <- RunUMAP(seuOject_nodoub_withHarmony, reduction = "harmony", dims = 1:40)
seuOject_nodoub_withHarmony <- FindNeighbors(seuOject_nodoub_withHarmony, reduction = "harmony", dims = 1:40)
seuOject_nodoub_withHarmony@graphs
seuOject_nodoub_withHarmony <- FindClusters(seuOject_nodoub_withHarmony,
                                                 resolution = c(0.1,0.2,0.3,0.4,0.5,0.6, 0.7, 0.8, 0.9, 1),
                                                 algorithm = 1, n.iter = 1000,graph.name = "integrated_snn")

DimPlot(object = seuOject_nodoub_withHarmony, reduction = "umap", group.by = "integrated_snn_res.0.3", label = TRUE, pt.size = 0.5) + theme(legend.position="right")
dittoDimPlot(seuOject_nodoub_withHarmony, var = "integrated_snn_res.0.3", reduction.use = "umap", size = 0.5,opacity = 1,show.others = T,do.label = TRUE)
#DimPlot(object = seuOject_nodoub_withHarmony, reduction = "umap", group.by = "integrated_snn_res.0.2", label = TRUE, pt.size = 0.5, split.by = "Sample") + theme(legend.position="right")
#FeaturePlot(object = seuOject_nodoub_withHarmony, features = c("C1qc","C1qb"))
resolutions <- seq(0.1, 1, by = 0.1)
dir.create("Pancreas/output/Harmony_UMAP_by_resolution", showWarnings = FALSE)
### For loop to check resolutions
for (res in resolutions) {
  res_col <- paste0("integrated_snn_res.", res)
  Idents(seuOject_nodoub_withHarmony) <- seuOject_nodoub_withHarmony[[res_col]][,1]
  p <- DimPlot(seuOject_nodoub_withHarmony, reduction = "umap", label = TRUE) +ggtitle(paste("Resolution", res)) +theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(filename = paste0("Pancreas/output/Harmony_UMAP_by_resolution/UMAP_res_", res, ".png"),plot = p,width = 8, height = 6, dpi = 300)
  print(p)}


#' 
## ----Saving--------------------------------------------------------------------------------------------------------
save(seuObject_integrated, file = "E:/Radiation_sc_paper_output/Seurat_integrated.RData")
save(seuObject_slim_nodoub, file = "Single_cell_analysis_radiation_paper/output/Seurat_nodoub_integrated.RData")
save(seuOject_nodoub_withHarmony, file = "Pancreas/output/Harmony_Seurat_nodoub_integrated.RData")

#' 
#' 
#' # DownStream analysis
#' #---------------------
## ------------------------------------------------------------------------------------------------------------------
## Running wilcoxox.Seurat function on celltype2
DefaultAssay(seuObject_slim_nodoub)<-"RNA"
Idents(seuObject_slim_nodoub)<- seuObject_slim_nodoub$celltype2

all_markers_clustID <- wilcoxauc.Seurat(seuObject_slim_nodoub, group_by ='celltype2')
unique(all_markers_clustID$group)
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>% dplyr::filter(padj < 0.05, logFC > 0.3)
#top20 <- presto::top_markers(all_markers.Sign, n = 20, auc_min = 0.5, pval_max = 0.05)

unique(all_markers.Sign$group)

toppData <-  toppFun(all_markers.Sign,topp_categories = NULL,cluster_col = "group", gene_col = "feature", p_val_col = "padj",logFC_col = "logFC", pval_cutoff = 0.05, min_genes = 10,max_genes = 500,max_results = 50)

toppPlot(toppData,category = "GeneOntologyMolecularFunction",num_terms = 10, p_val_adj = "BH", p_val_display = "log",save = TRUE, save_dir = "Refinment/pseudobulk_moderate_GO",width = 5, height = 6)

toppPlot(toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,p_val_adj = "BH",p_val_display = "log",save = TRUE,save_dir = "Refinment/pseudobulk_moderate_GO",width = 5,height = 6)

save(toppData, file ="Refinment/pseudobulk_moderate_GO/topData.RData")

## ----r session-info------------------------------------------------------------------------------------------------
sessioninfo::session_info()

#' 
#' 
## ------------------------------------------------------------------------------------------------------------------
knitr::purl(
  "E:/Blanco_Lab/Single_cell_analysis_radiation_paper/Radiation_workflow.Rmd",output = "Radiation_workflow.R",documentation = 2)

si <- capture.output(sessionInfo())
si_comment <- paste0("# ", si)
write(c("\n\n# Session Information\n",si_comment),file = "Radiation_workflow.R",append = TRUE)

#' 


# Session Information

# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
#   LAPACK version 3.12.1
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] scds_1.24.0                 BiocSingular_1.26.1         scDblFinder_1.22.0          cowplot_1.2.0              
#  [5] data.table_1.18.2.1         leiden_0.4.3.1              clustree_0.5.1              ggraph_2.2.2               
#  [9] Seurat_5.5.0                SeuratObject_5.4.0          sp_2.2-1                    muscat_1.22.0              
# [13] broom_1.0.12                magrittr_2.0.4              speckle_1.8.0               ggpubr_0.6.3               
# [17] BiocParallel_1.44.0         scater_1.36.0               scuttle_1.18.0              SingleCellExperiment_1.30.1
# [21] SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.3        
# [25] IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.54.1         generics_0.1.4             
# [29] MatrixGenerics_1.20.0       matrixStats_1.5.0           emmeans_2.0.3               ggrepel_0.9.8              
# [33] lubridate_1.9.5             forcats_1.0.1               stringr_1.6.0               dplyr_1.2.1                
# [37] purrr_1.2.2                 readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1               
# [41] ggplot2_4.0.2               tidyverse_2.0.0            
# 
# loaded via a namespace (and not attached):
#   [1] dichromat_2.0-0.1        progress_1.2.3           goftest_1.2-3            Biostrings_2.78.0       
#   [5] vctrs_0.7.2              spatstat.random_3.4-5    digest_0.6.39            png_0.1-9               
#   [9] corpcor_1.6.10           shape_1.4.6.1            deldir_2.0-4             parallelly_1.46.1       
#  [13] MASS_7.3-65              reshape2_1.4.5           httpuv_1.6.17            foreach_1.5.2           
#  [17] withr_3.0.2              xfun_0.55                survival_3.8-3           memoise_2.0.1           
#  [21] ggbeeswarm_0.7.3         Seqinfo_1.0.0            zoo_1.8-15               GlobalOptions_0.1.4     
#  [25] gtools_3.9.5             pbapply_1.7-4            Formula_1.2-5            prettyunits_1.2.0       
#  [29] promises_1.5.0           otel_0.2.0               httr_1.4.8               rstatix_0.7.3           
#  [33] restfulr_0.0.16          globals_0.19.1           fitdistrplus_1.2-6       rstudioapi_0.18.0       
#  [37] UCSC.utils_1.4.0         miniUI_0.1.2             curl_7.0.0               ScaledMatrix_1.18.0     
#  [41] polyclip_1.10-7          GenomeInfoDbData_1.2.14  SparseArray_1.10.8       xtable_1.8-8            
#  [45] doParallel_1.0.17        evaluate_1.0.5           S4Arrays_1.10.1          hms_1.1.4               
#  [49] irlba_2.3.7              colorspace_2.1-2         ROCR_1.0-12              reticulate_1.46.0       
#  [53] spatstat.data_3.1-9      lmtest_0.9-40            later_1.4.8              viridis_0.6.5           
#  [57] lattice_0.22-7           spatstat.geom_3.7-3      future.apply_1.20.2      scattermore_1.2         
#  [61] XML_3.99-0.23            RcppAnnoy_0.0.23         pillar_1.11.1            nlme_3.1-168            
#  [65] iterators_1.0.14         caTools_1.18.3           compiler_4.5.1           beachmat_2.26.0         
#  [69] RSpectra_0.16-2          stringi_1.8.7            tensor_1.5.1             minqa_1.2.8             
#  [73] GenomicAlignments_1.44.0 plyr_1.8.9               crayon_1.5.3             abind_1.4-8             
#  [77] BiocIO_1.18.0            blme_1.0-7               locfit_1.5-9.12          graphlayouts_1.2.3      
#  [81] sandwich_3.1-1           codetools_0.2-20         GetoptLong_1.1.1         plotly_4.12.0           
#  [85] remaCor_0.0.20           mime_0.13                splines_4.5.1            circlize_0.4.18         
#  [89] Rcpp_1.1.1               fastDummies_1.7.6        knitr_1.51               here_1.0.2              
#  [93] clue_0.3-68              lme4_2.0-1               listenv_0.10.1           Rdpack_2.6.6            
#  [97] ggsignif_0.6.4           estimability_1.5.1       Matrix_1.7-3             statmod_1.5.1           
# [101] tzdb_0.5.0               fANCOVA_0.6-1            tweenr_2.0.3             pkgconfig_2.0.3         
# [105] tools_4.5.1              cachem_1.1.0             RhpcBLASctl_0.23-42      rbibutils_2.4.1         
# [109] viridisLite_0.4.3        numDeriv_2016.8-1.1      fastmap_1.2.0            rmarkdown_2.31          
# [113] scales_1.4.0             grid_4.5.1               ica_1.0-3                Rsamtools_2.24.1        
# [117] patchwork_1.3.2          coda_0.19-4.1            dotCall64_1.2            carData_3.0-6           
# [121] RANN_2.6.2               farver_2.1.2             reformulas_0.4.4         aod_1.3.3               
# [125] tidygraph_1.3.1          mgcv_1.9-3               yaml_2.3.12              rtracklayer_1.68.0      
# [129] cli_3.6.6                lifecycle_1.0.5          uwot_0.2.4               glmmTMB_1.1.14          
# [133] mvtnorm_1.3-6            bluster_1.18.0           sessioninfo_1.2.3        backports_1.5.1         
# [137] timechange_0.4.0         gtable_0.3.6             rjson_0.2.23             ggridges_0.5.7          
# [141] progressr_0.19.0         pROC_1.19.0.1            parallel_4.5.1           limma_3.66.0            
# [145] jsonlite_2.0.0           edgeR_4.8.2              RcppHNSW_0.6.0           bitops_1.0-9            
# [149] xgboost_3.2.1.1          Rtsne_0.17               spatstat.utils_3.2-3     BiocNeighbors_2.2.0     
# [153] metapod_1.16.0           dqrng_0.4.1              spatstat.univar_3.2-0    pbkrtest_0.5.5          
# [157] lazyeval_0.2.3           shiny_1.13.0             htmltools_0.5.9          sctransform_0.4.3       
# [161] rappdirs_0.3.4           glue_1.8.0               spam_2.11-3              XVector_0.50.0          
# [165] RCurl_1.98-1.18          rprojroot_2.1.1          scran_1.36.0             gridExtra_2.3           
# [169] EnvStats_3.1.0           boot_1.3-31              igraph_2.2.3             variancePartition_1.38.1
# [173] TMB_1.9.21               R6_2.6.1                 DESeq2_1.48.2            gplots_3.3.0            
# [177] cluster_2.1.8.2          nloptr_2.2.1             DelayedArray_0.36.0      tidyselect_1.2.1        
# [181] vipor_0.4.7              ggforce_0.5.0            car_3.1-5                future_1.70.0           
# [185] rsvd_1.0.5               KernSmooth_2.23-26       S7_0.2.1                 htmlwidgets_1.6.4       
# [189] ComplexHeatmap_2.26.0    RColorBrewer_1.1-3       rlang_1.1.7              spatstat.sparse_3.2-0   
# [193] spatstat.explore_3.8-0   lmerTest_3.2-1           beeswarm_0.4.0          
