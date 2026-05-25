#' ---
#' title: "TME reclustering"
#' author: "Carlos Alfaro"
#' date: "2025-07-09"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ------------------------------------------------------------------------------------------------------------------
library(SingleCellExperiment);library(Seurat);library(tidyverse);library(Matrix);library(scales);library(cowplot);library(RCurl);library(openxlsx);library(knitr);library(monocle3);library(SeuratWrappers);library(scToppR);library(SeuratExtend);library(presto)

#' 
## ----Subset--------------------------------------------------------------------------------------------------------
load("Single_cell_analysis_radiation_paper/output/Seurat_nodoub_integrated.RData")
DimPlot(seuObject_slim_nodoub, reduction= "umap", group.by  = "celltype2")

## Idents
Idents(seuObject_slim_nodoub) <- seuObject_slim_nodoub$celltype2
DefaultAssay(seuObject_slim_nodoub)<- "RNA"

## Subset TME
unique(seuObject_slim_nodoub$celltype2)
TME <-subset(x = seuObject_slim_nodoub, 
             idents = c("TAM macrophages/ferritin", "T lymphocytes", "Neutrophils", "Endothelial cells","BAM macrophages/monocytic","Choroid plexus epithelium","Oligodendrocyte lineage","Astrocytes","DAM Microglia"))

DimPlot(TME, reduction = "umap",label = TRUE, pt.size = 0.5)

## ----Integrating workflow------------------------------------------------------------------------------------------
# Integration workflow
seuObject_split <- SplitObject(TME, split.by = "SampleID")
options(future.globals.maxSize = 8000 * 1024^2)  # 2GB
cat("SC_transfomr is processing now...")
for (i in 1:length(seuObject_split)) {
  message("Running SC Transform on : ", names(seuObject_split)[i])
  seuObject_split[[i]] <- SCTransform(seuObject_split[[i]],vars.to.regress = c("nUMI", "pMito", "pRibo"),verbose = FALSE) 
  gc()}
# Select integration features across the groups
integ_features <- SelectIntegrationFeatures(object.list = seuObject_split, nfeatures = 3000)
head(integ_features,10)
# Preparing the SCT objects for integration
cat("Preparing the SCT object...")
seuObject_split <- PrepSCTIntegration(object.list = seuObject_split,anchor.features = integ_features)
cat("Done")
### Finding the anchors
cat("Finding integration anchors....")
start_time <- Sys.time()
integ_anchors <- FindIntegrationAnchors(object.list = seuObject_split,normalization.method = "SCT",anchor.features = integ_features)
end_time <- Sys.time()
cat("Finished...!")
# Integrate the data sets into a single Seurat object
### This step is very memory   consuming so lets go to sequential
library(future)
cat("Integrating dataset...")
seuObject_integrated <- IntegrateData(anchorset = integ_anchors,new.assay.name = "integrated",normalization.method = "SCT",
  dims = 1:50,k.weight = 100,sd.weight = 1,eps = 0.5,verbose = TRUE)
cat('done...')
# Set the integrated assay as default
DefaultAssay(seuObject_integrated) <- "integrated"
# Dimensionality reduction
cat("Starting dimensionality reduction...")
seuObject_integrated <- RunPCA(seuObject_integrated,features = NULL,weight.by.var = TRUE,ndims.print = 1:5,nfeatures.print = 30,npcs = 50,reduction.name = "pca")

#pdf("Astrocytes_reclustering/PCA_dimensions.pdf", width = 9, height = 7)
ElbowPlot(seuObject_integrated, ndims = 50)
VizDimLoadings(seuObject_integrated, dims = 1:2, reduction = "pca")
DimHeatmap(seuObject_integrated, dims = 1:15, cells = 800, balanced = TRUE)
#invisible(dev.off())
# Using the first 25 component analysis in here
seuObject_integrated <- FindNeighbors(object=seuObject_integrated,reduction = "pca",dims = 1:30,nn.eps = 0.5)
seuObject_integrated <- FindClusters(seuObject_integrated,resolution = seq(0.1, 1.2, by = 0.1),algorithm = 1, n.iter = 1000)

# UMAP viz
set.seed(7081998)
seuObject_integrated <- RunUMAP(seuObject_integrated, dims = 1:30, reduction = "pca")
resolutions <- seq(0.1, 1.2, by = 0.1)
dir.create("Refinment/TME/UMAP_by_resolution", showWarnings = FALSE)
### For loop to check resolutions
for (res in resolutions) {
  res_col <- paste0("integrated_snn_res.", res)
  Idents(seuObject_integrated) <- seuObject_integrated[[res_col]][,1]
  p <- DimPlot(seuObject_integrated, reduction = "umap", label = TRUE) +ggtitle(paste("Resolution", res)) +theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(filename = paste0("Refinment/TME/UMAP_by_resolution/UMAP_res_", res, ".png"),plot = p,width = 8, height = 6, dpi = 300)
  print(p)}

# Normalize
DefaultAssay(seuObject_integrated) <- "RNA"
seuObject_integrated <- NormalizeData(object = seuObject_integrated,normalization.method = "LogNormalize",scale.factor = 10000)

## ----Find Markers--------------------------------------------------------------------------------------------------
save(seuObject_integrated, file = "Refinment/TME/seuObject_integrated.RData")

## Presto markers
## Identity resolution I want 
seuObject_integrated<- JoinLayers(seuObject_integrated)
Idents(seuObject_integrated) <- "integrated_snn_res.0.6"
## Running Wilcoxau analysis \
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
all_markers_clustID <- wilcoxauc.Seurat(seuObject_integrated, group_by ='integrated_snn_res.0.6')
unique(all_markers_clustID$group)
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>% dplyr::filter(padj < 0.05, logFC > 0.3)
top20 <- presto::top_markers(all_markers.Sign, n = 20, auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(all_markers.Sign,
                     file = "Refinment/TME/0.6PrestoByCluster_Filteredmarkers_padjLT05_logfcGT0.xlsx",
                     colNames = TRUE,rowNames = FALSE,borders = "columns",sheetName="Markers")
openxlsx::write.xlsx(top20,
                     file = "Refinment/TME/0.6PrestoByCluster_Top20.xlsx",
                     colNames = TRUE,rowNames = FALSE,borders = "columns",sheetName="Markers")

#' 
#' 
## ----Annotation----------------------------------------------------------------------------------------------------
load("Refinment/TME/seuObject_integrated.RData")
DefaultAssay(seuObject_integrated)<- "RNA"
### More specific
labels <- c("0" = "Bergmann Astrocytes",
            "1" = "Homeostatic microglia",
            "2" = "MHC-II BAMs",
            "3" = "Undetermined myeloid cells",
            "4" = "GNP cells",
            "5" = "DAM1",
            "6" = "Inflammatory neutrophils", 
            "7" = "Lipid macrophages",
            "8" = "OPCs", 
            "9" = "Ccl3 microglia", 
            "10" = "Endothelial cells",
            "11" = "Reactive astrocytes A2",
            "12" = "Metallothionein astrocytes A1",
            "13" = "cDC2",
            "14" = "Fibroblasts",
            "15" = "Cd3 T lymphocytes",
            "16" = "Mature oligodendrocytes",
            "17" = "Newly-formed oligodendrocytes",
            "18" = "Choroid plexus epithelium",
            "19" = "Differentiating GNP cells",
            "20" = "Migratory mature DCs"
            )
Idents(seuObject_integrated) <- "integrated_snn_res.0.6"
#DimPlot(seuObject_integrated, group.by = "celltype")
seuObject_integrated@meta.data$celltype3<- labels[as.character(Idents(seuObject_integrated))]
DimPlot(seuObject_integrated, group.by = "celltype3")

## ----Dotplot-------------------------------------------------------------------------------------------------------
markers_by_cluster <- list(
  `0`  = c("Fabp7","Gpr37l1","Hopx","Slc1a3"),
  `1`  = c("Hexb","Fcrls","C1qa","Selenop"),
  `2`  = c("Cd74","H2-Aa","Ms4a7","Tgfbi"),
  `3`  = c("Ftl1","Fth1","Fabp5","Tmsb4x"),
  `4`  = c("Sfrp1","Nfib","Ccnd2","Igfbpl1"),
  `5`  = c("Cst7","Trem2","Apoe","Ctsd"),
  `6`  = c("S100a8","S100a9","Il1b","Cxcl2"),
  `7`  = c("Spp1","Plin2","Gpnmb","Stab1"),
  `8`  = c("Pdgfra","Cspg5","Ptprz1","Lhfpl3"),
  `9`  = c("Ccl3","Ccl4","Cx3cr1","Csf1r"),
  `10` = c("Cldn5","Ly6c1","Flt1","Pglyrp1"),
  `11` = c("Clu","Mlc1","Cryab","Ntrk2"),
  `12` = c("Mt1","Mt2","Timp1","Cp"),
  `13` = c("Napsa","Plbd1","H2-Ab1","Klrd1"),
  `14` = c("Dcn","Lum","Col1a2","Vtn"),
  `15` = c("Cd3d","Cd3g","Trbc2","Ccl5"),
  `16` = c("Plp1","Mbp","Mag","Cldn11"),
  `17` = c("Bcas1","Olig1","Olig2","S100a1"),
  `18` = c("Ttr","Folr1","Kl","Kcnj13"),
  `19` = c("Neurod1","Stmn2","Sox11","Zic1"),
  `20` = c("Ccr7","Fscn1","Tmem123","Bcl2a1b")
)
features <- unique(unlist(markers_by_cluster, use.names = FALSE))
Idents(seuObject_integrated)<- "celltype3"
#Idents(seuObject_integrated) <- factor(Idents(seuObject_integrated), levels = as.character(0:20))
p <- DotPlot(seuObject_integrated, features = features,cols = "RdBu",dot.scale = 4, cluster.idents = FALSE) +
     RotatedAxis() +theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),axis.text.y = element_text(size = 8)) + labs(x = NULL, y = NULL, title = "SmoA1 TME — 4 markers per cluster")
print(p)

#' 
## ----Seurat Extend-------------------------------------------------------------------------------------------------
Idents(seuObject_integrated)<- "celltype3"

pdf("Refinment/TME/Sox2_increase.pdf", width = 10, height = 8)
DimPlot2(seuObject_integrated, features = c("celltype3", "Sox2"), split.by = "Status", ncol = 1)
invisible(dev.off())


pdf("Refinment/TME/UMAP_general.pdf", width = 10, height = 8)
DimPlot2(seuObject_integrated, label = F , box = TRUE, label.color = "black", repel = TRUE, theme =  NoAxes()) + theme_umap_arrows()
invisible(dev.off())

#DefaultAssay(seuObject_integrated)<-"SCT"
#SCpubr::do_FeaturePlot(sample = seuObject_integrated,features = "Aqp4",order = T,label = F,enforce_symmetry = TRUE,min.cutoff = 0,max.cutoff = 4,split.by = "Status")

### another way 
table(seuObject_integrated$celltype3, seuObject_integrated$Status) %>%prop.table(margin = 2) %>%  round(3) * 100
ClusterDistrBar(origin = seuObject_integrated$Status, cluster = seuObject_integrated$celltype3)
ClusterDistrBar(origin = seuObject_integrated$Status, cluster = seuObject_integrated$celltype3, percent = FALSE)
ClusterDistrBar(origin = seuObject_integrated$Status, cluster = seuObject_integrated$celltype3, rev = TRUE, normalize = TRUE)

pdf("Refinment/TME/proportionsV1.pdf", width = 13, height = 8)
ClusterDistrPlot(origin = seuObject_integrated$SampleID,cluster = seuObject_integrated$celltype3,condition = seuObject_integrated$Status, hide.ns = F, stat.method = "t.test", cols = c("red","blue","black"))
invisible(dev.off())


## ----Gene ontology-------------------------------------------------------------------------------------------------
DefaultAssay(seuObject_integrated)<- "RNA"
Idents(seuObject_integrated)<- "celltype3"
seuObject_integrated<- JoinLayers(seuObject_integrated)
all_markers_clustID <- wilcoxauc.Seurat(seuObject_integrated, group_by ='celltype3')
unique(all_markers_clustID$group)
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>% dplyr::filter(padj < 0.05, logFC > 0.3)
top20 <- presto::top_markers(all_markers.Sign, n = 20, auc_min = 0.5, pval_max = 0.05)

unique(all_markers.Sign$group)
toppData <-  toppFun(all_markers.Sign,topp_categories = NULL,cluster_col = "group",gene_col = "feature",p_val_col = "padj", pval_cutoff = 0.05,logFC_col = "logFC",min_genes = 10,max_genes = 500,max_results = 50)

toppPlot(toppData, category = "GeneOntologyMolecularFunction",num_terms = 10,p_val_adj = "BH", p_val_display = "log",save = TRUE,save_dir = "Refinment/TME/GO/pseudobulk_moderate_GO",width = 5,height = 6)

toppPlot(toppData,category = "GeneOntologyBiologicalProcess",num_terms = 10,p_val_adj = "BH", p_val_display = "log",save = TRUE,save_dir = "Refinment/TME/GO/pseudobulk_moderate_GO",width = 5,height = 6)

save(toppData, file ="Refinment/TME/GO/topData.RData")
save(seuObject_integrated, file = "Refinment/TME/seuObject_integrated.RData")

## ----r session-info------------------------------------------------------------------------------------------------
sessioninfo::session_info()

#' 
## ------------------------------------------------------------------------------------------------------------------
knitr::purl(
  "E:/Blanco_Lab/Single_cell_analysis_radiation_paper/TME_workflow.Rmd",output = "TME_workflow.R",documentation = 2)

si <- capture.output(sessionInfo())
si_comment <- paste0("# ", si)
write(c("\n\n# Session Information\n",si_comment),file = "TME_workflow.R",append = TRUE)



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
#  [1] presto_1.0.0                data.table_1.18.2.1         Rcpp_1.1.1                  SeuratExtend_1.2.5         
#  [5] SeuratExtendData_0.3.0      scToppR_0.99.4              SeuratWrappers_0.4.0        monocle3_1.4.26            
#  [9] knitr_1.51                  openxlsx_4.2.8.1            RCurl_1.98-1.18             cowplot_1.2.0              
# [13] scales_1.4.0                Matrix_1.7-3                lubridate_1.9.5             forcats_1.0.1              
# [17] stringr_1.6.0               dplyr_1.2.1                 purrr_1.2.2                 readr_2.2.0                
# [21] tidyr_1.3.2                 tibble_3.3.1                ggplot2_4.0.2               tidyverse_2.0.0            
# [25] Seurat_5.5.0                SeuratObject_5.4.0          sp_2.2-1                    SingleCellExperiment_1.30.1
# [29] SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.3        
# [33] IRanges_2.44.0              S4Vectors_0.48.0            BiocGenerics_0.54.1         generics_0.1.4             
# [37] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.23        splines_4.5.1           later_1.4.8             bitops_1.0-9           
#   [5] R.oo_1.27.1             polyclip_1.10-7         fastDummies_1.7.6       lifecycle_1.0.5        
#   [9] httr2_1.2.2             Rdpack_2.6.6            globals_0.19.1          lattice_0.22-7         
#  [13] MASS_7.3-65             magrittr_2.0.4          plotly_4.12.0           rmarkdown_2.31         
#  [17] yaml_2.3.12             remotes_2.5.0           httpuv_1.6.17           otel_0.2.0             
#  [21] sctransform_0.4.3       spam_2.11-3             zip_2.3.3               sessioninfo_1.2.3      
#  [25] spatstat.sparse_3.2-0   reticulate_1.46.0       pbapply_1.7-4           minqa_1.2.8            
#  [29] RColorBrewer_1.1-3      abind_1.4-8             Rtsne_0.17              R.utils_2.13.0         
#  [33] rappdirs_0.3.4          GenomeInfoDbData_1.2.14 ggrepel_0.9.8           irlba_2.3.7            
#  [37] listenv_0.10.1          spatstat.utils_3.2-3    goftest_1.2-3           RSpectra_0.16-2        
#  [41] spatstat.random_3.4-5   fitdistrplus_1.2-6      parallelly_1.46.1       codetools_0.2-20       
#  [45] DelayedArray_0.36.0     tidyselect_1.2.1        UCSC.utils_1.4.0        farver_2.1.2           
#  [49] viridis_0.6.5           lme4_2.0-1              spatstat.explore_3.8-0  jsonlite_2.0.0         
#  [53] progressr_0.19.0        ggridges_0.5.7          survival_3.8-3          tools_4.5.1            
#  [57] ica_1.0-3               glue_1.8.0              gridExtra_2.3           SparseArray_1.10.8     
#  [61] xfun_0.55               withr_3.0.2             BiocManager_1.30.27     fastmap_1.2.0          
#  [65] boot_1.3-31             digest_0.6.39           rsvd_1.0.5              timechange_0.4.0       
#  [69] R6_2.6.1                mime_0.13               scattermore_1.2         tensor_1.5.1           
#  [73] dichromat_2.0-0.1       spatstat.data_3.1-9     R.methodsS3_1.8.2       httr_1.4.8             
#  [77] htmlwidgets_1.6.4       S4Arrays_1.10.1         uwot_0.2.4              pkgconfig_2.0.3        
#  [81] gtable_0.3.6            lmtest_0.9-40           S7_0.2.1                XVector_0.50.0         
#  [85] htmltools_0.5.9         dotCall64_1.2           png_0.1-9               spatstat.univar_3.2-0  
#  [89] reformulas_0.4.4        rstudioapi_0.18.0       tzdb_0.5.0              reshape2_1.4.5         
#  [93] nlme_3.1-168            nloptr_2.2.1            zoo_1.8-15              KernSmooth_2.23-26     
#  [97] parallel_4.5.1          miniUI_0.1.2            pillar_1.11.1           grid_4.5.1             
# [101] vctrs_0.7.2             RANN_2.6.2              promises_1.5.0          xtable_1.8-8           
# [105] cluster_2.1.8.2         evaluate_1.0.5          cli_3.6.6               compiler_4.5.1         
# [109] rlang_1.1.7             future.apply_1.20.2     plyr_1.8.9              stringi_1.8.7          
# [113] viridisLite_0.4.3       deldir_2.0-4            lazyeval_0.2.3          spatstat.geom_3.7-3    
# [117] RcppHNSW_0.6.0          hms_1.1.4               patchwork_1.3.2         future_1.70.0          
# [121] shiny_1.13.0            rbibutils_2.4.1         ROCR_1.0-12             igraph_2.2.3           
