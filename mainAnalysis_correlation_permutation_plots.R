library(Seurat)
library(ggplot2)
library(dplyr)
library(dittoSeq)
library(ComplexHeatmap)
library(circlize)
library(arrow)  
library(Matrix)
library(tidyr)
library(tibble)
library(RColorBrewer)

seurat_obj <- LoadXenium("/data/lowe/reyesj3/pdac_loh/Xenium/rawdata/20231117__193818__2023-11-17_AF-2693_SNV/output-XETG00065__0011655__Region_3__20231117__193845")

## adding the qv score and counts from nuclei 
transcripts <- read.csv("transcript.csv")
qv_summary <- transcripts %>%
  filter(cell_id != "UNASSIGNED") %>%  
  group_by(cell_id) %>%              
  summarise(mean_qv = mean(qv, na.rm = TRUE)) 

filtered_transcripts <- subset(transcripts, overlaps_nucleus == 1 & cell_id != "UNASSIGNED")
filtered_counts <- aggregate(
  x = list(count = filtered_transcripts$feature_name),
  by = list(gene = filtered_transcripts$feature_name, cell = filtered_transcripts$cell_id),
  FUN = length
)

metadata <- seurat_obj@meta.data
cell_ids <- rownames(metadata)
filtered_counts <- filtered_counts[filtered_counts$cell %in% rownames(metadata), ]
valid_cell_ids <- unique(filtered_counts$cell)
seurat_obj <- subset(seurat_obj, cells = valid_cell_ids)
genes <- unique(filtered_counts$gene)
cells <- unique(filtered_counts$cell)
filtered_sparse_matrix <- sparseMatrix(
  i = match(filtered_counts$gene, genes),
  j = match(filtered_counts$cell, cells),
  x = filtered_counts$count,
  dimnames = list(genes, cells)
)
seurat_obj[["Nuclei"]] <- CreateAssayObject(counts = filtered_sparse_matrix)
seurat_obj <- AddMetaData(seurat_obj, metadata = qv_summary)

## analysis

VlnPlot(seurat_obj, features = c("mean_qv", "nCount_Xenium","nFeature_Xenium"), group.by="UPAR_label")
cancer <- subset(seurat_obj, subset = mean_qv > 25 & nCount_Nuclei > 10 & nFeature_Nuclei > 3)

counts<- cancer@assays$Nuclei$counts %>% as.matrix() %>% t()
scaling_factor <- mean(cancer@meta.data$nCount_Nuclei)
norm <- Matrix::Diagonal(x = scaling_factor/cancer@meta.data$nCount_Nuclei,names=colnames(counts)) %*% counts
rownames(norm) <- rownames(counts)
norm_sparse <- as(norm, "dgCMatrix")
norm_sparse_t <- norm_sparse  %>% t()
cancer[["RNA_norm"]] <- CreateAssayObject(data = norm_sparse_t)
DefaultAssay(cancer) <- "RNA_norm"

norm_data <- cancer@assays$RNA_norm@data
norm_data_log <- norm_data
norm_data_log@x <- log1p(norm_data@x) 
cancer[["RNA_norm_log"]] <- CreateAssayObject(data = norm_data_log)

cancer <- FindVariableFeatures(cancer)
cancer <- ScaleData(cancer)
cancer <- RunPCA(cancer, assay = "RNA_norm", reduction.name = "pca.xenium")
cancer <- FindNeighbors(cancer, assay = "RNA_norm", reduction = "pca.xenium", dims = 1:30)
cancer <- FindClusters(cancer, cluster.name = "seurat_clusters", resolution = 0.3)
cancer <- RunUMAP(cancer, assay= "RNA_norm", reduction = "pca.xenium", reduction.name = "umap.xenium",
 return.model = T, dims = 1:30)

# and then fibroblasts, immune cells, PDAC was subsetted normalized scaled and labeled 
# individually according to their FindAllMarkers results carried back to the main UMAP.

#this is run for the dotplot + added canonical markers 
markers <- FindAllMarkers(cancer, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, assay="RNA_norm_log") 
significant <- markers %>% filter(p_val_adj < 0.05) %>% arrange(cluster)

#Figure 1g
DefaultBoundary(cancer[["fov"]]) <- "segmentation"
colors <- c("Basal_PDAC" = "red", 
            "Classical_PDAC" = "orange",
            "uparCAF" = "blue",
            "UPAR+MacroMono"= "green")
ImageDimPlot(cancer, coord.fixed=T, border.size=NA, axes=F, cols=colors)

## correlation heatmap 
## after running the distance script you will have a all_final_distance_results.csv file 
#per cell you have a cellular neighborhood and get the percentage of cell types in each neighborhood

raw_distances <- read.csv("/Users/elif/Desktop/pHD/1st year/Rotations/Lowe Lab/BB209019_tumor_untreated/analysis/only_nuclei_analysis/correlation_heatmap/all_final_distance_results.csv")
all_final_distance_results <- data.table(raw_distances)
cell_type_columns <- c("Tcells", "NaiveB", "iCAF", "ActiveB", "Macrophage", 
                       "Tregs", "Endothelial", "VascularSM", "Mast", "Mixed_upar.", 
                       "Glial", "uPAR+Cycling", "Acinar_fibro", "myCAF", "Islet", 
                       "CircularSM", "Basal_PDAC", "uparCAF", "Glial_fibro", 
                       "Classical_PDAC", "UPAR.MacroMono", "Acinar")
filtered_cell_type_columns <- c(
   "UPAR.MacroMono","Mixed_upar.","uPAR+Cycling",
   "uparCAF", "Basal_PDAC","Classical_PDAC"
)

final_result <- all_final_distance_results %>% filter(distance == 250)

## combine Cycling and Classical PDAC
final_result$target_cell_type <- gsub("Cycling_PDAC", "Classical_PDAC",  final_result$target_cell_type)
final_result$Classical_PDAC <- final_result$Classical_PDAC + final_result$Cycling_PDAC
final_result$Cycling_PDAC <- NULL  


cell_percentage_data <- final_result[, ..cell_type_columns]

correlation_matrix <- cor(cell_percentage_data, method = "pearson", use = "complete.obs")
hc <- hclust(as.dist(1 - correlation_matrix)) #hierarchical ordering for the names
ordered_names <- rownames(correlation_matrix)[hc$order]

correlation_matrix <- correlation_matrix[ordered_names, ordered_names]
correlation_data <- reshape2::melt(correlation_matrix)
correlation_data$Var1 <- factor(correlation_data$Var1, levels = ordered_names)
correlation_data$Var2 <- factor(correlation_data$Var2, levels = ordered_names)

correlation_data <- correlation_data %>%
  mutate(
    Var1 = recode(Var1, "Mixed_upar." = "uPAR+Mixed", "uparCAF" = "uPAR+CAF", "UPAR.MacroMono" = "uPAR+Myeloid","uPAR+Cycling" = "uPAR+Cycling"),
    Var2 = recode(Var2, "Mixed_upar." = "uPAR+Mixed", "uparCAF" = "uPAR+CAF", "UPAR.MacroMono" = "uPAR+Myeloid","uPAR+Cycling" = "uPAR+Cycling")
  )

## run this before plotting for the small heatmap ##
excluded_types <- c("Tcells", "NaiveB", "iCAF", "ActiveB", "Macrophage", 
                       "Tregs", "Endothelial", "VascularSM", "Mast",
                       "Glial", "Acinar_fibro", "myCAF", "Islet", 
                       "Glial_fibro", "Acinar","CircularSM") 

correlation_data<- correlation_data[
!(correlation_data$Var1 %in% excluded_types | correlation_data$Var2 %in% excluded_types), 
]
## end ## 

heatmap_plot <- ggplot(correlation_data, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(midpoint = 0, low = "blue", high = "red", mid = "white", name = "Correlation") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=9),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
## permutation analysis 

final_result <- all_final_distance_results[distance == 250] #select the distance of interest
final_result$target_cell_type <- gsub("Cycling_PDAC", "Classical_PDAC",  final_result$target_cell_type) ## combining classical with cycling since same type
final_result$Classical_PDAC <- final_result$Classical_PDAC + final_result$Cycling_PDAC
final_result$Cycling_PDAC <- NULL  #combine cycling PDAC with Classical PDAC
basal_pdac_cells <- final_result[target_cell_type == "Basal_PDAC"] #select the cell type of interest. 
classical_pdac_cells <- final_result[target_cell_type == "Classical_PDAC"] #select the cell type of interest.

num_iterations <- 10000
max_cell_types <- c()

set.seed(123)  
for (i in 1:num_iterations) {

  sampled_cells <- basal_pdac_cells[sample(.N, size = 10)] #pick 10 basal PDAC cells
  
  unwanted_columns <- c("X", "Basal_PDAC", "cell_id", "target_cell_type", "distance") #remove cell itself
  sampled_data <- sampled_cells[, !names(sampled_cells) %in% unwanted_columns, with = FALSE]

  average_percentages <- colMeans(sampled_data, na.rm = TRUE) # get the mean of per cell type percentage
  max_cell_type <- names(which.max(average_percentages)) # get the top result, the most common cell type within averages
  max_cell_types <- c(max_cell_types,max_cell_type)
} #repeat this 10000 times

frequency_table <- as.data.table(table(max_cell_types))
setnames(frequency_table, c("max_cell_types", "frequency"))

gg_basal <- ggplot(frequency_table, aes(x = max_cell_types, y = frequency, fill = max_cell_types)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +  # Apply custom colors
  labs(
    title = "",
    x = "Cell Type",
    y = "Frequency",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x = element_blank()
  ) + NoLegend()

## classical
num_iterations <- 10000
max_cell_types <- c()

set.seed(123)  
for (i in 1:num_iterations) {

  sampled_cells <- classical_pdac_cells[sample(.N, size = 10)] #pick 10 basal PDAC cells
  
  unwanted_columns <- c("X", "Classical_PDAC", "cell_id", "target_cell_type", "distance") #remove cell itself
  sampled_data <- sampled_cells[, !names(sampled_cells) %in% unwanted_columns, with = FALSE]

  average_percentages <- colMeans(sampled_data, na.rm = TRUE) # get the mean of per cell type percentage
  max_cell_type <- names(which.max(average_percentages)) # get the top result, the most common cell type within averages
  max_cell_types <- c(max_cell_types,max_cell_type)
} #repeat this 10000 times

frequency_table <- as.data.table(table(max_cell_types))
setnames(frequency_table, c("max_cell_types", "frequency"))

gg_classical <- ggplot(frequency_table, aes(x = max_cell_types, y = frequency, fill = max_cell_types)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = cell_type_colors) +  # Apply custom colors
  labs(
    title = "",
    x = "Cell Type",
    y = "Frequency",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.text.x = element_blank()
  ) + NoLegend()


cell_type_colors <- c(
"Tcells" = "#BF9000",        
    "NaiveB" = "#33FF57",         
    "Fibroblast" = "#00A2F0",     
    "ActiveB" = "#FF33A8",        
    "Macrophage" = "#FFC300",    
    "Tregs" = "#DA3D22",        
    "Endothelial" = "#AFABAB",   
    "VascularSM" = "#FF9CEF",     
    "Mast" = "#002060",          
    "Mixed_upar." = "#0096B1",    
    "Glial" = "#2ECC71",         
    "uPAR+Cycling" = "#1300FF",  
    "Islet" = "#9B59B6",         
    "CircularSM" = "#FF9CEF",     
    "PDAC" = "#E67E22",           
    "UPAR.MacroMono" = "#DA3B91",
    "Acinar" = "#C0392B",      
    "myCAF" = "#90BF72",
    "iCAF" = "#905F72",
    "uparCAF" = "#00A2F0",
    "Acinar_fibro" = "#FFFFFF",
    "Glial_fibro" = "#FFFFFF",
    "Classical_PDAC" = "#E67E22", 
    "Basal_PDAC" = "#00FFF7", 
    "Cycling_PDAC" = "#FFFF00"
  )
