
library(Seurat)
library(dplyr)
library(ggplot2)
library(UCell)
library(ggpubr)
library(ggrepel)

cancer <- readRDS("paper_object.rds")

## adding senescence scores to the cells 

# we have created this list by overlapping the scores existing in our lab and literature with the xenium panel genes 
senescence_signatures <- read.csv("correct_list_senescent_markers_overlapwXenium.csv")
senescence_signatures <- senescence_signatures[,c("SeneUP_ovlp_TCGAvsGTExUP.minpct0.4", "SENESCopedia_Core245", "Tasdemir_Senescence_UP")]
senescence_signatures <- senescence_signatures[
  apply(senescence_signatures != "" & !is.na(senescence_signatures), 1, any),
]
data <- cancer@assays$RNA_norm@data
scores <- ScoreSignatures_UCell(data, features=senescence_signatures, maxRank=53)
cancer <- AddMetaData(cancer, scores) # adding UCell scores back to the seurat object 

## plotting the scores 

SpatialColors <- colorRampPalette(colors = rev(brewer.pal(n = 11, name = "Spectral")))

module_scores <- c(
  "SeneUP_ovlp_TCGAvsGTExUP.minpct0.4_UCell",
  "SENESCopedia_Core245_UCell",
  "Tasdemir_Senescence_UP_UCell"
)

DefaultBoundary(cancer[["fov"]]) <- "segmentation" ## get the segmentation masks for better visualization 
my_palette <- colorRampPalette(brewer.pal(11, "Spectral"))(100)

##plotting the scores on the tissue 
for (score in module_scores) { 

  # Get 1st and 99th percentile for color scaling
  normalized_quantiles <- quantile(cancer@meta.data[[score]], probs = c(0.01, 0.99), na.rm = TRUE)

  pdf(paste0(score, "_Score.pdf"), height = 30, width = 30)

  gg <- ImageFeaturePlot(
    object = cancer,
    features = score,
    border.size = NA
  ) & scale_fill_gradientn(
    limits = c(0, normalized_quantiles[[2]]),
    colours = rev(my_palette)
  )

  plot(gg)
  dev.off()
}

## plotting individual myeloid and fibroblast boxplot in the figure

myeloids <- subset(cancer, idents=c("Macrophage","UPAR+MacroMono"))
fibroblasts <- subset(cancer, idents=c("myCAF","uPAR+CAF","iCAF"))

## myeloid box plot socre 

pla_expr <- FetchData(myeloids, vars = "PLAUR", slot = "data")[,1]
myeloids$uPAR_binary <- ifelse(pla_expr > 0, "uPAR-positive", "uPAR-negative") # set the uPAR+ cells 

dittoBarPlot(myeloids, var = "uPAR_binary", group.by = "PaperMainClusters")

ucell_cols <- grep("_UCell$", colnames(myeloids@meta.data), value = TRUE)
meta_long <- myeloids@meta.data %>%
  select(all_of(ucell_cols), uPAR_binary) %>%
  pivot_longer(cols = all_of(ucell_cols), names_to = "Signature", values_to = "Score")

output_dir <- "myeloid_posNeg_boxplot_PDAC"
dir.create(output_dir, showWarnings = FALSE)

custom_colors <- c("uPAR-positive" = "#4DAF4A",  # green
                   "uPAR-negative" = "#E41A1C")  # red-ish

ucell_cols <- grep("_UCell$", colnames(cancer@meta.data), value = TRUE)
for (signature in ucell_cols) {
  plot_data <- myeloids@meta.data %>%
    select(Score = all_of(signature), uPAR_binary)

  p <- ggplot(plot_data, aes(x = uPAR_binary, y = Score, fill = uPAR_binary)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.3, color = "black") +
    geom_jitter(width = 0.1, size = 0.3, alpha = 0.3, color = "grey") +
    scale_fill_manual(values = custom_colors) +
    stat_compare_means(
      method = "wilcox.test",  
      label = "p.signif",
      comparisons = list(c("uPAR-positive", "uPAR-negative")),
      label.y = max(plot_data$Score, na.rm = TRUE) + 0.03
    ) +
    theme_minimal(base_size = 10) +
    labs(
      title = signature,
      x = NULL,
      y = "UCell Score"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid = element_blank()
    )

  ggsave(
    filename = file.path(output_dir, paste0(signature, ".pdf")),
    plot = p,
    width = 5,
    height = 5
  )
}

## fibroblast boxplot score 

pla_expr <- FetchData(fibroblasts, vars = "PLAUR", slot = "data")[,1]
fibroblasts$uPAR_binary <- ifelse(pla_expr > 0, "uPAR-positive", "uPAR-negative")
dittoBarPlot(fibroblasts, var = "uPAR_binary", group.by = "PaperMainClusters")

ucell_cols <- grep("_UCell$", colnames(fibroblasts@meta.data), value = TRUE)

meta_long <- fibroblasts@meta.data %>%
  select(all_of(ucell_cols), uPAR_binary) %>%
  pivot_longer(cols = all_of(ucell_cols), names_to = "Signature", values_to = "Score")

output_dir <- "fibrob_posNeg_boxplot_PDAC"
dir.create(output_dir, showWarnings = FALSE)

custom_colors <- c("uPAR-positive" = "#4DAF4A",  # green
                   "uPAR-negative" = "#E41A1C")  # red-ish


for (signature in ucell_cols) {
  plot_data <- fibroblasts@meta.data %>%
    select(Score = all_of(signature), uPAR_binary)

  p <- ggplot(plot_data, aes(x = uPAR_binary, y = Score, fill = uPAR_binary)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.3, color = "black") +
    geom_jitter(width = 0.1, size = 0.3, alpha = 0.3, color = "grey") +
    scale_fill_manual(values = custom_colors) +
    stat_compare_means(
      method = "wilcox.test",  
      label = "p.signif",
      comparisons = list(c("uPAR-positive", "uPAR-negative")),
      label.y = max(plot_data$Score, na.rm = TRUE) + 0.03
    ) +
    theme_minimal(base_size = 10) +
    labs(
      title = signature,
      x = NULL,
      y = "UCell Score"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid = element_blank()
    )

  ggsave(
    filename = file.path(output_dir, paste0(signature, ".pdf")),
    plot = p,
    width = 5,
    height = 5
  )
}

## correlation plot between uPAR expression and senescence score 

pla_expr <- FetchData(pdac_data, vars = "PLAUR", slot = "data")[, 1]
pdac_data$PLAUR_expr <- pla_expr # add as a metadata column 

# Compute mean senescence and PLAUR by cell types
df_celltype_means <- pdac_data@meta.data %>%
  group_by(PaperSubclusters) %>%
  summarise(
    mean_senescence = mean(SeneUP_ovlp_TCGAvsGTExUP.minpct0.4_UCell, na.rm = TRUE),
    mean_PLAUR = mean(PLAUR_expr, na.rm = TRUE),
    n_cells = n()
  )


# Calculate correlation
cor_test <- cor.test(df_celltype_means$mean_senescence, df_celltype_means$mean_PLAUR)
r_value <- round(cor_test$estimate, 2)
p_value <- signif(cor_test$p.value, 2)
cor_label <- paste0("R = ", r_value, ", p = ", p_value)

ggplot(df_celltype_means, aes(x = mean_senescence, y = mean_PLAUR)) +
  geom_point(color = "#4DAF4A", size = 2.5) +
  geom_text_repel(aes(label = PaperSubclusters), size = 4, max.overlaps = Inf, box.padding = 0.3) +
  annotate("text", 
           x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5, 
           label = cor_label, size = 3.2, fontface = "italic") +
  theme_minimal(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
  ) +
  labs(
    title = "Senescence vs PLAUR Expression (by Cell Type)",
    x = "Mean SeneUP_ovlp_TCGAvsGTExUP.minpct0.4_UCell Score",
    y = "Mean PLAUR Expression"
  ) 

## senescence score dotplots 

dotplot_genes <- c(senescence_signatures$SeneUP_ovlp_TCGAvsGTExUP.minpct0.4,senescence_signatures$SENESCopedia_Core245, senescence_signatures$Tasdemir_Senescence_UP)
gg1 <- DotPlot(cancer, features = unique(c(dotplot_genes, "CDKN2A","CDKN1A")), 
                    assay = "RNA_norm",
                    scale= F, cluster=T) +
            ggtitle(paste("senescence score genes")) +  
            xlab("Marker Genes") +
            ylab("") +
            scale_color_gradient2(low = "#01A9DB", mid = "white", high = "#8A084B", midpoint = 0, name = "Expression") +
            RotatedAxis() +
            theme(axis.text.y = element_text(size = rel(1), angle = 0, hjust = 1, vjust = 0.5, face = "bold.italic"),
                  axis.text.x = element_text(size = rel(1)),
                  legend.position = "right", legend.box = "vertical", legend.direction = "vertical")

my_palette <- colorRampPalette(brewer.pal(11, "Spectral"), direction=-1)(100)
gg2 <- DotPlot(cancer, features = c("SeneUP_ovlp_TCGAvsGTExUP.minpct0.4_UCell","SENESCopedia_Core245_UCell","Tasdemir_Senescence_UP_UCell"), 
                    assay = "RNA_norm",
                    scale= F, cluster=T) +
            ggtitle(paste("senescence score genes")) +  
            xlab("Marker Genes") +
            ylab("") +
            scale_color_gradientn(
            colours = rev(my_palette),
            name = "Score") +
            RotatedAxis() +
            theme(axis.text.y = element_text(size = rel(1), angle = 0, hjust = 1, vjust = 0.5, face = "bold.italic"),
                  axis.text.x = element_text(size = rel(1)),
                  legend.position = "right", legend.box = "vertical", legend.direction = "vertical") +coord_flip()


## senescence score over fibroblast types

caf_cols <- c(
  "myCAF"    = "#90BF72",
  "iCAF"     = "#905F72",
  "uPAR+CAF" = "#00A2F0"
)

p <- VlnPlot(
  cancer,features = c("SeneUP_ovlp_TCGAvsGTExUP.minpct0.4_UCell","SENESCopedia_Core245_UCell","Tasdemir_Senescence_UP_UCell"),
  idents   = c("myCAF", "iCAF", "uPAR+CAF"),
  group.by = "PaperSubclusters",
  pt.size  = 0.001
) &
  scale_fill_manual(values = caf_cols) &
  geom_boxplot(
    width = 0.1,          # thinner than violin
    outlier.shape = NA,    # hide outliers
    color = "white",
    alpha = 0.1
  ) &
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
