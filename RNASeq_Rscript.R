#### Load libraries ####
library(edgeR)
library(limma)
library(fgsea)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(DESeq2)
library(dplyr)
library(AnnotationDbi)
library(ggrepel)
library(fgsea)
library(msigdbr)

#### Load data ####
counts_path <- paste0("path/to/file/gene_symbol_counts_matrix.csv")
metadata_path <- paste0("path/to/file/metadata.csv")
counts <- read.csv(counts_path, row.names = 1)
metadata <- read.csv(metadata_path)

#### Create DESeq2 object for size factor estimation ####
# Remove 0 values
counts <- counts %>%
  filter_all(all_vars(. != 0))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Type)
dds <- estimateSizeFactors(dds)

# Apply Variance Stabilizing Transformation
vst_data <- vst(dds, blind = FALSE)
vst_counts <- assay(vst_data)
write.csv(vst_counts, file = "path/to/file/vst_counts.csv")
vst_path<-paste0("path/to/file/vst_counts.csv")
vst_counts<-read.csv(vst_path)

# Calculate FPKM
fpkm <- counts(dds, normalized = TRUE) / (sizeFactors(dds) * 1e-6)

# Write FPKM values to an Excel file
write.csv(fpkm, file = paste0("path/to/file/fpkm_values.csv"))

#### Create sample-to-sample correlation plot ####
sample_cor <- cor(vst_counts)
pheatmap(sample_cor, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         main = "Sample-to-Sample Correlation Matrix",
         fontsize = 12,
         fontsize_row = 12,
         fontsize_col = 12, 
         fontsize_title = 16)

#### Create heatmap ####
# Calculate variances and get top 50 genes by variance
variances <- apply(vst_counts, 1, var)
ordered_indices <- order(variances, decreasing = TRUE)
top_genes <- head(ordered_indices, 50)

# Create a hierarchical heatmap for the top 50 genes
pheatmap(log2(vst_counts[top_genes, ] + 1), 
         scale = "row", 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean",
         fontsize = 14,
         fontsize_row = 12,
         fontsize_col = 12)

#### PCA plot ####
# Perform PCA
pcaData <- plotPCA(vst_data, intgroup = "Type", returnData = TRUE)

# Calculate percentage of variance
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create PCA plot
ggplot(pcaData, aes(PC1, PC2, color = Type)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()+
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.text = element_text(size = 14), 
    axis.title = element_text(size = 16)
  )

#### Prepare DGE object ####
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~ Type, data = metadata)

# Estimate dispersion and test for differential expression
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
res <- glmQLFTest(fit, coef = 2) # Adjust 'coef' based on your condition placement in design

#### Get top genes for DE analysis ####

# Extract significant genes
sig_genes <- topTags(res, n = Inf)$table
sig_genes <- sig_genes[sig_genes$FDR < 0.05 & sig_genes$PValue < 0.05, ]
sig_gene_fc$id <- row.names(sig_gene_fc)

#Sort by p value
sorted_order = order(sig_genes[,4], decreasing=FALSE)
sig_genes = sig_genes[sorted_order,]

# Make columns for -log10p and direction
sig_genes$gene = rownames(sig_genes)
sig_genes$mlog10p = -log10(sig_genes$PValue)
sig_genes$Direction <- ifelse(sig_genes$logFC > 0, "Upregulated", "Downregulated")

# Write log2 FC normalized values to an Excel file
write.csv(sig_genes, file = paste0("path/to/file/sig_genes.csv"))
sig_path <- paste0("path/to/file/sig_genes.csv")
sig_genes <- read.csv(sig_path, row.names=1)

# upregulated & downregulated genes
sig_up = subset(sig_genes, PValue < 0.05 & logFC > 0)
sig_down = subset(sig_genes, PValue < 0.05 & logFC < 0)
write.csv(sig_up, file = paste0("path/to/file/sig_up.csv"))
write.csv(sig_down, file = paste0("path/to/file/sig_down.csv"))
sig_up_top5 = sig_up[1:5,]
sig_down_top5=sig_down[1:5,]

##### Top 5 upregulated and 5 downregulated genes #####

sig_genes[1:10,]
top10 <- c("GENES")
subset_gene_data <- sig_genes %>% filter(gene %in% top10)
write.csv(subset_gene_data, file = paste0("path/to/file/top_10_genes.csv"))

##### Genes that encode for protein hubs #####

# Get expression level data for genes that encode for protein hubs found through STRING and Cytoscape
protein_genes <- c("GENES")
protein_gene_data <- sig_genes %>% filter(gene %in% protein_genes)
write.csv(protein_gene_data, file = paste0("path/to/file/top10proteindata.csv"))

#### Volcano plot ####
ggplot(sig_genes, aes(x = logFC, y = mlog10p, colour = Direction)) +
  geom_point(alpha = 0.5) +
  xlim(c(-5, 5)) +
  scale_color_manual(values = c("Downregulated" = "cyan4", "Upregulated" = "deeppink4")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16), 
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)   
  ) +
  geom_label_repel(data = sig_up_top5, aes(label = gene), show.legend = FALSE) +
  geom_label_repel(data = sig_down_top5, aes(label = gene), show.legend = FALSE)
ggsave(filename = paste0("path/to/file/volcano_plot.png"),
       width = 10, height = 8, units = "in", dpi = 300)

#### GO Pathways ####

# If gene names are used instead of ensembl ids, use this code to transform
ensembl_ids <- mapIds(org.Hs.eg.db, 
                      keys = sig_gene_ids, 
                      column = "ENSEMBL", 
                      keytype = "SYMBOL", 
                      multiVals = "first")

# Convert ENSEMBL IDs to Entrez IDs as required by clusterProfiler
entrez_ids <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Some genes might not map and will be NA; filter them out
entrez_ids <- entrez_ids[!is.na(entrez_ids$ENTREZID), ]

# Perform GO enrichment analysis
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",  # Biological Process. Change to "CC" or "MF" as needed
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05)

# Visualize top GO terms
barplot(go_results, showCategory = 10)

# Write GO enrichment results to a file
write.csv(as.data.frame(go_results), file = paste0("path/to/file/go_enrichment_results.csv"))

#### GSEA ####
ranked_genes <- sig_genes$logFC
names(ranked_genes) <- rownames(sig_genes)

# Sort the ranked genes in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)
# Load Hallmark gene sets
msigdbr_species <- "Homo sapiens"
msigdbr_category <- "H"  # Hallmark gene sets
msigdbr_sets <- msigdbr(species = msigdbr_species, category = msigdbr_category)

# Convert to list format required by fgsea
pathways <- split(msigdbr_sets$gene_symbol, msigdbr_sets$gs_name)
fgsea_results <- fgsea(pathways = pathways, stats = ranked_genes, nperm = 1000)

# View results
fgsea_results <- fgsea_results[order(fgsea_results$padj), ]
print(fgsea_results)
top_pathways_up <- fgsea_results[ES > 0,][order(padj),][1:5,]
top_pathways_down <- fgsea_results[ES < 0,][order(padj),][1:5,]
top_pathways <- rbind(top_pathways_up, top_pathways_down)

# Plot pathways
ggplot(top_pathways, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES > 0)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Top GSEA Pathways") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, face = "bold"), 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


