# loading the libraries
library(dplyr)
library(tidyr)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(PoiClaClu)
library(RColorBrewer)
library(GenomicRanges)
library(clusterProfiler)
# load the data
rse <- readRDS("data/EwS.rds")

# some data preprocessing
rse$condition <- ifelse(rse$condition=="shCTR","control","knockout")
fixed_rownames <- sapply(strsplit(rownames(rse),".",fixed = T), function(x) x[1])
rownames(rse) <- fixed_rownames

# creating DESeqdataset object

dds <- DESeqDataSet(rse, design = ~condition)
dds <- dds[rowSums(counts(dds)) >=10,] #filter low count genes

dds <- estimateSizeFactors(dds)

# getting the normalised counts
normalized_counts <- counts(dds, normalized = TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

# QC checks
## transform counts for the sample variance viz.
rld <- rlog(dds, blind= TRUE)
# 1. PCA plot
pc <-plotPCA(rld, intgroup = "condition")
ggsave("results/PCA plot.png",plot = pc)

# 2. Hierarchical clustering
## creating matrix from the transformed rld object for visualisation

rld_matrix <- assay(rld)
rld_cor <- cor(rld_matrix)
## creating annotation df 
annotation <- data.frame(condition = rse$condition)
rownames(annotation)<-colnames(rld_matrix)

pheatmap(rld_cor,annotation_col = annotation)

### perform differential gene analysis
dds <- DESeq(dds)

#dispersion plot
disp_plot <-plotDispEsts(dds)

# obtain the results
res <-results(dds)

#MA plot
ma<-plotMA(res,ylim = c(-4,4))


res_df <- data.frame(res)
# histogram of p values
p_hist <-ggplot(res_df, aes(x=pvalue))+
  geom_histogram(bins = 100)

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
#subset significant genes
sig_res <- filter(res_df, padj<padj.cutoff& abs(log2FoldChange)> lfc.cutoff)
# Adding gene annotation to the table with significant genes
anno <- AnnotationDbi::select(org.Hs.eg.db,rownames(sig_res),
                              columns = c("ENSEMBL","ENTREZID","SYMBOL","GENENAME"),
                              keytype = "ENSEMBL")
results <- cbind(ENSEMBL = rownames(sig_res),sig_res)
anno_results <- left_join(as.data.frame(results),anno)
# saving it in the form of table
write.csv(anno_results, file = "data/Ewing_Sarcoma_sigGenes.csv", sep = ",", col.names = TRUE,row.names = TRUE)
###visualisation


###Volcano Plot###
volcano_plot <-EnhancedVolcano(anno_results, x = "log2FoldChange", y = "padj", lab = anno_results$SYMBOL,
                title = "DESeq2 results",
                subtitle = "Differential expression",
                border = "full")

### Set a color palette
heat_colors <- brewer.pal(8, "YlOrRd")
DE_genes<-normalized_counts[anno_results$ENSEMBL,]

DE_heatmap<-pheatmap(DE_genes,
         annotation_col = annotation,
         color=heat_colors,
         show_rownames = F,
         cluster_rows = T,
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20,
         main = "differentially expressed genes")
ggsave("results/heatmap of DEgenes.png",plot=DE_heatmap)

# Subset the data frame to remove rows with NA values
anno_results <- anno_results[rowSums(is.na(anno_results)) == 0, ]
anno_results <- anno_results[which(duplicated(anno_results$ENTREZID)==F),]

# extract the fold changes
foldChanges <- anno_results$log2FoldChange
names(foldChanges)<- anno_results$ENTREZID
## Sort fold changes in decreasing order
foldChanges <- sort(foldChanges, decreasing = TRUE)
foldChanges <- setNames(as.data.frame(foldChanges), "log2FoldChange") %>%
  cbind(ENTREZID = rownames(.))

# Perform KEGG pathway enrichment analysis
gene_list <- foldChanges$ENTREZID
kegg_enrich<-enrichKEGG(
  gene_list,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500)
kegg_res <- kegg_enrich@result %>%
  as.data.frame()%>%
  arrange(pvalue)

#Writing results to CSV
#write.csv(kegg_res, "results/KEGGEnrichmentAnalysis.csv")
#Figure to summarize Results (Figure to Summarize Enrichment Analysis?)



# Plot the top 10 enriched pathways
top_pathways <- head(kegg_res, 10)

# Create the bar plot
bar_plot <- ggplot(top_pathways, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "red") +
  labs(title = "Top 10 Enriched KEGG Pathways",
       x = "Pathways",
       y = "-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("results/KEGG enrichment analysis barplot.png",plot=bar_plot)
 