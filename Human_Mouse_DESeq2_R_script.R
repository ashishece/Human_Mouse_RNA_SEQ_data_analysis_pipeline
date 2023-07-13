setwd("C:/Users/hp/OneDrive/Desktop")
getwd()
library(DESeq2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)

countData <- read.csv('counts.csv', header = TRUE, sep = ",")
head(countData)
tail(countData,10)

metaData <- read.csv('metadata.csv', header = TRUE, sep = ",")
metaData


dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~dex, tidy=TRUE)
dds$dex
dds$dex <- relevel(dds$dex, ref = "control")

dds <- DESeq(dds)

results <- results(dds)
head(results(dds))

summary(results)

resdata <- merge(as.data.frame(results), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'ENSEMBL'
write.csv(resdata, file ="results-with-normalized.csv")



columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

keys(org.Hs.eg.db, keytype="SYMBOL")[1:10]

select(org.Hs.eg.db, keys=resdata$ENSEMBL,
       keytype = "ENSEMBL",columns=c("SYMBOL","GENENAME")
)

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=resdata$ENSEMBL,
                              columns=c("SYMBOL","GENENAME"),
                              keytype="ENSEMBL")

anno <- AnnotationDbi::select(org.Hs.eg.db,keys=resdata$ENSEMBL,
                              columns=c("ENSEMBL","SYMBOL","GENENAME"),
                              keytype="ENSEMBL") %>% 
  filter(!duplicated(SYMBOL))

dim(anno)


anno <- dplyr::rename(anno, GeneID = SYMBOL)
results_annotated <- left_join(resdata, anno,by="ENSEMBL")

head(results_annotated) 

write.csv(results_annotated, file="Sample_deseq2_normalized_annotated.csv")


vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="dex")


library(EnhancedVolcano)


res <- results(dds, contrast = c('dex', 'control', 'diseased'))
res <- lfcShrink(dds, contrast=c('dex', 'control', 'diseased'), res=res, type= 'normal')


#volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')



EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'diseased versus control',
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)



#Hierarchical clustering (Cluster Dendrogram)
data <- read.csv("counts_dendogram.csv", header=TRUE, sep = ",")
dist<- dist(t(data), method="euclidean")
dist
hc <- hclust(dist)
plot(hc)


dend <- dist %>% scale %>% dist %>% 
  hclust %>% as.dendrogram %>%
  set("branches_k_color", k=3) %>% set("branches_lwd", 1.2) %>%
  set("labels_colors") %>% set("labels_cex", c(.9,1.2)) %>% 
  set("leaves_pch", 19) %>% set("leaves_col", c("blue", "red"))
# plot the dend in usual "base" plotting engine:
plot(dend)


ggd1 <- as.ggdend(dend)
ggplot(ggd1) 
ggplot(ggd1, horiz = TRUE, theme = NULL) 
ggplot(ggd1, theme = theme_minimal()) 
ggplot(ggd1, labels = FALSE) + 
  scale_y_reverse(expand = c(0.2, 0)) +
  coord_polar(theta="x")



#Heatmap
library(pheatmap)
library(tidyverse)
countfile<- read.csv("heatmap.csv", header = TRUE, sep = ",")
head(countfile)

countfile %>%
  dplyr::select(1:10) %>%
  column_to_rownames("Gene") -> heatmap_data
head(heatmap_data)

heatmap_data %>%
  pheatmap()

heatmap_data %>%
  log2() -> heatmap_data_log2
head(heatmap_data_log2)
heatmap_data_log2 %>%
  pheatmap()

heatmap_data_log2 - rowMeans((heatmap_data_log2)) ->
  heatmap_data_meanSubtract
head(heatmap_data_meanSubtract)

heatmap_data_meanSubtract %>%
  pheatmap()

heatmap_data_meanSubtract/rowSds(as.matrix(heatmap_data_log2)) ->
  heatmap_data_zscores
heatmap_data_zscores %>%
  pheatmap()