cuttag <- read.csv("chip_seeker/Nrf2.csv")
deep <- read.csv("chip_seeker/ids.csv")
cuttag$Var1 <- toupper(cuttag$Var1)
a<- cuttag[cuttag$Var1 %in% deep$SYMBOL,]
a$Var1 <- tolower(a$Var1)
a$Var1 <- stringr::str_to_title(a$Var1)
##############
library(DESeq2)
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(forcats)
library(org.Mm.eg.db)
library(pheatmap)
library(ggrepel)
library(ggplot2)
library(igraph)
library(enrichplot)
BiocManager::install("clusterProfiler",force = TRUE)
ids <- bitr(a$Var1,'SYMBOL','ENTREZID','org.Mm.eg.db')
ego_ALL <- enrichGO(gene = ids$ENTREZID,
                    OrgDb=org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE) 
write.csv(ego_ALL,"ego_nrf2_brain.csv")

p5 <- ggplot(ego_ALL@result[ego_ALL@result$Description %in% target,], 
             aes(x = Count, y = Description, color = -log(pvalue), size = Count)) + 
  geom_point() + 
  scale_color_gradientn(
    name = "pvalue", 
    values = seq(0, 1, 0.5), 
    colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  scale_x_continuous(limits = c(10, 120)) +
  geom_point(shape = 21, color = "black", stroke = 1, alpha = 1) + 
  scale_size_continuous(name = "Count", range = c(4, 8), breaks = c(10, 40, 80, 120, 160)) +  
  theme_minimal() +  
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    axis.line = element_line(colour = "black"),  
    panel.background = element_rect(fill = "white"), 
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),
    axis.text = element_text(family = "Arial", size = 12), 
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank() 
  )

print(p5)

kk <- enrichKEGG(gene = ids$ENTREZID,keyType = "kegg",organism= "mmu", qvalueCutoff = 1, pvalueCutoff=1)
kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
kk1 <- as.data.frame(kk)
write.csv(kk,"kk2.csv")
#########cuttag
cuttag_up <- read.csv("chip_seeker/cbd_cuttag_up.csv")
cuttag_up <- cuttag_up %>%
  dplyr::filter(abs(distanceToTSS) < 3000)
cuttag_down <- read.csv("chip_seeker/cbd_cuttag_down.csv")
cuttag_down <- cuttag_down %>%
  dplyr::filter(abs(distanceToTSS) < 3000)
all <- rbind(cuttag_down,cuttag_up)


ids <- bitr(all$geneId,'SYMBOL','ENTREZID','org.Mm.eg.db')
ego_ALL <- enrichGO(gene = ids$ENTREZID,
                    OrgDb=org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE) 
write.csv(as.data.frame(ego_ALL),"ego_ALL.CSV")

cuttag_down <- cuttag_down[cuttag_down$geneId %in% a$Var1,]
cuttag_up <-  cuttag_up[cuttag_up$geneId %in% a$Var1,]
##################
nrf2 <- read.csv("../step3_RNA-seq/CBD_DEG.csv")
NRF2 <- nrf2 %>%
  dplyr::filter(P.Value < 0.05)

NRF2_UP <- NRF2%>%
  dplyr::filter( logFC>0)
NRF2_DOWM <- NRF2%>%
  dplyr::filter( logFC<0)


cuttag_down <- cuttag_down[cuttag_down$geneId %in% NRF2_DOWM$Row.names,]
cuttag_up <-  cuttag_up[cuttag_up$geneId %in% NRF2_UP$Row.names,]
all <- rbind(cuttag_down,cuttag_up)


ids <- bitr(cuttag_up$geneId,'SYMBOL','ENTREZID','org.Mm.eg.db')
ego_ALL <- enrichGO(gene = ids$ENTREZID,
                    OrgDb=org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE) 
write.csv(as.data.frame(ego_ALL),"ego_ALL1.CSV")


p5 <- ggplot(ego_ALL@result[ego_ALL@result$Description %in% target,], 
             aes(x = Count, y = Description, color = -log(pvalue), size = Count)) + 
  geom_point() + 
  scale_color_gradientn(
    name = "pvalue", 
    values = seq(0, 1, 0.5), 
    colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 100)) +
  scale_x_continuous(limits = c(0, 4)) +
  geom_point(shape = 21, color = "black", stroke = 1, alpha = 1) + 
  scale_size_continuous(name = "Count", range = c(4, 8), breaks = c(10, 40, 80, 120, 160)) +  
  theme_minimal() +  
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.line = element_line(colour = "black"), 
    panel.background = element_rect(fill = "white"),  
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title.x = element_blank(),  
    axis.title.y = element_blank(),  
    axis.text = element_text(family = "Arial", size = 12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()  
  )
p5
###########

my_colors <- c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25") 

p1 <- pheatmap(data_matrix,  
               color = colorRampPalette(my_colors)(100),
               border_color = "black",  
               scale = "row", 
               cluster_rows =T, 
               cluster_cols = T, 
               legend = TRUE, 
               legend_breaks = c(-1, 0, 1), 
               legend_labels = c("low","","high"),
               show_rownames = FALSE, 
               show_colnames = TRUE, 
               kmeans_k = 250,
               clustering_method = "ward.D",
               annotation_colors = annotation_colors,
               fontsize = 8
)
