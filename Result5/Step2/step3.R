library(limma)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(dplyr)
library(GseaVis)
#########
b <- read.csv("output_first_three_columns_with_gene_names.csv")
b$id <- c(1:length(rownames(b)))
b <-b[,c(4,5)]
colnames(b) <- c("ID","id")
library(data.table)
c <-fread("GPL570-55999.txt",)
c <- c[ c$ID%in% b$ID,]
c <- c[,c(1,11)]
c$`Gene Symbol`<- sub("///.*", "", c$`Gene Symbol`)
E  <- merge(b,c,by="ID")
rownames(E) <- E$ID

#######
pre <- read.csv("pre_final.csv",header = F)
pre[1,1] <- "name"
pre2 <- t(pre)
colnames(pre2) <- pre2[1,]
pre2 <- pre2[-1,]
pre2 <- as.data.frame(pre2)
rownames(pre2) <- pre2$rid
pre2 <- pre2[,-1]
pre2 <- merge(E,pre2,by="row.names")

pre2$`Gene Symbol` <- gsub(" ","",pre2$`Gene Symbol`)

pre2 <- pre2[!duplicated(pre2$`Gene Symbol`),]
rownames(pre2) <- pre2$`Gene Symbol`
pre2 <- pre2[,-c(1:5)]


DES_all3 <-pre2 

#########
col <- c(rep("CBD",8),rep("DMSO",8))
a<- apply(DES_all3,2,function(x) as.numeric(as.character(x)))
rownames(a) <- rownames(DES_all3)
duplicated(c4)
library(limma)
group <- as.matrix(col)
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))

rownames(design) = colnames(a)
cc<- paste0("CBD", " - ", "DMSO")
contrast.matrix <- makeContrasts(cc, levels = design) 


fit <- lmFit(a,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG_ot <- topTable(fit2, adjust.method="fdr", coef=1,sort.by="logFC",n = Inf)
DEG_ot$ID <- rownames(DEG_ot)
library(ggplot2)
library(ggrepel)
DEG_ot <- as.data.frame(DEG_ot)
DEG_ot_label <- DEG_ot[DEG_ot$P.Value < 3e-4,]

DEG_ot$log2fc <- DEG_ot$logFC * log2(10)

ggplot(DEG_ot, aes(log2fc, -log10(P.Value))) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "#999999") +
  geom_point(aes(size = -log10(P.Value), color = -log10(P.Value))) +
  scale_color_gradientn(values = seq(0, 1, 0.1),
                        colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  scale_size_continuous(range = c(1, 3)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1, 0.7),
        legend.justification = c(0, 1)) +
  guides(col = guide_colourbar(title = "-Log10(p-value)"),
         size = "none") +
  geom_text_repel(data = DEG_ot_label , aes(label = DEG_ot_label [,7]), color = "black", size = 3) +
  xlab("Log2FC") +
  ylab("-Log10(p-value)")