library(data.table)
Rawcounts <- read.csv("FC_CBD.csv")
Rawcounts <- na.omit(Rawcounts)
library(ggplot2)
counts <- Rawcounts[,7:ncol(Rawcounts)] 
rownames(counts) <- Rawcounts$Geneid 
geneid_efflen <- subset(Rawcounts,select = c("Geneid","Length"))
library(dplyr)
geneid_efflen <- subset(Rawcounts,select = c("Geneid","Length"))
library(dplyr)

colnames(geneid_efflen) <- c("geneid","efflen")
geneid_efflen_fc <- geneid_efflen

dim(geneid_efflen)
efflen <- geneid_efflen[match(rownames(counts),
                              geneid_efflen$geneid),
                        "efflen"]



counts2TPM <- function(count=count, efflength=efflen){
  RPK <- count/(efflength/1000)      
  PMSC_rpk <- sum(RPK)/1e6        
  RPK/PMSC_rpk                    
}  
counts <- as.matrix(counts)
rownames(counts) <- Rawcounts$Geneid 
tpm <- as.data.frame(apply(counts,2,counts2TPM))
colSums(tpm)
colnames(tpm) <- colnames(counts)
rownames(tpm) <- rownames(counts)
write.csv(tpm, file = 'tpm.csv',row.names = T)



library(limma)
DES <- tpm
DES <- log10(DES+1)
col<- data.frame( x = colnames(DES),condition = c(rep("model",3),rep("cbd",3)) )
group <- col$condition
group <- as.matrix(group)
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(DES)
contrast.matrix <- makeContrasts("cbd-model", levels = design) 
fit <- lmFit(DES,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG_ot <- topTable(fit2, adjust.method="fdr", coef=1,sort.by="logFC",n = Inf)
DEG_ot$GEN <- rownames(DEG_ot) 
DEG_ot <- na.omit(DEG_ot)
deg <- merge(DES,DEG_ot,by="row.names")
write.csv(deg,"CBD_DEG.csv")
