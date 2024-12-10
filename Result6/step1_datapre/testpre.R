library(GenomicFeatures)
library(Homo.sapiens)
library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(clusterProfiler)
library(dplyr)
library(rtracklayer)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(gkmSVM)

genome <- BSgenome.Hsapiens.UCSC.hg19
all_genes <- suppressMessages(genes(Homo.sapiens))

all_gene_TSS <- resize(all_genes, 1)
promoter <- promoters(all_gene_TSS, upstream = 201, downstream = 0)
p <- as.data.frame(promoter)
p$GENEID <- rownames(p)
ids <- bitr(p$GENEID,'ENTREZID','SYMBOL','org.Hs.eg.db')
p2 <- data.frame(ENTREZID = p1$GENEID,SYMBOL=p1$GENEID)
ids <- rbind(ids,p2)
rownames(ids) <- ids$ENTREZID
pp <- cbind(p,ids)

sig <- p[,1:3]
sig$start <- sig$start
write.table(sig,"all4.bed",row.names = F,col.names = F,sep="\t",quote = F)
bed_file <- import("all4.bed", format = "BED")
sequences <- getSeq(genome, bed_file)

write_sequences_to_fasta <- function(sequences, bed_file, filepath) {
  fasta_lines <- character()
  for (i in seq_along(sequences)) {
    chr <- bed_file[i, "seqnames"]
    start <- bed_file[i, "start"]
    end <- bed_file[i, "end"]
    seq_name <- paste0(chr, ":", start, "-", end)
    fasta_lines <- c(fasta_lines, paste0(">", seq_name), as.character(sequences[[i]]))
  }
  writeLines(fasta_lines, con = filepath)
}


write_sequences_to_fasta(sequences, sig, filepath = "promoter_sequences1.fa")
promoter <- promoters(all_gene_TSS, upstream = 101, downstream = 100)
p <- as.data.frame(promoter)
p$GENEID <- rownames(p)
ids <- bitr(p$GENEID,'ENTREZID','SYMBOL','org.Hs.eg.db')
p2 <- data.frame(ENTREZID = p1$GENEID,SYMBOL=p1$GENEID)
ids <- rbind(ids,p2)
rownames(ids) <- ids$ENTREZID
pp <- cbind(p,ids)

sig <- p[,1:3]
sig$start <- sig$start
write.table(sig,"all4.bed",row.names = F,col.names = F,sep="\t",quote = F)
bed_file <- import("all4.bed", format = "BED")
sequences <- getSeq(genome, bed_file)

write_sequences_to_fasta <- function(sequences, bed_file, filepath) {
  fasta_lines <- character()
  for (i in seq_along(sequences)) {
    chr <- bed_file[i, "seqnames"]
    start <- bed_file[i, "start"]
    end <- bed_file[i, "end"]
    seq_name <- paste0(chr, ":", start, "-", end)
    fasta_lines <- c(fasta_lines, paste0(">", seq_name), as.character(sequences[[i]]))
  }
  writeLines(fasta_lines, con = filepath)
}

write_sequences_to_fasta(sequences, sig, filepath = "promoter_sequences2.fa")

promoter <- promoters(all_gene_TSS, upstream = 1, downstream = 200)
p <- as.data.frame(promoter)
p$GENEID <- rownames(p)
ids <- bitr(p$GENEID,'ENTREZID','SYMBOL','org.Hs.eg.db')

rownames(ids) <- ids$ENTREZID
pp <- cbind(p,ids)

sig <- p[,1:3]
sig$start <- sig$start
write.table(sig,"all4.bed",row.names = F,col.names = F,sep="\t",quote = F)
bed_file <- import("all4.bed", format = "BED")
sequences <- getSeq(genome, bed_file)

write_sequences_to_fasta <- function(sequences, bed_file, filepath) {
  fasta_lines <- character()
  for (i in seq_along(sequences)) {
    chr <- bed_file[i, "seqnames"]
    start <- bed_file[i, "start"]
    end <- bed_file[i, "end"]
    seq_name <- paste0(chr, ":", start, "-", end)
    fasta_lines <- c(fasta_lines, paste0(">", seq_name), as.character(sequences[[i]]))
  }
  writeLines(fasta_lines, con = filepath)
}

write_sequences_to_fasta(sequences, sig, filepath = "promoter_sequences3.fa")
##################
data1 <- as.data.frame(readDNAStringSet("promoter_sequences1.fa"))
data2 <- as.data.frame(readDNAStringSet("promoter_sequences2.fa"))
data3 <- as.data.frame(readDNAStringSet("promoter_sequences3.fa"))
rownames(data1) <- rownames(sig)
rownames(data2) <- rownames(sig)
rownames(data3) <- rownames(sig)
data1$gene <- rownames(data1) 
data2$gene <- rownames(data2) 
data3$gene <- rownames(data3) 
data_all <- rbind(data1,data2,data3)
##############
data_all2 <- data_all
data_all2$gene <- 0
colnames(data_all2) <- c("sequence","labels")
write.table(data_all2,"test.txt",col.names = F,row.names = F,sep="\t",quote = F)