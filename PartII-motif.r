BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("org.Mm.eg.db", character.only = TRUE)
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("ChIPseeker")
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(enrichplot)
library("org.Mm.eg.db", character.only = TRUE)
library("DESeq2")
library(ggplot2)
library(pheatmap)
library(dplyr)
library(EnhancedVolcano)

ATAC_trx <- read.delim("~/Desktop/ATAC-seq/Trx_combined_tss.txt", comment.char="#")
ATAC_trx$Location <- paste(ATAC_trx$Chr,ATAC_trx$Start,sep="_")
ATAC_trx_fc <- ATAC_trx[,5:10]
row.names(ATAC_trx_fc) <- ATAC_trx$Location

info_ATAC_trx <- data.frame(X = c("Trx_WT1","Trx_WT2","Trx_WT3","Trx_Mut1","Trx_Mut2","Trx_Mut3"), 
                             Condition = c("WT","WT","WT","Mut","Mut","Mut") )
dds_ATAC_trx <- DESeqDataSetFromMatrix(countData = ATAC_trx_fc, colData = info_ATAC_trx, design = ~Condition)
keep_ATAC_trx <- rowSums2(counts(dds_ATAC_trx)) >= 300
dds_ATAC_trx <- dds_ATAC_trx[keep_ATAC_trx,]
DE_ATAC_trx <- DESeq(dds_ATAC_trx)
plotMA(results(DE_ATAC_trx, contrast = c('Condition','Mut','WT'), alpha = 0.05), ylim =c(-5,5),main="ATAC_trx: Mut vs WT")

nCounts_ATAC_trx <- counts(DE_ATAC_trx,normalized = T)
write.csv(nCounts_ATAC_trx,"~/Desktop/ATAC-seq/ATACseq_trx_normCnt.txt")
Cnt_ATAC_trx <- read.csv("~/Desktop/ATAC-seq/ATACseq_trx_normCnt.txt", row.names=1)

#PCA
tnormCounts_ATAC_trx <- t(nCounts_ATAC_trx)
pca_ATAC_trx <- prcomp(tnormCounts_ATAC_trx, scale. = TRUE)
screeplot(pca_ATAC_trx, type = "l", main = "Screeplot for ATAC_germ")
abline(h=14000, col = 'red', lty =2)

summary(pca_ATAC_trx)
pca_ATAC_trx$x
pca_ATAC_trx <- data.frame(info_ATAC_trx[,2], pca_ATAC_trx$x[,1:2])
colnames(pca_ATAC_trx) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_ATAC_trx, aes(PC1, PC2, group=Condition)) + geom_point(size=6, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15,16))+
  scale_color_manual(values = c("red3", "blue4"))+
  labs(title="ATAC_trx" ,x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-220,250) + ylim(-250, 150)

res_ATAC_trx <- results(DE_ATAC_trx , contrast = c("Condition", "Mut", "WT"), alpha = 0.05, cooksCutoff = FALSE)
res_ATAC_trx <- res_ATAC_trx[order(res_ATAC_trx$padj),]
write.csv(res_ATAC_trx, "~/Desktop/ATAC-seq/ATACseq_trx_pval.txt")
Pval_ATAC_trx <- read.csv("~/Desktop/ATAC-seq/ATACseq_trx_pval.txt", row.names = 1)
row.names(Pval_ATAC_trx) <- Pval_ATAC_trx$Location

row.names(ATAC_trx) <- ATAC_trx$Location
Pval_ATAC_trx <- na.omit(Pval_ATAC_trx)
Pval_ATAC_trx <- merge(Pval_ATAC_trx, ATAC_trx ,by=0)
Pval_ATAC_trx <- Pval_ATAC_trx[,-1]
Pval_ATAC_trx <- Pval_ATAC_trx[order(Pval_ATAC_trx$padj),]
write.csv(Pval_ATAC_trx, "~/Desktop/ATAC-seq/ATACseq_trx_pval.txt")
Pval_ATAC_trx <- read.csv("~/Desktop/ATAC-seq/ATACseq_trx_pval.txt", row.names = 1)

EnhancedVolcano(Pval_ATAC_trx, lab = NA, x='log2FoldChange',y='padj',
                title = 'ATAC_trx: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','blue'), colAlpha = 0.5,
                FCcutoff = 1, legendPosition = 'right') + xlim(-3,3)+ ylim(0,6)

Genelist_ATAC_trx <- Pval_ATAC_trx$log2FoldChange
names(Genelist_ATAC_trx) <- Pval_ATAC_trx$Entrez.ID
Genelist_ATAC_trx <- na.omit(Genelist_ATAC_trx)
Genelist_ATAC_trx <- sort(Genelist_ATAC_trx, decreasing = TRUE)
Gse_ATAC_trx <- gseGO(geneList = Genelist_ATAC_trx, ont = "BP", keyType = "ENSEMBL", 
                       minGSSize = 50, maxGSSize = 400, pvalueCutoff = 0.05, verbose = TRUE,
                       OrgDb = "org.Mm.eg.db",pAdjustMethod = "BH",eps=0)
Gsea_ATAC_trx <- as.data.frame(Gse_ATAC_trx)
require(DOSE)
dotplot(Gse_ATAC_trx, font.size=8, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("GSEA: ATAC_germ")

Overlap <- read_csv("~/Dropbox/Mac/Desktop/ATAC-seq/Overlap.csv")
Overlap_DDR <- read_csv("~/Dropbox/Mac/Desktop/ATAC-seq/Overlap_DDR.csv")
Overlap_HSC <- read_csv("~/Dropbox/Mac/Desktop/ATAC-seq/Overlap_HSC.csv")
Overlap_Remodel <- read_csv("~/Dropbox/Mac/Desktop/ATAC-seq/Overlap_Remodel.csv")
Overlap$Location <- paste(Overlap$Chr,Overlap$Start,sep="_")
row.names(Overlap) <- Overlap$Location
Overlap_DDR$Location <- paste(Overlap_DDR$Chr,Overlap_DDR$Start,sep="_")
row.names(Overlap_DDR) <- Overlap_DDR$Location
Overlap_HSC$Location <- paste(Overlap_HSC$Chr,Overlap_HSC$Start,sep="_")
row.names(Overlap_HSC) <- Overlap_HSC$Location
Overlap_Remodel$Location <- paste(Overlap_Remodel$Chr,Overlap_Remodel$Start,sep="_")
row.names(Overlap_Remodel) <- Overlap_Remodel$Location

ATAC_trx_norm <- merge(Cnt_ATAC_trx, ATAC_trx ,by=0)
ATAC_trx_norm <- ATAC_trx_norm[!duplicated(ATAC_trx_norm$Gene.Name),]
row.names(ATAC_trx_norm)<-ATAC_trx_norm$Entrez.ID
ATAC_germ_norm <- merge(Cnt_ATAC_germ, ATAC_germline ,by=0)

row.names(Overlap)<- Overlap$Entrez.ID
Overlap_DDR <- Overlap_DDR[,-1]
row.names(Overlap_DDR)<- Overlap_DDR$converted_alias
row.names(Overlap_HSC)<- Overlap_HSC$Entrez.ID
row.names(Overlap_Remodel)<- Overlap_Remodel$Entrez.ID

ATAC_Cnt_overlap <- merge(ATAC_trx_norm,ATAC_germ_norm, by=0)
ATAC_Cnt_overlap <- merge(ATAC_Cnt_overlap, Overlap, by=0)
ATAC_Cnt_overlap <- ATAC_Cnt_overlap[,-1:-3]
ATAC_Cnt_overlap <- ATAC_Cnt_overlap[,-7:-21]
row.names(ATAC_Cnt_overlap)<-ATAC_Cnt_overlap$Entrez.ID
ATAC_Cnt_overlap <- ATAC_Cnt_overlap[,1:10]
ATAC_Cnt_overlap$Entrez.ID <- row.names(ATAC_Cnt_overlap)

ATAC_Cnt_overlap_DDR <- merge(ATAC_Cnt_overlap, Overlap_DDR,by=c('Entrez.ID','Entrez.ID'))
ATAC_Cnt_overlap_DDR <- ATAC_Cnt_overlap_DDR[,-12:-18]
ATAC_Cnt_overlap_HSC <- merge(ATAC_Cnt_overlap, Overlap_HSC,by=c('Entrez.ID','Entrez.ID'))
ATAC_Cnt_overlap_HSC <- ATAC_Cnt_overlap_HSC[,-12:-19]
ATAC_Cnt_overlap_Remodel <- merge(ATAC_Cnt_overlap, Overlap_Remodel,by=c('Entrez.ID','Entrez.ID'))
ATAC_Cnt_overlap_Remodel <- ATAC_Cnt_overlap_Remodel[,-12:-19]
write.csv(ATAC_Cnt_overlap, "~/Desktop/ATAC-seq/ATAC_Cnt_overlap_all.csv")
write.csv(ATAC_Cnt_overlap_DDR, "~/Desktop/ATAC-seq/ATAC_Cnt_overlap_DDR.csv")
write.csv(ATAC_Cnt_overlap_HSC, "~/Desktop/ATAC-seq/ATAC_Cnt_overlap_HSC.csv")
write.csv(ATAC_Cnt_overlap_Remodel, "~/Desktop/ATAC-seq/ATAC_Cnt_overlap_Remodel.csv")

TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
ATAC_trx_wt_peak <- readPeakFile("~/Desktop/ATAC-seq/Peak/Trx_WT_summits.bed")
ATAC_trx_mut_peak <- readPeakFile("~/Desktop/ATAC-seq/Peak/Trx_Mut_summits.bed")

promoter <- getPromoters(TxDb=TxDb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(ATAC_germ_wt_peak, windows=promoter)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")
