library(DESeq2)
library(dplyr)
library(stringr)
library(gdata)
library(ggplot2)
library(venneuler)
library(WGCNA)
library(tidyverse)
library(pheatmap)
#browseVignettes("DESeq2")

setwd("C:/Users/P309374/Documents/PhD/Macquarie/LTF TEtranscripts/")

#full-factorial analysis for main effects (results described in main text and Table S2)
countData <- read.csv("Counts_final.csv",h=T,row.names=1)
colnames(countData)
#updatedCountData <- read.csv("merged_p1_Counts.csv",h=T,row.names=1)
#countData[colnames(updatedCountData[2:23])]<- updatedCountData[2:23]
#countData_TEs <- countData[str_detect(rownames(countData), "ID="), ]
colData <- read.csv("all_samples_metadata.csv",header=TRUE,row.names=1)
colData <- subset(colData, colData$tech_rep=="main")
colData$index<- match(rownames(colData),colnames(countData))
colData<-colData %>% arrange(index)
colData$taxon <- as.factor(colData$taxon)
colData$sex <- as.factor(colData$sex)
colData$tissue <- as.factor(colData$tissue)

#dataset partinioned into tissue type 
library(edgeR)
ovaries_meta <- subset(colData, colData$sex == "Female" & colData$tissue=="gonad")
ovaries_count <- countData[which(colnames(countData)%in% rownames(ovaries_meta))]

testes_meta <- subset(colData, colData$sex == "Male" & colData$tissue=="gonad")
testes_count <- countData[which(colnames(countData)%in% rownames(testes_meta))]
ovaries_meta$x <- ifelse(ovaries_meta$taxon == "acu", '1',
                         ifelse(ovaries_meta$taxon == "hec", '2',
                                ifelse(ovaries_meta$taxon == "cin", 3,
                                       ifelse(ovaries_meta$taxon == "acu-hec",4, 5))))

testes_meta$x <- ifelse(testes_meta$taxon == "acu", '1',
                        ifelse(testes_meta$taxon == "hec", '2',
                               ifelse(testes_meta$taxon == "cin", 3,
                                      ifelse(testes_meta$taxon == "acu-hec",4, 5))))
liver_meta <- subset(colData,colData$tissue=="liver")
liver_count <- countData[which(colnames(countData)%in% rownames(liver_meta))]
liver_meta$x <- ifelse(liver_meta$taxon == "acu" & liver_meta$sex=="Female", '1',
                       ifelse(liver_meta$taxon == "hec" & liver_meta$sex=="Female", '2',
                              ifelse(liver_meta$taxon == "cin" & liver_meta$sex=="Female", 3,
                                     ifelse(liver_meta$taxon == "acu-hec" & liver_meta$sex=="Female", 4, 
                                            ifelse(liver_meta$taxon == "hec-cin" & liver_meta$sex=="Male", 5,
                                                   ifelse(liver_meta$taxon == "acu" & liver_meta$sex=="Male", 6,
                                                          ifelse(liver_meta$taxon == "hec" & liver_meta$sex=="Male", 7,
                                                                 ifelse(liver_meta$taxon == "cin" & liver_meta$sex=="Male",8,
                                                                        ifelse(liver_meta$taxon == "acu-hec" & liver_meta$sex=="Male",9,10)))))))))

ovaries_meta$group <- paste0(ovaries_meta$x,"_", ovaries_meta$mtdna)
liver_meta$group <- ifelse(liver_meta$sex=="Female", paste0(liver_meta$x,"_",liver_meta$mtdna),liver_meta$x)
ovaries_meta$group <- as.factor(ovaries_meta$group)
liver_meta$group<- as.factor(liver_meta$group)

dds_ovaries <- DESeqDataSetFromMatrix(countData = ovaries_count,colData = ovaries_meta, design = ~ group)
counts(dds_ovaries)
keep_o <- rowSums(counts(dds_ovaries)) >= 10
dds_ovaries <- dds_ovaries[keep_o,]
dds_ovaries<- estimateSizeFactors(dds_ovaries)
dds_ovaries <- DESeq(dds_ovaries, modelMatrixType="standard")
resultsNames(dds_ovaries)


dds_testes <- DESeqDataSetFromMatrix(countData = testes_count,colData = testes_meta, design = ~ x)
counts(dds_testes)
keep_t <- rowSums(counts(dds_testes)) >= 10
dds_testes <- dds_testes[keep_t,]
dds_testes<- estimateSizeFactors(dds_testes)
dds_testes <- DESeq(dds_testes, modelMatrixType="standard")
resultsNames(dds_testes)


dds_liver <- DESeqDataSetFromMatrix(countData = liver_count,colData = liver_meta, design = ~ group)
counts(dds_liver)
keep_l <- rowSums(counts(dds_liver)) >= 10
dds_liver <- dds_liver[keep_l,]
dds_liver<- estimateSizeFactors(dds_liver)
dds_liver <- DESeq(dds_liver, modelMatrixType="standard")
resultsNames(dds_liver)
dim(liver_count)
#liver_count_filtered <- liver_count[1:69]
#liver_meta_filtered <- liver_meta[1:69,]
#liver_meta_filtered$group <- factor(paste0(liver_meta_filtered$taxon, liver_meta_filtered$sex))

#RUV
library(RUVSeq)
library(edgeR)
# Ovaries RUVseq
dge_o<- DGEList(counts = counts(dds_ovaries), samples = row.names(ovaries_meta), group = ovaries_meta$group)
dge_o <- calcNormFactors(dge_o, method = "TMM")
counts_RUV_o <- cpm(dge_o, log=F)
counts_RUV_o <- round(counts_RUV_o)

set<- newSeqExpressionSet(counts = counts_RUV_o, phenoData = data.frame(ovaries_meta$group, row.names= colnames(counts(dds_ovaries))))
keep_set <- rowSums(counts(set)) >= 10
set <- set[keep_set,]

library(RColorBrewer)
colors <- brewer.pal(11, "Spectral")
x<- ovaries_meta$group
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[set$ovaries_meta.group])
plotPCA(set, col=colors[x], cex=1.2)

#Get empirical control genes
design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]
res <- residuals(fit, type="deviance")

#RUV
set1 <- RUVr(set, empirical, k=1, res)
plotRLE(set1, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=1")
plotPCA(set1, col=colors[x], main="PCA - k=1")

set2 <- RUVr(set, empirical, k=2, res)
plotRLE(set2, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=2")
plotPCA(set2, col=colors[x], main="PCA - k=2")

set3 <- RUVr(set, empirical, k=3, res)
plotRLE(set2, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=3")
plotPCA(set2, col=colors[x], main="PCA - k=3")

set4 <- RUVr(set, empirical, k=4, res)
plotRLE(set4, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=4")
plotPCA(set4, col=colors[x], main="PCA - k=4")

set5 <- RUVr(set, empirical, k=5, res)
plotRLE(set5, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=5")
plotPCA(set5, col=colors[x], main="PCA - k=5")

set6 <- RUVr(set, empirical, k=6, res)
plotRLE(set6, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=6")
plotPCA(set6, col=colors[x], main="PCA - k=6")

set7 <- RUVr(set, empirical, k=7, res)
plotRLE(set7, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=7")
plotPCA(set7, col=colors[x], main="PCA - k=7")

set8 <- RUVr(set, empirical, k=8, res)
plotRLE(set8, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=8")
plotPCA(set8, col=colors[x], main="PCA - k=8")

set9 <- RUVr(set, empirical, k=9, res)
plotRLE(set9, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=9")
plotPCA(set9, col=colors[x], main="PCA - k=9")

set10 <- RUVr(set, empirical, k=10, res)
plotRLE(set10, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=10")
plotPCA(set10, col=colors[x], main="PCA - k=10")

set11 <- RUVr(set, empirical, k=11, res)
plotRLE(set11, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=11")
plotPCA(set11, col=colors[x], main="PCA - k=11")

#design <- model.matrix(~x+W_1+W_2+W_3+W_4+W_5+W_6+W_7+W_8+W_9+W_10+W_11, data = pData(set11))
#y <- DGEList(counts=counts(set11), group=x)
#y <- calcNormFactors(y, method="upperquartile")
#y <- estimateGLMCommonDisp(y, design)
#y <- estimateGLMTagwiseDisp(y, design)

#fit <- glmFit(y, design)
#lrt <- glmLRT(fit, coef=3)
#topTags_o<- topTags(lrt, n=100000L, p.value =0.05)
#dim(topTags_o$table)
#topTags_o_TEs<- topTags_o$table[str_detect(row.names(topTags_o$table), "ID="), ]
#length(which(topTags_o_TEs$logFC > 0))
### 3938 upregulated TEs
#length(which(topTags_o_TEs$logFC < 0))
### 269 downregulated TEs

### Load w chromosome insertions ###
W_chrTEs <- read.table("chrW_TEs.txt", fill=T)
sum(is.na(W_chrTEs$V13))
# No NAs in W chr transcript IDs ##

## Try DESeq2
dds_o_RUV <- DESeqDataSetFromMatrix(countData = counts(set11),
                                    colData = pData(set11),
                                    design = ~ 0+ovaries_meta.group + W_1+W_2+W_3+W_4+W_5+W_6+W_7+W_8+W_9+W_10+W_11)
dds_o_RUV <- DESeq(dds_o_RUV)
resultsNames(dds_o_RUV)
contrast_hec.cin_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","1_acu","2_hec"), alpha = 0.01)
table(contrast_hec.cin_RUV_o$padj < 0.01)
contrast_hec.cin_RUV_o <- subset(contrast_hec.cin_RUV_o, padj < 0.01)
contrast_hec.cin_RUV_o_TEs<- contrast_hec.cin_RUV_o[str_detect(contrast_hec.cin_RUV_o@rownames, "RND"), ]
#acu_hec_WTEs_o<- subset(contrast_hec.cin_RUV_o_TEs, unlist(strsplit(contrast_hec.cin_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin_RUV_o_TEs@rownames),3)],
                            contrast_hec.cin_RUV_o_TEs$log2FoldChange)      
class_acu_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin_RUV_o_TEs@rownames),3)],
                            contrast_hec.cin_RUV_o_TEs$log2FoldChange)
acu_hec_up_o<- length(which(contrast_hec.cin_RUV_o_TEs$log2FoldChange > 0))
### 13 upregulated TEs
acu_hec_dw_o<- (length(which(contrast_hec.cin_RUV_o_TEs$log2FoldChange < 0)))*-1
### 20 downregulated TEs
acu_hec_up_o_W<- length(which(acu_hec_WTEs_o$log2FoldChange > 0))
acu_hec_dw_o_W<- (length(which(acu_hec_WTEs_o$log2FoldChange < 0)))*-1

contrast_hec.cin2_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","2_hec","3_cin"), alpha = 0.01)
table(contrast_hec.cin2_RUV_o$padj < 0.01)
contrast_hec.cin2_RUV_o <- subset(contrast_hec.cin2_RUV_o, padj < 0.01)
contrast_hec.cin2_RUV_o_TEs<- contrast_hec.cin2_RUV_o[str_detect(contrast_hec.cin2_RUV_o@rownames, "RND"), ]
#hec_cin_WTEs_o<- subset(contrast_hec.cin2_RUV_o_TEs, unlist(strsplit(contrast_hec.cin2_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin2_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin2_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec_cin_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin2_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin2_RUV_o_TEs@rownames)])[seq(2,2*length(contrast_hec.cin2_RUV_o_TEs@rownames),3)],
                            contrast_hec.cin2_RUV_o_TEs$log2FoldChange)      
class_hec_cin_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin2_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin2_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin2_RUV_o_TEs@rownames),3)],
                            contrast_hec.cin2_RUV_o_TEs$log2FoldChange)
hec_cin_up_o<- length(which(contrast_hec.cin2_RUV_o_TEs$log2FoldChange > 0))
### 4 upregulated TEs
hec_cin_dw_o<- length(which(contrast_hec.cin2_RUV_o_TEs$log2FoldChange < 0))*-1
### 11 downregulated TEs
hec_cin_up_o_W<- length(which(hec_cin_WTEs_o$log2FoldChange > 0))
hec_cin_dw_o_W<- (length(which(hec_cin_WTEs_o$log2FoldChange < 0)))*-1



contrast_hec.cin3_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","5_hec","2_hec"), alpha = 0.01)
table(contrast_hec.cin3_RUV_o$padj < 0.01)
contrast_hec.cin3_RUV_o <- subset(contrast_hec.cin3_RUV_o, padj < 0.01)
contrast_hec.cin3_RUV_o_TEs<- contrast_hec.cin3_RUV_o[str_detect(contrast_hec.cin3_RUV_o@rownames, "RND"), ]
#hec.cin_hec_WTEs_o<- subset(contrast_hec.cin3_RUV_o_TEs, unlist(strsplit(contrast_hec.cin3_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin3_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin3_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec.cin_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin3_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin3_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin3_RUV_o_TEs@rownames),3)],
                                contrast_hec.cin3_RUV_o_TEs$log2FoldChange)      
class_hec.cin_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin3_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin3_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin3_RUV_o_TEs@rownames),3)],
                                contrast_hec.cin3_RUV_o_TEs$log2FoldChange)
hec.cin_hec_up_o<- length(which(contrast_hec.cin3_RUV_o_TEs$log2FoldChange > 0))
### 151 upregulated TEs
hec.cin_hec_dw_o<- -1*length(which(contrast_hec.cin3_RUV_o_TEs$log2FoldChange < 0))
### 24 downregulated TEs
hec.cin_hec_up_o_W<- length(which(hec.cin_hec_WTEs_o$log2FoldChange > 0))
hec.cin_hec_dw_o_W<- (length(which(hec.cin_hec_WTEs_o$log2FoldChange < 0)))*-1


contrast_hec.cin4_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","5_hec","3_cin"), alpha = 0.01)
table(contrast_hec.cin4_RUV_o$padj < 0.01)
contrast_hec.cin4_RUV_o <- subset(contrast_hec.cin4_RUV_o, padj < 0.01)
contrast_hec.cin4_RUV_o_TEs<- contrast_hec.cin4_RUV_o[str_detect(contrast_hec.cin4_RUV_o@rownames, "RND"), ]
#hec.cin_cin_WTEs_o<- subset(contrast_hec.cin4_RUV_o_TEs, unlist(strsplit(contrast_hec.cin4_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin4_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin4_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec.cin_cin_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin4_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin4_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin4_RUV_o_TEs@rownames),3)],
                                contrast_hec.cin4_RUV_o_TEs$log2FoldChange)      
class_hec.cin_cin_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin4_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin4_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin4_RUV_o_TEs@rownames),3)],
                                contrast_hec.cin4_RUV_o_TEs$log2FoldChange)
hec.cin_cin_up_o<- length(which(contrast_hec.cin4_RUV_o_TEs$log2FoldChange > 0))
### 167 upregulated TEs
hec.cin_cin_dw_o<- -1*length(which(contrast_hec.cin4_RUV_o_TEs$log2FoldChange < 0))
### 18 downregulated TEs
hec.cin_cin_up_o_W<- length(which(hec.cin_cin_WTEs_o$log2FoldChange > 0))
hec.cin_cin_dw_o_W<- (length(which(hec.cin_cin_WTEs_o$log2FoldChange < 0)))*-1


contrast_hec.cin5_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","4_hec","2_hec"), alpha = 0.01)
table(contrast_hec.cin5_RUV_o$padj < 0.01)
contrast_hec.cin5_RUV_o <- subset(contrast_hec.cin5_RUV_o, padj < 0.01)
contrast_hec.cin5_RUV_o_TEs<- contrast_hec.cin5_RUV_o[str_detect(contrast_hec.cin5_RUV_o@rownames, "RND"), ]
#acu.hec_hec_WTEs_o<- subset(contrast_hec.cin5_RUV_o_TEs, unlist(strsplit(contrast_hec.cin5_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin5_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin5_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu.hec_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin5_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin5_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin5_RUV_o_TEs@rownames),3)],
                                contrast_hec.cin5_RUV_o_TEs$log2FoldChange)      
class_acu.hec_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin5_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin5_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin5_RUV_o_TEs@rownames),3)],
                                contrast_hec.cin5_RUV_o_TEs$log2FoldChange)
acu.hec_hec_up_o<-length(which(contrast_hec.cin5_RUV_o_TEs$log2FoldChange > 0))
### 130 upregulated TEs
acu.hec_hec_dw_o<- -1*length(which(contrast_hec.cin5_RUV_o_TEs$log2FoldChange < 0))
### 52 downregulated TEs
acu.hec_hec_up_o_W<- length(which(acu.hec_hec_WTEs_o$log2FoldChange > 0))
acu.hec_hec_dw_o_W<- (length(which(acu.hec_hec_WTEs_o$log2FoldChange < 0)))*-1

contrast_hec.cin6_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","4_hec","1_acu"), alpha = 0.01)
table(contrast_hec.cin6_RUV_o$padj < 0.01)
contrast_hec.cin6_RUV_o <- subset(contrast_hec.cin6_RUV_o, padj < 0.01)
contrast_hec.cin6_RUV_o_TEs<- contrast_hec.cin6_RUV_o[str_detect(contrast_hec.cin6_RUV_o@rownames, "RND"), ]
#acu.hec_acu_WTEs_o<- subset(contrast_hec.cin6_RUV_o_TEs, unlist(strsplit(contrast_hec.cin6_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin6_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin6_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu.hec_acu_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin6_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin6_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin6_RUV_o_TEs@rownames),3)],
                                contrast_hec.cin6_RUV_o_TEs$log2FoldChange)      
class_acu.hec_acu_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin6_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin6_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin6_RUV_o_TEs@rownames),3)],
                                contrast_hec.cin6_RUV_o_TEs$log2FoldChange)
acu.hec_acu_up_o<- length(which(contrast_hec.cin6_RUV_o_TEs$log2FoldChange > 0))
### 100 upregulated TEs
acu.hec_acu_dw_o<- -1*length(which(contrast_hec.cin6_RUV_o_TEs$log2FoldChange < 0))
### 32 downregulated TEs
acu.hec_acu_up_o_W<- length(which(acu.hec_acu_WTEs_o$log2FoldChange > 0))
acu.hec_acu_dw_o_W<- (length(which(acu.hec_acu_WTEs_o$log2FoldChange < 0)))*-1


contrast_hec.cin7_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","1_acu","3_cin"), alpha = 0.01)
table(contrast_hec.cin7_RUV_o$padj < 0.01)
contrast_hec.cin7_RUV_o <- subset(contrast_hec.cin7_RUV_o, padj < 0.01)
contrast_hec.cin7_RUV_o_TEs<- contrast_hec.cin7_RUV_o[str_detect(contrast_hec.cin7_RUV_o@rownames, "RND"), ]
#acu_cin_WTEs_o<- subset(contrast_hec.cin7_RUV_o_TEs, unlist(strsplit(contrast_hec.cin7_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin7_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin7_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu_cin_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin7_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin7_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin7_RUV_o_TEs@rownames),3)],
                            contrast_hec.cin7_RUV_o_TEs$log2FoldChange)      
class_acu_cin_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin7_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin7_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin7_RUV_o_TEs@rownames),3)],
                            contrast_hec.cin7_RUV_o_TEs$log2FoldChange)
acu_cin_up_o<- length(which(contrast_hec.cin7_RUV_o_TEs$log2FoldChange > 0))
### 24 upregulated TEs
acu_cin_dw_o<--1*length(which(contrast_hec.cin7_RUV_o_TEs$log2FoldChange < 0))
### 21 downregulated TEs
acu_cin_up_o_W<- length(which(acu_cin_WTEs_o$log2FoldChange > 0))
acu_cin_dw_o_W<- (length(which(acu_cin_WTEs_o$log2FoldChange < 0)))*-1


contrast_hec.cin8_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","5_cin","2_hec"), alpha = 0.01)
table(contrast_hec.cin8_RUV_o$padj < 0.01)
contrast_hec.cin8_RUV_o <- subset(contrast_hec.cin8_RUV_o, padj < 0.01)
contrast_hec.cin8_RUV_o_TEs<- contrast_hec.cin8_RUV_o[str_detect(contrast_hec.cin8_RUV_o@rownames, "RND"), ]
#hec.cinWcin_hec_WTEs_o<- subset(contrast_hec.cin8_RUV_o_TEs, unlist(strsplit(contrast_hec.cin8_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin8_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin8_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec.cinWcin_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin8_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin8_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin8_RUV_o_TEs@rownames),3)],
                                    contrast_hec.cin8_RUV_o_TEs$log2FoldChange)      
class_hec.cinWcin_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin8_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin8_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin8_RUV_o_TEs@rownames),3)],
                                    contrast_hec.cin8_RUV_o_TEs$log2FoldChange)
hec.cinWcin_hec_up_o<- length(which(contrast_hec.cin8_RUV_o_TEs$log2FoldChange > 0))
### 146 upregulated TEs
hec.cinWcin_hec_dw_o<- -1*length(which(contrast_hec.cin8_RUV_o_TEs$log2FoldChange < 0))
### 25 downregulated TEs
hec.cinWcin_hec_up_o_W<- length(which(hec.cinWcin_hec_WTEs_o$log2FoldChange > 0))
hec.cinWcin_hec_dw_o_W<- (length(which(hec.cinWcin_hec_WTEs_o$log2FoldChange < 0)))*-1


contrast_hec.cin9_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","5_cin","3_cin"), alpha = 0.01)
table(contrast_hec.cin9_RUV_o$padj < 0.01)
contrast_hec.cin9_RUV_o <- subset(contrast_hec.cin9_RUV_o, padj < 0.01)
contrast_hec.cin9_RUV_o_TEs<- contrast_hec.cin9_RUV_o[str_detect(contrast_hec.cin9_RUV_o@rownames, "RND"), ]
#hec.cinWcin_cin_WTEs_o<- subset(contrast_hec.cin9_RUV_o_TEs, unlist(strsplit(contrast_hec.cin9_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin9_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin9_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec.cinWcin_cin_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin9_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin9_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin9_RUV_o_TEs@rownames),3)],
                                    contrast_hec.cin9_RUV_o_TEs$log2FoldChange)      
class_hec.cinWcin_cin_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin9_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin9_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin9_RUV_o_TEs@rownames),3)],
                                    contrast_hec.cin9_RUV_o_TEs$log2FoldChange)
hec.cinWcin_cin_up_o<- length(which(contrast_hec.cin9_RUV_o_TEs$log2FoldChange > 0))
### 160 upregulated TEs
hec.cinWcin_cin_dw_o<- -1*length(which(contrast_hec.cin9_RUV_o_TEs$log2FoldChange < 0))
### 20 downregulated TEs
hec.cinWcin_cin_up_o_W<- length(which(hec.cinWcin_cin_WTEs_o$log2FoldChange > 0))
hec.cinWcin_cin_dw_o_W<- (length(which(hec.cinWcin_cin_WTEs_o$log2FoldChange < 0)))*-1


contrast_hec.cin10_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","4_acu","2_hec"), alpha = 0.01)
table(contrast_hec.cin10_RUV_o$padj < 0.01)
contrast_hec.cin10_RUV_o <- subset(contrast_hec.cin10_RUV_o, padj < 0.01)
contrast_hec.cin10_RUV_o_TEs<- contrast_hec.cin10_RUV_o[str_detect(contrast_hec.cin10_RUV_o@rownames, "RND"), ]
#acu.hecWacu_hec_WTEs_o<- subset(contrast_hec.cin10_RUV_o_TEs, unlist(strsplit(contrast_hec.cin10_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin10_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin10_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu.hecWacu_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin10_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin10_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin10_RUV_o_TEs@rownames),3)],
                                    contrast_hec.cin10_RUV_o_TEs$log2FoldChange)      
class_acu.hecWacu_hec_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin10_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin10_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin10_RUV_o_TEs@rownames),3)],
                                    contrast_hec.cin10_RUV_o_TEs$log2FoldChange)
acu.hecWacu_hec_up_o<-length(which(contrast_hec.cin10_RUV_o_TEs$log2FoldChange > 0))
### 132 upregulated TEs
acu.hecWacu_hec_dw_o<- -1*length(which(contrast_hec.cin10_RUV_o_TEs$log2FoldChange < 0))
### 38 downregulated TEs
acu.hecWacu_hec_up_o_W<- length(which(acu.hecWacu_hec_WTEs_o$log2FoldChange > 0))
acu.hecWacu_hec_dw_o_W<- (length(which(acu.hecWacu_hec_WTEs_o$log2FoldChange < 0)))*-1


contrast_hec.cin11_RUV_o<-results(dds_o_RUV, contrast = c("ovaries_meta.group","4_acu","1_acu"), alpha = 0.01)
table(contrast_hec.cin11_RUV_o$padj < 0.01)
contrast_hec.cin11_RUV_o <- subset(contrast_hec.cin11_RUV_o, padj < 0.01)
contrast_hec.cin11_RUV_o_TEs<- contrast_hec.cin11_RUV_o[str_detect(contrast_hec.cin11_RUV_o@rownames, "RND"), ]
#acu.hecWacu_acu_WTEs_o<- subset(contrast_hec.cin11_RUV_o_TEs, unlist(strsplit(contrast_hec.cin11_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin11_RUV_o_TEs@rownames)])[seq(1,4*length(contrast_hec.cin11_RUV_o_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu.hecWacu_acu_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin11_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin11_RUV_o_TEs@rownames)])[seq(2,3*length(contrast_hec.cin11_RUV_o_TEs@rownames),3)],
                                    contrast_hec.cin11_RUV_o_TEs$log2FoldChange)      
class_acu.hecWacu_acu_TEs_o<- cbind(unlist(strsplit(contrast_hec.cin11_RUV_o_TEs@rownames, ":")[1:length(contrast_hec.cin11_RUV_o_TEs@rownames)])[seq(3,3*length(contrast_hec.cin11_RUV_o_TEs@rownames),3)],
                                    contrast_hec.cin11_RUV_o_TEs$log2FoldChange)
acu.hecWacu_acu_up_o<- length(which(contrast_hec.cin11_RUV_o_TEs$log2FoldChange > 0))
### 92 upregulated TEs
acu.hecWacu_acu_dw_o<- -1*length(which(contrast_hec.cin11_RUV_o_TEs$log2FoldChange < 0))
### 20 downregulated TEs
acu.hecWacu_acu_up_o_W<- length(which(acu.hecWacu_acu_WTEs_o$log2FoldChange > 0))
acu.hecWacu_acu_dw_o_W<- (length(which(acu.hecWacu_acu_WTEs_o$log2FoldChange < 0)))*-1


## Plot DETEs for ovaries###
library(data.table)
fams_o <- list(sup_f_acu_hec_TEs_o,sup_f_hec_cin_TEs_o, sup_f_acu_cin_TEs_o, sup_f_acu.hec_acu_TEs_o,
               sup_f_acu.hec_hec_TEs_o, sup_f_acu.hecWacu_acu_TEs_o,sup_f_acu.hecWacu_hec_TEs_o,
               sup_f_hec.cin_cin_TEs_o,sup_f_hec.cin_hec_TEs_o,sup_f_hec.cinWcin_cin_TEs_o,
               sup_f_hec.cinWcin_hec_TEs_o)
fams_o <- lapply(fams_o, as.data.frame)
fams_o_up<- lapply(fams_o, subset, V2 > 0)
fams_o_dw<- lapply(fams_o, subset, V2 < 0)
fams_o_up <- bind_rows(fams_o_up, .id = "df_id")
fams_o_dw <- bind_rows(fams_o_dw, .id = "df_id")

class_o <- list(class_acu_hec_TEs_o,class_hec_cin_TEs_o, class_acu_cin_TEs_o, class_acu.hec_acu_TEs_o,
                class_acu.hec_hec_TEs_o, class_acu.hecWacu_acu_TEs_o,class_acu.hecWacu_hec_TEs_o,
                class_hec.cin_cin_TEs_o,class_hec.cin_hec_TEs_o,class_hec.cinWcin_cin_TEs_o,
                class_hec.cinWcin_hec_TEs_o)
class_o <- lapply(class_o, as.data.frame)
class_o_up<- lapply(class_o, subset, V2 > 0)
class_o_dw<- lapply(class_o, subset, V2 < 0)
class_o_up <- bind_rows(class_o_up, .id = "df_id")
class_o_dw <- bind_rows(class_o_dw, .id = "df_id")

up_o <-rbind(acu_hec_up_o, hec_cin_up_o, acu_cin_up_o, acu.hec_acu_up_o, acu.hec_hec_up_o,
             acu.hecWacu_acu_up_o,acu.hecWacu_hec_up_o,hec.cin_cin_up_o, hec.cin_hec_up_o,
             hec.cinWcin_cin_up_o,hec.cinWcin_hec_up_o)
dw_o <-rbind(acu_hec_dw_o, hec_cin_dw_o, acu_cin_dw_o, acu.hec_acu_dw_o, acu.hec_hec_dw_o,
             acu.hecWacu_acu_dw_o,acu.hecWacu_hec_dw_o,hec.cin_cin_dw_o, hec.cin_hec_dw_o,
             hec.cinWcin_cin_dw_o,hec.cinWcin_hec_dw_o)
contrasts<- c("acu_hec","hec_cin","acu_cin","acu.hec_acu(hecW)","acu.hec_hec(hecW)",
              "acu.hec_acu(acuW)","acu.hec_hec(acuW)","hec.cin_cin(hecW)",
              "hec.cin_hec(hecW)","hec.cin_cin(cinW)","hec.cin_hec(cinW)")
level_order_o<- contrasts
for(i in 1:length(fams_o_up$df_id)){fams_o_up$contrast[i]<-contrasts[as.numeric(fams_o_up$df_id[i])]}
for(i in 1:length(fams_o_dw$df_id)){fams_o_dw$contrast[i]<-contrasts[as.numeric(fams_o_dw$df_id[i])]}
for(i in 1:length(class_o_up$df_id)){class_o_up$contrast[i]<-contrasts[as.numeric(class_o_up$df_id[i])]}
for(i in 1:length(class_o_dw$df_id)){class_o_dw$contrast[i]<-contrasts[as.numeric(class_o_dw$df_id[i])]}
fams_o_up_counts <- fams_o_up %>% group_by(contrast, V1) %>% summarise(counts=n())
fams_o_dw_counts <- fams_o_dw %>% group_by(contrast, V1) %>% summarise(counts=n()*-1)
fams_o_counts <- rbind(fams_o_up_counts,fams_o_dw_counts)

class_o_up_counts <- class_o_up %>% group_by(contrast, V1) %>% summarise(counts=n())
class_o_dw_counts <- class_o_dw %>% group_by(contrast, V1) %>% summarise(counts=n()*-1)
class_o_counts <- rbind(class_o_up_counts,class_o_dw_counts)

DETEs_o <- as.data.table(cbind(contrasts, up_o[,1],dw_o[,1]))
DETEs_o$V2 <- as.numeric(DETEs_o$V2)
DETEs_o$V3 <- as.numeric(DETEs_o$V3)
colnames(DETEs_o)[2]<- "Upregulated"
colnames(DETEs_o)[3]<- "Downregulated"
DETEs_o_melt<- melt.data.table(DETEs_o, id.vars='contrasts')
ovaries_plot<- ggplot(DETEs_o_melt, aes(x=factor(contrasts, levels = level_order_o), y=value)) +
  geom_bar(stat='identity', aes(fill=variable))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Ovaries TE expression")+theme_bw()+
  #geom_bracket(xmin = 4, xmax = 11, label = "Hybrids", y.position = 180)+
  labs(fill="Regulation")+ geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-75,200), breaks = seq(-100, 200, 50), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

up_o_w <-rbind(acu_hec_up_o_W, hec_cin_up_o_W, acu_cin_up_o_W, acu.hec_acu_up_o_W, acu.hec_hec_up_o_W,
               acu.hecWacu_acu_up_o_W,acu.hecWacu_hec_up_o_W,hec.cin_cin_up_o_W, hec.cin_hec_up_o_W,
               hec.cinWcin_cin_up_o_W,hec.cinWcin_hec_up_o_W)
dw_o_w <-rbind(acu_hec_dw_o_W, hec_cin_dw_o_W, acu_cin_dw_o_W, acu.hec_acu_dw_o_W, acu.hec_hec_dw_o_W,
               acu.hecWacu_acu_dw_o_W,acu.hecWacu_hec_dw_o_W,hec.cin_cin_dw_o_W, hec.cin_hec_dw_o_W,
               hec.cinWcin_cin_dw_o_W,hec.cinWcin_hec_dw_o_W)
contrasts<- c("acu_hec","hec_cin","acu_cin","acu.hec_acu(hecW)","acu.hec_hec(hecW)",
              "acu.hec_acu(acuW)","acu.hec_hec(acuW)","hec.cin_cin(hecW)",
              "hec.cin_hec(hecW)","hec.cin_cin(cinW)","hec.cin_hec(cinW)")
level_order_o_W<- contrasts
DETEs_o_W <- as.data.table(cbind(contrasts, up_o_w[,1],dw_o_w[,1]))
DETEs_o_W$V2 <- as.numeric(DETEs_o_W$V2)
DETEs_o_W$V3 <- as.numeric(DETEs_o_W$V3)
colnames(DETEs_o_W)[2]<- "Upregulated"
colnames(DETEs_o_W)[3]<- "Downregulated"
DETEs_o_melt_W<- melt.data.table(DETEs_o_W, id.vars='contrasts')
DETEs_o_melt$chr <- "all"
DETEs_o_melt_W$chr <- "W"
DETEs_o_all <- rbind(DETEs_o_melt, DETEs_o_melt_W)
DETEs_o_all$chr<- as.factor(DETEs_o_all$chr)
DETEs_o_all$cat <- paste0(DETEs_o_all$variable," ",DETEs_o_all$chr)

#ovaries_plot_W<- ggplot(DETEs_o_melt_W, aes(x=factor(contrasts, levels = level_order_o), y=value)) +
geom_bar(stat='identity',aes(fill=variable))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Ovaries W chrom. TE expression")+theme_bw()+
  #geom_bracket(xmin = 4, xmax = 11, label = "Hybrids", y.position = 180)+
  labs(fill="Regulation")+ geom_hline(yintercept = 0)+  scale_y_continuous(limits = c(-75,200), breaks = seq(-100, 200, 50), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ovaries_plot_fams<- ggplot(fams_o_counts, aes(x=factor(contrast, levels = level_order_o), y=counts)) +
  geom_bar(stat='identity')+geom_col(color="black", linewidth=0.1,aes(fill=V1))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Ovaries TE superfamilies")+theme_bw()+
  #geom_bracket(xmin = 4, xmax = 11, label = "Hybrids", y.position = 180)+
  labs(fill="Superfamily")+geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-75,200), breaks = seq(-100, 200, 50), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ovaries_plot_class<- ggplot(class_o_counts, aes(x=factor(contrast, levels = level_order_o), y=counts)) +
  geom_bar(stat='identity')+geom_col(color="black", linewidth=0.1,aes(fill=V1))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Ovaries TE orders")+theme_bw()+
  #geom_bracket(xmin = 4, xmax = 11, label = "Hybrids", y.position = 180)+
  labs(fill="Order")+geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-75,200), breaks = seq(-100, 200, 50), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Testes dataset RUVseq
dge_t <- DGEList(counts = counts(dds_testes), samples = row.names(testes_meta), group = testes_meta$taxon)
dge_t <- calcNormFactors(dge_t, method = "TMM")
counts_RUV_t <- cpm(dge_t, log=F)
counts_RUV_t <- round(counts_RUV_t)

set_t<- newSeqExpressionSet(counts = counts_RUV_t, phenoData = data.frame(testes_meta$taxon, row.names= colnames(counts(dds_testes))))
keep_set_t <- rowSums(counts(set_t)) >= 10
set_t <- set_t[keep_set_t,]
x<- testes_meta$taxon
plotRLE(set_t, outline=FALSE, ylim=c(-4, 4), col=colors[set_t$testes_meta.taxon])
plotPCA(set_t, col=colors[x], cex=1.2)

#Get empirical control genes
design <- model.matrix(~x, data=pData(set_t))
y <- DGEList(counts=counts(set_t), group=x)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

top <- topTags(lrt, n=nrow(set_t))$table
empirical <- rownames(set_t)[which(!(rownames(set_t) %in% rownames(top)[1:5000]))]
res <- residuals(fit, type="deviance")

#RUV
set1_t <- RUVr(set_t, empirical, k=1, res)
plotRLE(set1_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=1")
plotPCA(set1_t, col=colors[x], main="PCA - k=1")

set2_t <- RUVr(set_t, empirical, k=2, res)
plotRLE(set2_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=2")
plotPCA(set2_t, col=colors[x], main="PCA - k=2")

set3_t <- RUVr(set_t, empirical, k=3, res)
plotRLE(set3_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=3")
plotPCA(set3_t, col=colors[x], main="PCA - k=3")

set4_t <- RUVr(set_t, empirical, k=4, res)
plotRLE(set4_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=4")
plotPCA(set4_t, col=colors[x], main="PCA - k=4")

set5_t <- RUVr(set_t, empirical, k=5, res)
plotRLE(set5_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=5")
plotPCA(set5_t, col=colors[x], main="PCA - k=5")

set6_t <- RUVr(set_t, empirical, k=6, res)
plotRLE(set6_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=6")
plotPCA(set6_t, col=colors[x], main="PCA - k=6")

set7_t <- RUVr(set_t, empirical, k=7, res)
plotRLE(set7_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=7")
plotPCA(set7_t, col=colors[x], main="PCA - k=7")

set8_t <- RUVr(set_t, empirical, k=8, res)
plotRLE(set8_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=8")
plotPCA(set8_t, col=colors[x], main="PCA - k=8")

set9_t <- RUVr(set_t, empirical, k=9, res)
plotRLE(set9_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=9")
plotPCA(set9_t, col=colors[x], main="PCA - k=9")

set10_t <- RUVr(set_t, empirical, k=10, res)
plotRLE(set10_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=10")
plotPCA(set10_t, col=colors[x], main="PCA - k=10")

set11_t <- RUVr(set_t, empirical, k=11, res)
plotRLE(set11_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=11")
plotPCA(set11_t, col=colors[x], main="PCA - k=11")

set12_t <- RUVr(set_t, empirical, k=12, res)
plotRLE(set12_t, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=12")
plotPCA(set12_t, col=colors[x], main="PCA - k=12")

#design <- model.matrix(~x+W_1+W_2+W_3+W_4+W_5+W_6+W_7+W_8+W_9+W_10+W_11, data = pData(set11))
#y <- DGEList(counts=counts(set11), group=x)
#y <- calcNormFactors(y, method="upperquartile")
#y <- estimateGLMCommonDisp(y, design)
#y <- estimateGLMTagwiseDisp(y, design)

dds_t_RUV <- DESeqDataSetFromMatrix(countData = counts(set12_t)[,1:29],
                                    colData = pData(set12_t)[1:29,],
                                    design = ~ 0+testes_meta.taxon + W_1+W_2+W_3+W_4+W_5+W_6+W_7+W_8+W_9+W_10+W_11+W_12)
dds_t_RUV <- DESeq(dds_t_RUV)
resultsNames(dds_t_RUV)
dds_t_RUV@design

contrast_hec.cin_RUV_t<-results(dds_t_RUV, contrast = c("testes_meta.taxon","acu","hec"), alpha = 0.01)
table(contrast_hec.cin_RUV_t$padj < 0.01)
contrast_hec.cin_RUV_t <- subset(contrast_hec.cin_RUV_t, padj < 0.01)
contrast_hec.cin_RUV_t_TEs<- contrast_hec.cin_RUV_t[str_detect(contrast_hec.cin_RUV_t@rownames, "RND"), ]
acu_hec_WTEs_t<- subset(contrast_hec.cin_RUV_t_TEs, unlist(strsplit(contrast_hec.cin_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin_RUV_t_TEs@rownames)])[seq(1,4*length(contrast_hec.cin_RUV_t_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin_RUV_t_TEs <- subset(contrast_hec.cin_RUV_t_TEs, !contrast_hec.cin_RUV_t_TEs %in% acu_hec_WTEs_t)
sup_f_acu_hec_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin_RUV_t_TEs@rownames)])[seq(3,4*length(contrast_hec.cin_RUV_t_TEs@rownames),4)],
                            contrast_hec.cin_RUV_t_TEs$log2FoldChange)      
class_acu_hec_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin_RUV_t_TEs@rownames)])[seq(4,4*length(contrast_hec.cin_RUV_t_TEs@rownames),4)],
                            contrast_hec.cin_RUV_t_TEs$log2FoldChange)
acu_hec_up_t<- length(which(contrast_hec.cin_RUV_t_TEs$log2FoldChange > 0))
### 3 upregulated TEs
acu_hec_dw_t<- -1*length(which(contrast_hec.cin_RUV_t_TEs$log2FoldChange < 0))
### 8 downregulated TEs
error_rate_5 <- dim(acu_hec_WTEs_t)[1]/dim(contrast_hec.cin_RUV_t_TEs)[1]

contrast_hec.cin2_RUV_t<-results(dds_t_RUV, contrast = c("testes_meta.taxon","hec","cin"), alpha = 0.01)
table(contrast_hec.cin2_RUV_t$padj < 0.01)
contrast_hec.cin2_RUV_t <- subset(contrast_hec.cin2_RUV_t, padj < 0.01)
contrast_hec.cin2_RUV_t_TEs<- contrast_hec.cin2_RUV_t[str_detect(contrast_hec.cin2_RUV_t@rownames, "RND"), ]
hec_cin_WTEs_t<- subset(contrast_hec.cin2_RUV_t_TEs, unlist(strsplit(contrast_hec.cin2_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin2_RUV_t_TEs@rownames)])[seq(1,4*length(contrast_hec.cin2_RUV_t_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin2_RUV_t_TEs <- subset(contrast_hec.cin2_RUV_t_TEs, !contrast_hec.cin2_RUV_t_TEs %in% hec_cin_WTEs_t)
sup_f_hec_cin_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin2_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin2_RUV_t_TEs@rownames)])[seq(3,4*length(contrast_hec.cin2_RUV_t_TEs@rownames),4)],
                            contrast_hec.cin2_RUV_t_TEs$log2FoldChange)      
class_hec_cin_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin2_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin2_RUV_t_TEs@rownames)])[seq(4,4*length(contrast_hec.cin2_RUV_t_TEs@rownames),4)],
                            contrast_hec.cin2_RUV_t_TEs$log2FoldChange)
hec_cin_up_t<- length(which(contrast_hec.cin2_RUV_t_TEs$log2FoldChange > 0))
### 8 upregulated TEs
hec_cin_dw_t<- -1*length(which(contrast_hec.cin2_RUV_t_TEs$log2FoldChange < 0))
### 54 downregulated TEs

error_rate_6 <- dim(hec_cin_WTEs_t)[1]/dim(contrast_hec.cin2_RUV_t_TEs)[1]

contrast_hec.cin3_RUV_t<-results(dds_t_RUV, contrast = c("testes_meta.taxon","hec.cin","cin"), alpha = 0.01)
table(contrast_hec.cin3_RUV_t$padj < 0.01)
contrast_hec.cin3_RUV_t <- subset(contrast_hec.cin3_RUV_t, padj < 0.01)
contrast_hec.cin3_RUV_t_TEs<- contrast_hec.cin3_RUV_t[str_detect(contrast_hec.cin3_RUV_t@rownames, "RND"), ]
hec.cin_cin_WTEs_t<- subset(contrast_hec.cin3_RUV_t_TEs, unlist(strsplit(contrast_hec.cin3_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin3_RUV_t_TEs@rownames)])[seq(1,4*length(contrast_hec.cin3_RUV_t_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin3_RUV_t_TEs <- subset(contrast_hec.cin3_RUV_t_TEs, !contrast_hec.cin3_RUV_t_TEs %in% hec.cin_cin_WTEs_t)
sup_f_hec.cin_cin_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin3_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin3_RUV_t_TEs@rownames)])[seq(3,4*length(contrast_hec.cin3_RUV_t_TEs@rownames),4)],
                                contrast_hec.cin3_RUV_t_TEs$log2FoldChange)      
class_hec.cin_cin_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin3_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin3_RUV_t_TEs@rownames)])[seq(4,4*length(contrast_hec.cin3_RUV_t_TEs@rownames),4)],
                                contrast_hec.cin3_RUV_t_TEs$log2FoldChange)
hec.cin_cin_up_t<- length(which(contrast_hec.cin3_RUV_t_TEs$log2FoldChange > 0))
### 198 upregulated TEs
hec.cin_cin_dw_t<- -1*length(which(contrast_hec.cin3_RUV_t_TEs$log2FoldChange < 0))
### 7 downregulated TEs

error_rate_1 <- dim(hec.cin_cin_WTEs_t)[1]/(dim(hec.cin_cin_WTEs_t)[1] + dim(contrast_hec.cin3_RUV_t_TEs)[1])


contrast_hec.cin4_RUV_t<-results(dds_t_RUV, contrast = c("testes_meta.taxon","hec.cin","hec"), alpha = 0.01)
table(contrast_hec.cin4_RUV_t$padj < 0.01)
contrast_hec.cin4_RUV_t <- subset(contrast_hec.cin4_RUV_t, padj < 0.01)
contrast_hec.cin4_RUV_t_TEs<- contrast_hec.cin4_RUV_t[str_detect(contrast_hec.cin4_RUV_t@rownames, "RND"), ]
hec.cin_hec_WTEs_t<- subset(contrast_hec.cin4_RUV_t_TEs, unlist(strsplit(contrast_hec.cin4_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin4_RUV_t_TEs@rownames)])[seq(1,4*length(contrast_hec.cin4_RUV_t_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin4_RUV_t_TEs <- subset(contrast_hec.cin4_RUV_t_TEs, !contrast_hec.cin4_RUV_t_TEs %in% hec.cin_hec_WTEs_t)
sup_f_hec.cin_hec_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin4_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin4_RUV_t_TEs@rownames)])[seq(3,4*length(contrast_hec.cin4_RUV_t_TEs@rownames),4)],
                                contrast_hec.cin4_RUV_t_TEs$log2FoldChange)      
class_hec.cin_hec_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin4_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin4_RUV_t_TEs@rownames)])[seq(4,4*length(contrast_hec.cin4_RUV_t_TEs@rownames),4)],
                                contrast_hec.cin4_RUV_t_TEs$log2FoldChange)
hec.cin_hec_up_t<- length(which(contrast_hec.cin4_RUV_t_TEs$log2FoldChange > 0))
### 251 upregulated TEs
hec.cin_hec_dw_t<- -1*length(which(contrast_hec.cin4_RUV_t_TEs$log2FoldChange < 0))
### 5 downregulated TEs

error_rate_2 <- dim(hec.cin_hec_WTEs_t)[1]/(dim(hec.cin_hec_WTEs_t)[1] + dim(contrast_hec.cin4_RUV_t_TEs)[1])

contrast_hec.cin5_RUV_t<-results(dds_t_RUV, contrast = c("testes_meta.taxon","acu.hec","hec"), alpha = 0.01)
table(contrast_hec.cin5_RUV_t$padj < 0.01)
contrast_hec.cin5_RUV_t <- subset(contrast_hec.cin5_RUV_t, padj < 0.01)
contrast_hec.cin5_RUV_t_TEs<- contrast_hec.cin5_RUV_t[str_detect(contrast_hec.cin5_RUV_t@rownames, "RND"), ]
acu.hec_hec_WTEs_t<- subset(contrast_hec.cin5_RUV_t_TEs, unlist(strsplit(contrast_hec.cin5_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin5_RUV_t_TEs@rownames)])[seq(1,4*length(contrast_hec.cin5_RUV_t_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin5_RUV_t_TEs <- subset(contrast_hec.cin5_RUV_t_TEs, !contrast_hec.cin5_RUV_t_TEs %in% acu.hec_hec_WTEs_t)
sup_f_acu.hec_hec_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin5_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin5_RUV_t_TEs@rownames)])[seq(3,4*length(contrast_hec.cin5_RUV_t_TEs@rownames),4)],
                                contrast_hec.cin5_RUV_t_TEs$log2FoldChange)      
class_acu.hec_hec_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin5_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin5_RUV_t_TEs@rownames)])[seq(4,4*length(contrast_hec.cin5_RUV_t_TEs@rownames),4)],
                                contrast_hec.cin5_RUV_t_TEs$log2FoldChange)
acu.hec_hec_up_t<- length(which(contrast_hec.cin5_RUV_t_TEs$log2FoldChange > 0))
### 247 upregulated TEs
acu.hec_hec_dw_t<- -1*length(which(contrast_hec.cin5_RUV_t_TEs$log2FoldChange < 0))
### 12 downregulated TEs

error_rate_3 <- dim(acu.hec_hec_WTEs_t)[1]/(dim(acu.hec_hec_WTEs_t)[1] + dim(contrast_hec.cin5_RUV_t_TEs)[1])

contrast_hec.cin6_RUV_t<-results(dds_t_RUV, contrast = c("testes_meta.taxon","acu.hec","acu"), alpha = 0.01)
table(contrast_hec.cin6_RUV_t$padj < 0.01)
contrast_hec.cin6_RUV_t <- subset(contrast_hec.cin6_RUV_t, padj < 0.01)
contrast_hec.cin6_RUV_t_TEs<- contrast_hec.cin6_RUV_t[str_detect(contrast_hec.cin6_RUV_t@rownames, "RND"), ]
acu.hec_acu_WTEs_t<- subset(contrast_hec.cin6_RUV_t_TEs, unlist(strsplit(contrast_hec.cin6_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin6_RUV_t_TEs@rownames)])[seq(1,4*length(contrast_hec.cin6_RUV_t_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin6_RUV_t_TEs <- subset(contrast_hec.cin6_RUV_t_TEs, !contrast_hec.cin6_RUV_t_TEs %in% acu.hec_acu_WTEs_t)
sup_f_acu.hec_acu_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin6_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin6_RUV_t_TEs@rownames)])[seq(3,4*length(contrast_hec.cin6_RUV_t_TEs@rownames),4)],
                                contrast_hec.cin6_RUV_t_TEs$log2FoldChange)      
class_acu.hec_acu_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin6_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin6_RUV_t_TEs@rownames)])[seq(4,4*length(contrast_hec.cin6_RUV_t_TEs@rownames),4)],
                                contrast_hec.cin6_RUV_t_TEs$log2FoldChange)
acu.hec_acu_up_t<- length(which(contrast_hec.cin6_RUV_t_TEs$log2FoldChange > 0))
### 266 upregulated TEs
acu.hec_acu_dw_t<- -1*length(which(contrast_hec.cin6_RUV_t_TEs$log2FoldChange < 0))
### 18 downregulated TEs

error_rate_4 <- dim(acu.hec_acu_WTEs_t)[1]/(dim(acu.hec_acu_WTEs_t)[1] + dim(contrast_hec.cin6_RUV_t_TEs)[1])

contrast_hec.cin7_RUV_t<-results(dds_t_RUV, contrast = c("testes_meta.taxon","acu","cin"), alpha = 0.01)
table(contrast_hec.cin7_RUV_t$padj < 0.01)
contrast_hec.cin7_RUV_t <- subset(contrast_hec.cin7_RUV_t, padj < 0.01)
contrast_hec.cin7_RUV_t_TEs<- contrast_hec.cin7_RUV_t[str_detect(contrast_hec.cin7_RUV_t@rownames, "RND"), ]
acu_cin_WTEs_t<- subset(contrast_hec.cin7_RUV_t_TEs, unlist(strsplit(contrast_hec.cin7_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin7_RUV_t_TEs@rownames)])[seq(1,4*length(contrast_hec.cin7_RUV_t_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin7_RUV_t_TEs <- subset(contrast_hec.cin7_RUV_t_TEs, !contrast_hec.cin7_RUV_t_TEs %in% acu_cin_WTEs_t)
sup_f_acu_cin_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin7_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin7_RUV_t_TEs@rownames)])[seq(3,4*length(contrast_hec.cin7_RUV_t_TEs@rownames),4)],
                            contrast_hec.cin7_RUV_t_TEs$log2FoldChange)      
class_acu_cin_TEs_t<- cbind(unlist(strsplit(contrast_hec.cin7_RUV_t_TEs@rownames, ":")[1:length(contrast_hec.cin7_RUV_t_TEs@rownames)])[seq(4,4*length(contrast_hec.cin7_RUV_t_TEs@rownames),4)],
                            contrast_hec.cin7_RUV_t_TEs$log2FoldChange)
acu_cin_up_t<- length(which(contrast_hec.cin7_RUV_t_TEs$log2FoldChange > 0))
### 7 upregulated TEs
acu_cin_dw_t<- -1*length(which(contrast_hec.cin7_RUV_t_TEs$log2FoldChange < 0))
### 76 downregulated TEs

error_rate_7 <- dim(acu_cin_WTEs_t)[1]/dim(contrast_hec.cin7_RUV_t_TEs)[1]

mean_err_t_W <- mean(c(error_rate_1,error_rate_2,error_rate_3,error_rate_4,error_rate_5,
                       error_rate_6,error_rate_7))

## Plot DETEs for testes
fams_t <- list(sup_f_acu_hec_TEs_t,sup_f_hec_cin_TEs_t, sup_f_acu_cin_TEs_t, sup_f_acu.hec_acu_TEs_t,
               sup_f_acu.hec_hec_TEs_t,sup_f_hec.cin_cin_TEs_t,sup_f_hec.cin_hec_TEs_t)
fams_t <- lapply(fams_t, as.data.frame)
fams_t_up<- lapply(fams_t, subset, V2 > 0)
fams_t_dw<- lapply(fams_t, subset, V2 < 0)
fams_t_up <- bind_rows(fams_t_up, .id = "df_id")
fams_t_dw <- bind_rows(fams_t_dw, .id = "df_id")

class_t <- list(class_acu_hec_TEs_t,class_hec_cin_TEs_t, class_acu_cin_TEs_t, class_acu.hec_acu_TEs_t,
                class_acu.hec_hec_TEs_t,class_hec.cin_cin_TEs_t,class_hec.cin_hec_TEs_t)
class_t <- lapply(class_t, as.data.frame)
class_t_up<- lapply(class_t, subset, V2 > 0)
class_t_dw<- lapply(class_t, subset, V2 < 0)
class_t_up <- bind_rows(class_t_up, .id = "df_id")
class_t_dw <- bind_rows(class_t_dw, .id = "df_id")

up_t <-rbind(acu_hec_up_t, hec_cin_up_t, acu_cin_up_t, acu.hec_acu_up_t, acu.hec_hec_up_t,
             hec.cin_cin_up_t, hec.cin_hec_up_t)
dw_t <-rbind(acu_hec_dw_t, hec_cin_dw_t, acu_cin_dw_t, acu.hec_acu_dw_t, acu.hec_hec_dw_t,
             hec.cin_cin_dw_t, hec.cin_hec_dw_t)
contrasts<- c("acu_hec","hec_cin","acu_cin","acu.hec_acu","acu.hec_hec","hec.cin_cin",
              "hec.cin_hec")
level_order_t<- contrasts

for(i in 1:length(fams_t_up$df_id)){fams_t_up$contrast[i]<-contrasts[as.numeric(fams_t_up$df_id[i])]}
for(i in 1:length(fams_t_dw$df_id)){fams_t_dw$contrast[i]<-contrasts[as.numeric(fams_t_dw$df_id[i])]}
for(i in 1:length(class_t_up$df_id)){class_t_up$contrast[i]<-contrasts[as.numeric(class_t_up$df_id[i])]}
for(i in 1:length(class_t_dw$df_id)){class_t_dw$contrast[i]<-contrasts[as.numeric(class_t_dw$df_id[i])]}
fams_t_up_counts <- fams_t_up %>% group_by(contrast, V1) %>% summarise(counts=n())
fams_t_dw_counts <- fams_t_dw %>% group_by(contrast, V1) %>% summarise(counts=n()*-1)
fams_t_counts <- rbind(fams_t_up_counts,fams_t_dw_counts)

class_t_up_counts <- class_t_up %>% group_by(contrast, V1) %>% summarise(counts=n())
class_t_dw_counts <- class_t_dw %>% group_by(contrast, V1) %>% summarise(counts=n()*-1)
class_t_counts <- rbind(class_t_up_counts,class_t_dw_counts)
DETEs_t <- as.data.table(cbind(contrasts, up_t[,1],dw_t[,1]))
DETEs_t$V2 <- as.numeric(DETEs_t$V2)
DETEs_t$V3 <- as.numeric(DETEs_t$V3)
colnames(DETEs_t)[2]<- "Upregulated"
colnames(DETEs_t)[3]<- "Downregulated"
DETEs_t_melt<- melt.data.table(DETEs_t, id.vars='contrasts')
testes_plot<- ggplot(DETEs_t_melt, aes(x=factor(contrasts, levels = level_order_t), y=value)) +
  geom_bar(stat='identity', aes(fill=variable))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Testes TE expression")+theme_bw()+
  #geom_bracket(xmin = 4, xmax = 7, label = "Hybrids", y.position = 300)+
  labs(fill="Regulation")+geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-350,900), breaks = seq(-400, 1000, 200), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

testes_plot_fams<- ggplot(fams_t_counts, aes(x=factor(contrast, levels = level_order_t), y=counts)) +
  geom_bar(stat='identity')+geom_col(color="black", linewidth=0.1,aes(fill=V1))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Testes TE superfamilies")+theme_bw()+
  #geom_bracket(xmin = 4, xmax = 11, label = "Hybrids", y.position = 180)+
  labs(fill="Superfamily")+geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-350,900), breaks = seq(-400, 1000, 200), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

testes_plot_class<- ggplot(class_t_counts, aes(x=factor(contrast, levels = level_order_t), y=counts)) +
  geom_bar(stat='identity')+geom_col(color="black", linewidth=0.1,aes(fill=V1))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Testes TE orders")+theme_bw()+
  #geom_bracket(xmin = 4, xmax = 11, label = "Hybrids", y.position = 180)+
  labs(fill="Order")+geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-350,900), breaks = seq(-400, 1000, 200), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#### Liver data RUV#####
dge_l <- DGEList(counts = counts(dds_liver), samples = row.names(liver_meta), group = liver_meta$group)
dge_l <- calcNormFactors(dge_l, method = "TMM")
counts_RUV_l <- cpm(dge_l, log=F)
counts_RUV_l <- round(counts_RUV_l)

set_l<- newSeqExpressionSet(counts = counts_RUV_l, phenoData = data.frame(liver_meta$group, row.names= colnames(counts(dds_liver))))
keep_set_l <- rowSums(counts(set_l)) >= 10
set_l <- set_l[keep_set_l,]
x<- liver_meta$group
plotRLE(set_l, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set_l, col=colors[x], cex=1.2)

#Get empirical control genes
design <- model.matrix(~x, data=pData(set_l))
y <- DGEList(counts=counts(set_l), group=x)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)

top <- topTags(lrt, n=nrow(set_l))$table
empirical <- rownames(set_l)[which(!(rownames(set_l) %in% rownames(top)[1:5000]))]
res <- residuals(fit, type="deviance")

#RUV
set1_l <- RUVr(set_l, empirical, k=1, res)
plotRLE(set1_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=1")
plotPCA(set1_l, col=colors[x], main="PCA - k=1")

set2_l <- RUVr(set_l, empirical, k=2, res)
plotRLE(set2_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=2")
plotPCA(set2_l, col=colors[x], main="PCA - k=2")

set3_l <- RUVr(set_l, empirical, k=3, res)
plotRLE(set3_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=3")
plotPCA(set3_l, col=colors[x], main="PCA - k=3")

set4_l <- RUVr(set_l, empirical, k=4, res)
plotRLE(set4_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=4")
plotPCA(set4_l, col=colors[x], main="PCA - k=4")

set5_l <- RUVr(set_l, empirical, k=5, res)
plotRLE(set5_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=5")
plotPCA(set5_l, col=colors[x], main="PCA - k=5")

set6_l <- RUVr(set_l, empirical, k=6, res)
plotRLE(set6_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=6")
plotPCA(set6_l, col=colors[x], main="PCA - k=6")

set7_l <- RUVr(set_l, empirical, k=7, res)
plotRLE(set7_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=7")
plotPCA(set7_l, col=colors[x], main="PCA - k=7")

set8_l <- RUVr(set_l, empirical, k=8, res)
plotRLE(set8_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=8")
plotPCA(set8_l, col=colors[x], main="PCA - k=8")

set9_l <- RUVr(set_l, empirical, k=9, res)
plotRLE(set9_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=9")
plotPCA(set9_l, col=colors[x], main="PCA - k=9")

set10_l <- RUVr(set_l, empirical, k=10, res)
plotRLE(set10_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=10")
plotPCA(set10_l, col=colors[x], main="PCA - k=10")

set11_l <- RUVr(set_l, empirical, k=11, res)
plotRLE(set11_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=11")
plotPCA(set11_l, col=colors[x], main="PCA - k=11")

set12_l <- RUVr(set_l, empirical, k=12, res)
plotRLE(set12_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=12")
plotPCA(set12_l, col=colors[x], main="PCA - k=12")

set13_l <- RUVr(set_l, empirical, k=13, res)
plotRLE(set13_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=13")
plotPCA(set13_l, col=colors[x], main="PCA - k=13")

set14_l <- RUVr(set_l, empirical, k=14, res)
plotRLE(set14_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=14")
plotPCA(set14_l, col=colors[x], main="PCA - k=14")

set15_l <- RUVr(set_l, empirical, k=15, res)
plotRLE(set15_l, outline=F, ylim = c(-2, 2), col = colors[x], main="RLE - k=15")
plotPCA(set15_l, col=colors[x], main="PCA - k=15")

#design <- model.matrix(~x+W_1+W_2+W_3+W_4+W_5+W_6+W_7+W_8+W_9+W_10+W_11, data = pData(set11))
#y <- DGEList(counts=counts(set11), group=x)
#y <- calcNormFactors(y, method="upperquartile")
#y <- estimateGLMCommonDisp(y, design)
#y <- estimateGLMTagwiseDisp(y, design)

dds_l_RUV <- DESeqDataSetFromMatrix(countData = counts(set12_l),
                                    colData = pData(set12_l),
                                    design = ~ 0+liver_meta.group + W_1+W_2+W_3+W_4+W_5+W_6+W_7+W_8+W_9+W_10+W_11+W_12)
dds_l_RUV <- DESeq(dds_l_RUV)
resultsNames(dds_l_RUV)

contrast_hec.cin_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","1_acu","2_hec"), alpha = 0.01)
table(contrast_hec.cin_RUV_l$padj < 0.01)
contrast_hec.cin_RUV_l <- subset(contrast_hec.cin_RUV_l, padj < 0.01)
contrast_hec.cin_RUV_l_TEs<- contrast_hec.cin_RUV_l[str_detect(contrast_hec.cin_RUV_l@rownames, "RND"), ]
#acu_hec_WTEs_l<- subset(contrast_hec.cin_RUV_l_TEs, unlist(strsplit(contrast_hec.cin_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin_RUV_l_TEs@rownames),3)] %in% W_chrTEs$V13)
sup_f_acu_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin_RUV_l_TEs$log2FoldChange)      
class_acu_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin_RUV_l_TEs$log2FoldChange)

acu_hec_F_up_l<- length(which(contrast_hec.cin_RUV_l_TEs$log2FoldChange > 0))
### 0 upregulated TEs
acu_hec_F_dw_l<- -1*length(which(contrast_hec.cin_RUV_l_TEs$log2FoldChange < 0))
### 2 downregulated TEs
acu_hec_FW_up_l<- length(which(acu_hec_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
acu_hec_FW_dw_l<- -1*length(which(acu_hec_WTEs_l$log2FoldChange < 0))

contrast_hec.cin2_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","6","7"), alpha = 0.01)
table(contrast_hec.cin2_RUV_l$padj < 0.01)
contrast_hec.cin2_RUV_l <- subset(contrast_hec.cin2_RUV_l, padj < 0.01)
contrast_hec.cin2_RUV_l_TEs<- contrast_hec.cin2_RUV_l[str_detect(contrast_hec.cin2_RUV_l@rownames, "RND"), ]
#acu_hec_M_WTEs_l<- subset(contrast_hec.cin2_RUV_l_TEs, unlist(strsplit(contrast_hec.cin2_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin2_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin2_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin2_RUV_l_TEs <- subset(contrast_hec.cin2_RUV_l_TEs, !contrast_hec.cin2_RUV_l_TEs %in% acu_hec_M_WTEs_l)
sup_f_acu_hec_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin2_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin2_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin2_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin2_RUV_l_TEs$log2FoldChange)      
class_acu_hec_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin2_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin2_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin2_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin2_RUV_l_TEs$log2FoldChange)
acu_hec_M_up_l<- length(which(contrast_hec.cin2_RUV_l_TEs$log2FoldChange > 0))
### 0 upregulated TEs
acu_hec_M_dw_l<- -1*length(which(contrast_hec.cin2_RUV_l_TEs$log2FoldChange < 0))
### 4 downregulated TEs

#error_rate_l_1 <- dim(acu_hec_M_WTEs_l)[1]/dim(contrast_hec.cin2_RUV_l_TEs)[1]

contrast_hec.cin3_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","5","7"), alpha = 0.01)
table(contrast_hec.cin3_RUV_l$padj < 0.01)
contrast_hec.cin3_RUV_l <- subset(contrast_hec.cin3_RUV_l, padj < 0.01)
contrast_hec.cin3_RUV_l_TEs<- contrast_hec.cin3_RUV_l[str_detect(contrast_hec.cin3_RUV_l@rownames, "RND"), ]
#hec.cin_hec_M_WTEs_l<- subset(contrast_hec.cin3_RUV_l_TEs, unlist(strsplit(contrast_hec.cin3_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin3_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin3_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin3_RUV_l_TEs <- subset(contrast_hec.cin3_RUV_l_TEs, !contrast_hec.cin3_RUV_l_TEs %in% hec.cin_hec_M_WTEs_l)
sup_f_hec.cin_hec_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin3_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin3_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin3_RUV_l_TEs@rownames),3)],
                                  contrast_hec.cin3_RUV_l_TEs$log2FoldChange)      
class_hec.cin_hec_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin3_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin3_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin3_RUV_l_TEs@rownames),3)],
                                  contrast_hec.cin3_RUV_l_TEs$log2FoldChange)
hec.cin_hec_M_up_l<- length(which(contrast_hec.cin3_RUV_l_TEs$log2FoldChange > 0))
### 219 upregulated TEs
hec.cin_hec_M_dw_l<- -1*length(which(contrast_hec.cin3_RUV_l_TEs$log2FoldChange < 0))
### 8 downregulated TEs

#error_rate_l_2 <- dim(hec.cin_hec_M_WTEs_l)[1]/dim(contrast_hec.cin3_RUV_l_TEs)[1]

contrast_hec.cin4_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","10_hec","2_hec"), alpha = 0.01)
table(contrast_hec.cin4_RUV_l$padj < 0.01)
contrast_hec.cin4_RUV_l <- subset(contrast_hec.cin4_RUV_l, padj < 0.01)
contrast_hec.cin4_RUV_l_TEs<- contrast_hec.cin4_RUV_l[str_detect(contrast_hec.cin4_RUV_l@rownames, "RND"), ]
#hec.cinWhec_hec_WTEs_l<- subset(contrast_hec.cin4_RUV_l_TEs, unlist(strsplit(contrast_hec.cin4_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin4_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin4_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec.cinWhec_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin4_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin4_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin4_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin4_RUV_l_TEs$log2FoldChange)      
class_hec.cinWhec_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin4_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin4_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin4_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin4_RUV_l_TEs$log2FoldChange)
hec.cinWhec_hec_F_up_l<-length(which(contrast_hec.cin4_RUV_l_TEs$log2FoldChange > 0))
### 227 upregulated TEs
hec.cinWhec_hec_F_dw_l<- -1*length(which(contrast_hec.cin4_RUV_l_TEs$log2FoldChange < 0))
### 5 downregulated TEs
hec.cinWhec_hec_FW_up_l<- length(which(hec.cinWhec_hec_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
hec.cinWhec_hec_FW_dw_l<- -1*length(which(hec.cinWhec_hec_WTEs_l$log2FoldChange < 0))

contrast_hec.cin5_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","5","8"), alpha = 0.01)
table(contrast_hec.cin5_RUV_l$padj < 0.01)
contrast_hec.cin5_RUV_l <- subset(contrast_hec.cin5_RUV_l, padj < 0.01)
contrast_hec.cin5_RUV_l_TEs<- contrast_hec.cin5_RUV_l[str_detect(contrast_hec.cin5_RUV_l@rownames, "RND"), ]
#hec.cin_cin_M_WTEs_l<- subset(contrast_hec.cin5_RUV_l_TEs, unlist(strsplit(contrast_hec.cin5_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin5_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin5_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin5_RUV_l_TEs[174,] <- subset(contrast_hec.cin5_RUV_l_TEs, !contrast_hec.cin5_RUV_l_TEs %in% hec.cin_cin_M_WTEs_l)
sup_f_hec.cin_cin_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin5_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin5_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin5_RUV_l_TEs@rownames),3)],
                                  contrast_hec.cin5_RUV_l_TEs$log2FoldChange)      
class_hec.cin_cin_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin5_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin5_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin5_RUV_l_TEs@rownames),3)],
                                  contrast_hec.cin5_RUV_l_TEs$log2FoldChange)
hec.cin_cin_M_up_l<-length(which(contrast_hec.cin5_RUV_l_TEs$log2FoldChange > 0))
### 197 upregulated TEs
hec.cin_cin_M_dw_l<- -1*length(which(contrast_hec.cin5_RUV_l_TEs$log2FoldChange < 0))
### 3 downregulated TEs

#error_rate_l_3 <- dim(hec.cin_cin_M_WTEs_l)[1]/dim(contrast_hec.cin5_RUV_l_TEs)[1]

contrast_hec.cin6_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","10_hec","3_cin"), alpha = 0.01)
table(contrast_hec.cin6_RUV_l$padj < 0.01)
contrast_hec.cin6_RUV_l <- subset(contrast_hec.cin6_RUV_l, padj < 0.01)
contrast_hec.cin6_RUV_l_TEs<- contrast_hec.cin6_RUV_l[str_detect(contrast_hec.cin6_RUV_l@rownames, "RND"), ]
#hec.cinWhec_cin_WTEs_l<- subset(contrast_hec.cin6_RUV_l_TEs, unlist(strsplit(contrast_hec.cin6_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin6_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin6_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec.cinWhec_cin_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin6_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin6_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin6_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin6_RUV_l_TEs$log2FoldChange)      
class_hec.cinWhec_cin_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin6_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin6_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin6_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin6_RUV_l_TEs$log2FoldChange)
hec.cinWhec_cin_F_up_l<-length(which(contrast_hec.cin6_RUV_l_TEs$log2FoldChange > 0))
### 248 upregulated TEs
hec.cinWhec_cin_F_dw_l<- -1*length(which(contrast_hec.cin6_RUV_l_TEs$log2FoldChange < 0))
### 5 downregulated TEs
hec.cinWhec_cin_FW_up_l<- length(which(hec.cinWhec_cin_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
hec.cinWhec_cin_FW_dw_l<- -1*length(which(hec.cinWhec_cin_WTEs_l$log2FoldChange < 0))


contrast_hec.cin7_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","9","7"), alpha = 0.01)
table(contrast_hec.cin7_RUV_l$padj < 0.01)
contrast_hec.cin7_RUV_l <- subset(contrast_hec.cin7_RUV_l, padj < 0.01)
contrast_hec.cin7_RUV_l_TEs<- contrast_hec.cin7_RUV_l[str_detect(contrast_hec.cin7_RUV_l@rownames, "RND"), ]
acu.hec_hec_M_WTEs_l<- subset(contrast_hec.cin7_RUV_l_TEs, unlist(strsplit(contrast_hec.cin7_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin7_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin7_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin7_RUV_l_TEs <- subset(contrast_hec.cin7_RUV_l_TEs, !contrast_hec.cin7_RUV_l_TEs %in% acu.hec_hec_M_WTEs_l)
sup_f_acu.hec_hec_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin7_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin7_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin7_RUV_l_TEs@rownames),3)],
                                  contrast_hec.cin7_RUV_l_TEs$log2FoldChange)      
class_acu.hec_hec_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin7_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin7_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin7_RUV_l_TEs@rownames),3)],
                                  contrast_hec.cin7_RUV_l_TEs$log2FoldChange)
acu.hec_hec_M_up_l<- length(which(contrast_hec.cin7_RUV_l_TEs$log2FoldChange > 0))
### 121 upregulated TEs
acu.hec_hec_M_dw_l<- -1*length(which(contrast_hec.cin7_RUV_l_TEs$log2FoldChange < 0))
### 24 downregulated TEs

#error_rate_l_4 <- dim(acu.hec_hec_M_WTEs_l)[1]/dim(contrast_hec.cin7_RUV_l_TEs)[1]

contrast_hec.cin8_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","4_acu","2_hec"), alpha = 0.01)
table(contrast_hec.cin8_RUV_l$padj < 0.01)
contrast_hec.cin8_RUV_l <- subset(contrast_hec.cin8_RUV_l, padj < 0.01)
contrast_hec.cin8_RUV_l_TEs<- contrast_hec.cin8_RUV_l[str_detect(contrast_hec.cin8_RUV_l@rownames, "RND"), ]
acu.hecWacu_hec_WTEs_l<- subset(contrast_hec.cin8_RUV_l_TEs, unlist(strsplit(contrast_hec.cin8_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin8_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin8_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu.hecWacu_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin8_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin8_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin8_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin8_RUV_l_TEs$log2FoldChange)      
class_acu.hecWacu_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin8_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin8_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin8_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin8_RUV_l_TEs$log2FoldChange)
acu.hecWacu_hec_F_up_l<-length(which(contrast_hec.cin8_RUV_l_TEs$log2FoldChange > 0))
### 95 upregulated TEs
acu.hecWacu_hec_F_dw_l<- -1*length(which(contrast_hec.cin8_RUV_l_TEs$log2FoldChange < 0))
### 29 downregulated TEs
acu.hecWacu_hec_FW_up_l<- length(which(acu.hecWacu_hec_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
acu.hecWacu_hec_FW_dw_l<- -1*length(which(acu.hecWacu_hec_WTEs_l$log2FoldChange < 0))

contrast_hec.cin9_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","9","6"), alpha = 0.01)
table(contrast_hec.cin9_RUV_l$padj < 0.01)
contrast_hec.cin9_RUV_l <- subset(contrast_hec.cin9_RUV_l, padj < 0.01)
contrast_hec.cin9_RUV_l_TEs<- contrast_hec.cin9_RUV_l[str_detect(contrast_hec.cin9_RUV_l@rownames, "RND"), ]
acu.hec_acu_M_WTEs_l<- subset(contrast_hec.cin9_RUV_l_TEs, unlist(strsplit(contrast_hec.cin9_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin9_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin9_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin9_RUV_l_TEs <- subset(contrast_hec.cin9_RUV_l_TEs, !contrast_hec.cin9_RUV_l_TEs %in% acu.hec_acu_M_WTEs_l)
sup_f_acu.hec_acu_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin9_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin9_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin9_RUV_l_TEs@rownames),3)],
                                  contrast_hec.cin9_RUV_l_TEs$log2FoldChange)      
class_acu.hec_acu_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin9_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin9_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin9_RUV_l_TEs@rownames),3)],
                                  contrast_hec.cin9_RUV_l_TEs$log2FoldChange)
acu.hec_acu_M_up_l<-length(which(contrast_hec.cin9_RUV_l_TEs$log2FoldChange > 0))
### 117 upregulated TEs
acu.hec_acu_M_dw_l<- -1*length(which(contrast_hec.cin9_RUV_l_TEs$log2FoldChange < 0))
### 17 downregulated TEs

#error_rate_l_5 <- dim(acu.hec_acu_M_WTEs_l)[1]/dim(contrast_hec.cin9_RUV_l_TEs)[1]

contrast_hec.cin10_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","4_acu","1_acu"), alpha = 0.01)
table(contrast_hec.cin10_RUV_l$padj < 0.01)
contrast_hec.cin10_RUV_l <- subset(contrast_hec.cin10_RUV_l, padj < 0.01)
contrast_hec.cin10_RUV_l_TEs<- contrast_hec.cin10_RUV_l[str_detect(contrast_hec.cin10_RUV_l@rownames, "RND"), ]
acu.hecWacu_acu_WTEs_l<- subset(contrast_hec.cin10_RUV_l_TEs, unlist(strsplit(contrast_hec.cin10_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin10_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin10_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu.hecWacu_acu_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin10_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin10_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin10_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin10_RUV_l_TEs$log2FoldChange)      
class_acu.hecWacu_acu_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin10_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin10_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin10_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin10_RUV_l_TEs$log2FoldChange)
acu.hecWacu_acu_F_up_l<- length(which(contrast_hec.cin10_RUV_l_TEs$log2FoldChange > 0))
### 66 upregulated TEs
acu.hecWacu_acu_F_dw_l<- -1*length(which(contrast_hec.cin10_RUV_l_TEs$log2FoldChange < 0))
### 17 downregulated TEs
acu.hecWacu_acu_FW_up_l<- length(which(acu.hecWacu_acu_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
acu.hecWacu_acu_FW_dw_l<- -1*length(which(acu.hecWacu_acu_WTEs_l$log2FoldChange < 0))

contrast_hec.cin11_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","2_hec","3_cin"), alpha = 0.01)
table(contrast_hec.cin11_RUV_l$padj < 0.01)
contrast_hec.cin11_RUV_l <- subset(contrast_hec.cin11_RUV_l, padj < 0.01)
contrast_hec.cin11_RUV_l_TEs<- contrast_hec.cin11_RUV_l[str_detect(contrast_hec.cin11_RUV_l@rownames, "RND"), ]
hec_cin_WTEs_l<- subset(contrast_hec.cin11_RUV_l_TEs, unlist(strsplit(contrast_hec.cin11_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin11_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin11_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec_cin_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin11_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin11_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin11_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin11_RUV_l_TEs$log2FoldChange)      
class_hec_cin_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin11_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin11_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin11_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin11_RUV_l_TEs$log2FoldChange)
hec_cin_F_up_l<- length(which(contrast_hec.cin11_RUV_l_TEs$log2FoldChange > 0))
### 7 upregulated TEs
hec_cin_F_dw_l<- -1*length(which(contrast_hec.cin11_RUV_l_TEs$log2FoldChange < 0))
### 4 downregulated TEs
hec_cin_FW_up_l<- length(which(hec_cin_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
hec_cin_FW_dw_l<- -1*length(which(hec_cin_WTEs_l$log2FoldChange < 0))

contrast_hec.cin12_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","7","8"), alpha = 0.01)
table(contrast_hec.cin12_RUV_l$padj < 0.01)
contrast_hec.cin12_RUV_l <- subset(contrast_hec.cin12_RUV_l, padj < 0.01)
contrast_hec.cin12_RUV_l_TEs<- contrast_hec.cin12_RUV_l[str_detect(contrast_hec.cin12_RUV_l@rownames, "RND"), ]
hec_cin_M_WTEs_l<- subset(contrast_hec.cin12_RUV_l_TEs, unlist(strsplit(contrast_hec.cin12_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin12_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin12_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin12_RUV_l_TEs <- subset(contrast_hec.cin12_RUV_l_TEs, !contrast_hec.cin12_RUV_l_TEs %in% hec_cin_M_WTEs_l)
sup_f_hec_cin_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin12_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin12_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin12_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin12_RUV_l_TEs$log2FoldChange)      
class_hec_cin_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin12_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin12_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin12_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin12_RUV_l_TEs$log2FoldChange)
hec_cin_M_up_l<-length(which(contrast_hec.cin12_RUV_l_TEs$log2FoldChange > 0))
### 9 upregulated TEs
hec_cin_M_dw_l<- -1*length(which(contrast_hec.cin12_RUV_l_TEs$log2FoldChange < 0))
### 9 downregulated TEs

#error_rate_l_6 <- dim(hec_cin_M_WTEs_l)[1]/dim(contrast_hec.cin12_RUV_l_TEs)[1]

contrast_hec.cin13_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","1_acu","3_cin"), alpha = 0.01)
table(contrast_hec.cin13_RUV_l$padj < 0.01)
contrast_hec.cin13_RUV_l <- subset(contrast_hec.cin13_RUV_l, padj < 0.01)
contrast_hec.cin13_RUV_l_TEs<- contrast_hec.cin13_RUV_l[str_detect(contrast_hec.cin13_RUV_l@rownames, "RND"), ]
acu_cin_WTEs_l<- subset(contrast_hec.cin13_RUV_l_TEs, unlist(strsplit(contrast_hec.cin13_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin13_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin13_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu_cin_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin13_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin13_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin13_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin13_RUV_l_TEs$log2FoldChange)      
class_acu_cin_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin13_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin13_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin13_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin13_RUV_l_TEs$log2FoldChange)
acu_cin_F_up_l<- length(which(contrast_hec.cin13_RUV_l_TEs$log2FoldChange > 0))
### 7 upregulated TEs
acu_cin_F_dw_l<- -1*length(which(contrast_hec.cin13_RUV_l_TEs$log2FoldChange < 0))
### 11 downregulated TEs
acu_cin_FW_up_l<- length(which(acu_cin_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
acu_cin_FW_dw_l<- -1*length(which(acu_cin_WTEs_l$log2FoldChange < 0))

contrast_hec.cin14_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","6","8"), alpha = 0.01)
table(contrast_hec.cin14_RUV_l$padj < 0.01)
contrast_hec.cin14_RUV_l <- subset(contrast_hec.cin14_RUV_l, padj < 0.01)
contrast_hec.cin14_RUV_l_TEs<- contrast_hec.cin14_RUV_l[str_detect(contrast_hec.cin14_RUV_l@rownames, "RND"), ]
acu_cin_M_WTEs_l<- subset(contrast_hec.cin14_RUV_l_TEs, unlist(strsplit(contrast_hec.cin14_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin14_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin14_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
contrast_hec.cin14_RUV_l_TEs <- subset(contrast_hec.cin14_RUV_l_TEs, !contrast_hec.cin14_RUV_l_TEs %in% acu_cin_M_WTEs_l)
sup_f_acu_cin_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin14_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin14_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin14_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin14_RUV_l_TEs$log2FoldChange)      
class_acu_cin_M_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin14_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin14_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin14_RUV_l_TEs@rownames),3)],
                              contrast_hec.cin14_RUV_l_TEs$log2FoldChange)
acu_cin_M_up_l<- length(which(contrast_hec.cin14_RUV_l_TEs$log2FoldChange > 0))
### 5 upregulated TEs
acu_cin_M_dw_l<- -1*length(which(contrast_hec.cin14_RUV_l_TEs$log2FoldChange < 0))
### 13 downregulated TEs

#error_rate_l_7 <- dim(acu_cin_M_WTEs_l)[1]/dim(contrast_hec.cin14_RUV_l_TEs)[1]

contrast_hec.cin15_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","10_cin","2_hec"), alpha = 0.01)
table(contrast_hec.cin15_RUV_l$padj < 0.01)
contrast_hec.cin15_RUV_l <- subset(contrast_hec.cin15_RUV_l, padj < 0.01)
contrast_hec.cin15_RUV_l_TEs<- contrast_hec.cin15_RUV_l[str_detect(contrast_hec.cin15_RUV_l@rownames, "RND"), ]
hec.cinWcin_hec_WTEs_l<- subset(contrast_hec.cin15_RUV_l_TEs, unlist(strsplit(contrast_hec.cin15_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin15_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin15_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec.cinWcin_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin15_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin15_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin15_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin15_RUV_l_TEs$log2FoldChange)      
class_hec.cinWcin_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin15_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin15_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin15_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin15_RUV_l_TEs$log2FoldChange)
hec.cinWcin_hec_F_up_l<-length(which(contrast_hec.cin15_RUV_l_TEs$log2FoldChange > 0))
### 194 upregulated TEs
hec.cinWcin_hec_F_dw_l<- -1*length(which(contrast_hec.cin15_RUV_l_TEs$log2FoldChange < 0))
### 5 downregulated TEs
hec.cinWcin_hec_FW_up_l<- length(which(hec.cinWcin_hec_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
hec.cinWcin_hec_FW_dw_l<- -1*length(which(hec.cinWcin_hec_WTEs_l$log2FoldChange < 0))

contrast_hec.cin16_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","10_cin","3_cin"), alpha = 0.01)
table(contrast_hec.cin16_RUV_l$padj < 0.01)
contrast_hec.cin16_RUV_l <- subset(contrast_hec.cin16_RUV_l, padj < 0.01)
contrast_hec.cin16_RUV_l_TEs<- contrast_hec.cin16_RUV_l[str_detect(contrast_hec.cin16_RUV_l@rownames, "RND"), ]
hec.cinWcin_cin_WTEs_l<- subset(contrast_hec.cin16_RUV_l_TEs, unlist(strsplit(contrast_hec.cin16_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin16_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin16_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_hec.cinWcin_cin_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin16_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin16_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin16_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin16_RUV_l_TEs$log2FoldChange)      
class_hec.cinWcin_cin_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin16_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin16_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin16_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin16_RUV_l_TEs$log2FoldChange)
hec.cinWcin_cin_F_up_l<-length(which(contrast_hec.cin16_RUV_l_TEs$log2FoldChange > 0))
### 209 upregulated TEs
hec.cinWcin_cin_F_dw_l<- -1*length(which(contrast_hec.cin16_RUV_l_TEs$log2FoldChange < 0))
### 4 downregulated TEs
hec.cinWcin_cin_FW_up_l<- length(which(hec.cinWcin_cin_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
hec.cinWcin_cin_FW_dw_l<- -1*length(which(hec.cinWcin_cin_WTEs_l$log2FoldChange < 0))

contrast_hec.cin17_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","4_hec","2_hec"), alpha = 0.01)
table(contrast_hec.cin17_RUV_l$padj < 0.01)
contrast_hec.cin17_RUV_l <- subset(contrast_hec.cin17_RUV_l, padj < 0.01)
contrast_hec.cin17_RUV_l_TEs<- contrast_hec.cin17_RUV_l[str_detect(contrast_hec.cin17_RUV_l@rownames, "RND"), ]
acu.hecWhec_hec_WTEs_l<- subset(contrast_hec.cin17_RUV_l_TEs, unlist(strsplit(contrast_hec.cin17_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin17_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin17_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu.hecWhec_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin17_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin17_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin17_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin17_RUV_l_TEs$log2FoldChange)      
class_acu.hecWhec_hec_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin17_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin17_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin17_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin17_RUV_l_TEs$log2FoldChange)
acu.hecWhec_hec_F_up_l<-length(which(contrast_hec.cin17_RUV_l_TEs$log2FoldChange > 0))
### 332 upregulated TEs
acu.hecWhec_hec_F_dw_l<- -1*length(which(contrast_hec.cin17_RUV_l_TEs$log2FoldChange < 0))
### 72 downregulated TEs
acu.hecWhec_hec_FW_up_l<- length(which(acu.hecWhec_hec_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
acu.hecWhec_hec_FW_dw_l<- -1*length(which(acu.hecWhec_hec_WTEs_l$log2FoldChange < 0))

contrast_hec.cin18_RUV_l<-results(dds_l_RUV, contrast = c("liver_meta.group","4_hec","1_acu"), alpha = 0.01)
table(contrast_hec.cin18_RUV_l$padj < 0.01)
contrast_hec.cin18_RUV_l <- subset(contrast_hec.cin18_RUV_l, padj < 0.01)
contrast_hec.cin18_RUV_l_TEs<- contrast_hec.cin18_RUV_l[str_detect(contrast_hec.cin18_RUV_l@rownames, "RND"), ]
acu.hecWhec_acu_WTEs_l<- subset(contrast_hec.cin18_RUV_l_TEs, unlist(strsplit(contrast_hec.cin18_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin18_RUV_l_TEs@rownames)])[seq(1,4*length(contrast_hec.cin18_RUV_l_TEs@rownames),4)] %in% W_chrTEs$V13)
sup_f_acu.hecWhec_acu_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin18_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin18_RUV_l_TEs@rownames)])[seq(2,3*length(contrast_hec.cin18_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin18_RUV_l_TEs$log2FoldChange)      
class_acu.hecWhec_acu_F_TEs_l<- cbind(unlist(strsplit(contrast_hec.cin18_RUV_l_TEs@rownames, ":")[1:length(contrast_hec.cin18_RUV_l_TEs@rownames)])[seq(3,3*length(contrast_hec.cin18_RUV_l_TEs@rownames),3)],
                                      contrast_hec.cin18_RUV_l_TEs$log2FoldChange)
acu.hecWhec_acu_F_up_l<- length(which(contrast_hec.cin18_RUV_l_TEs$log2FoldChange > 0))
### 56 upregulated TEs
acu.hecWhec_acu_F_dw_l<- -1*length(which(contrast_hec.cin18_RUV_l_TEs$log2FoldChange < 0))
### 16 downregulated TEs
acu.hecWhec_acu_FW_up_l<- length(which(acu.hecWhec_acu_WTEs_l$log2FoldChange > 0))
### 0 upregulated TEs
acu.hecWhec_acu_FW_dw_l<- -1*length(which(acu.hecWhec_acu_WTEs_l$log2FoldChange < 0))

#mean_err_l_W <- mean(c(error_rate_l_1,error_rate_l_2,error_rate_l_3,error_rate_l_4,error_rate_l_5,
#                     error_rate_l_6,error_rate_l_7))

## Plot DETEs for liver ####

up_l <-rbind(acu_hec_F_up_l, acu_hec_M_up_l, hec_cin_F_up_l, hec_cin_M_up_l, acu_cin_F_up_l,
             acu_cin_M_up_l, acu.hecWacu_acu_F_up_l,acu.hecWhec_acu_F_up_l,acu.hec_acu_M_up_l,
             acu.hecWacu_hec_F_up_l,acu.hecWhec_hec_F_up_l, acu.hec_hec_M_up_l,
             hec.cinWhec_cin_F_up_l, hec.cinWcin_cin_F_up_l, hec.cin_cin_M_up_l,
             hec.cinWhec_hec_F_up_l,hec.cinWcin_hec_F_up_l, hec.cin_hec_M_up_l)
dw_l <-rbind(acu_hec_F_dw_l, acu_hec_M_dw_l, hec_cin_F_dw_l, hec_cin_M_dw_l, acu_cin_F_dw_l,
             acu_cin_M_dw_l, acu.hecWacu_acu_F_dw_l,acu.hecWhec_acu_F_dw_l,acu.hec_acu_M_dw_l,
             acu.hecWacu_hec_F_dw_l,acu.hecWhec_hec_F_dw_l, acu.hec_hec_M_dw_l,
             hec.cinWhec_cin_F_dw_l, hec.cinWcin_cin_F_dw_l, hec.cin_cin_M_dw_l,
             hec.cinWhec_hec_F_dw_l,hec.cinWcin_hec_F_dw_l, hec.cin_hec_M_dw_l)
contrasts_l<- c("acu_hec","acu_hec","hec_cin","hec_cin","acu_cin","acu_cin",
                "acu.hec(acuW)_acu","acu.hec(hecW)_acu","acu.hec_acu","acu.hec(acuW)_hec",
                "acu.hec(hecW)_hec","acu.hec_hec","hec.cin(hecW)_cin","hec.cin(cinW)_cin",
                "hec.cin_cin", "hec.cin(hecW)_hec","hec.cin(cinW)_hec","hec.cin_hec")
sex_l<- c("F","M","F","M","F","M","F","F","M","F","F","M","F","F","M","F","F","M")
level_order_l<- paste0(contrasts_l,sex_l)
DETEs_l <- as.data.table(cbind(contrasts_l, sex_l, up_l[,1],dw_l[,1]))
DETEs_l$V3 <- as.numeric(DETEs_l$V3)
DETEs_l$V4 <- as.numeric(DETEs_l$V4)
colnames(DETEs_l)[3]<- "Upregulated"
colnames(DETEs_l)[4]<- "Downregulated"
DETEs_l_melt<- melt.data.table(DETEs_l, id.vars=c('contrasts_l','sex_l'))
DETEs_l_melt$group <- paste0(DETEs_l_melt$contrasts_l, DETEs_l_melt$sex_l)

liver_plot<- ggplot(DETEs_l_melt, aes(x=factor(group, levels = level_order_l), y=value)) +
  geom_bar(stat='identity', aes(fill=variable))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Liver TE expression")+theme_bw()+
  #geom_bracket(xmin = 7, xmax = 8, label = "Hybrids", y.position = 1150)+
  facet_wrap(~sex_l,scale="free_x")+labs(fill="Regulation")+geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-50,250), breaks = seq(-50, 300, 50), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

up_l_W <-rbind(acu_hec_FW_up_l, hec_cin_FW_up_l, acu_cin_FW_up_l,
               acu.hecWacu_acu_FW_up_l,acu.hecWhec_acu_FW_up_l,
               acu.hecWacu_hec_FW_up_l,acu.hecWhec_hec_FW_up_l,
               hec.cinWhec_cin_FW_up_l, hec.cinWcin_cin_FW_up_l,
               hec.cinWhec_hec_FW_up_l,hec.cinWcin_hec_FW_up_l)
dw_l_W <-rbind(acu_hec_FW_dw_l, hec_cin_FW_dw_l, acu_cin_FW_dw_l,
               acu.hecWacu_acu_FW_dw_l,acu.hecWhec_acu_FW_dw_l,
               acu.hecWacu_hec_FW_dw_l,acu.hecWhec_hec_FW_dw_l,
               hec.cinWhec_cin_FW_dw_l, hec.cinWcin_cin_FW_dw_l,
               hec.cinWhec_hec_FW_dw_l,hec.cinWcin_hec_FW_dw_l)
contrasts_l_w<- c("acu_hec","hec_cin","acu_cin",
                  "acu.hec(acuW)_acu","acu.hec(hecW)_acu","acu.hec(acuW)_hec",
                  "acu.hec(hecW)_hec","hec.cin(hecW)_cin","hec.cin(cinW)_cin",
                  "hec.cin(hecW)_hec","hec.cin(cinW)_hec")
sex_l_w<- c("F","F","F","F","F","F","F","F","F","F","F")
level_order_l_w<- paste0(contrasts_l_w,sex_l_w)
DETEs_l_w <- as.data.table(cbind(contrasts_l_w, sex_l_w, up_l_W[,1],dw_l_W[,1]))
DETEs_l_w$V3 <- as.numeric(DETEs_l_w$V3)
DETEs_l_w$V4 <- as.numeric(DETEs_l_w$V4)
colnames(DETEs_l_w)[3]<- "Upregulated"
colnames(DETEs_l_w)[4]<- "Downregulated"
DETEs_l_melt_w<- melt.data.table(DETEs_l_w, id.vars=c('contrasts_l_w','sex_l_w'))
DETEs_l_melt_w$group <- paste0(DETEs_l_melt_w$contrasts_l_w, DETEs_l_melt_w$sex_l_w)

fams_l <- list(sup_f_acu_hec_F_TEs_l, sup_f_acu_hec_M_TEs_l, sup_f_hec_cin_F_TEs_l, sup_f_hec_cin_M_TEs_l, sup_f_acu_cin_F_TEs_l,
               sup_f_acu_cin_M_TEs_l, sup_f_acu.hecWacu_acu_F_TEs_l,sup_f_acu.hecWhec_acu_F_TEs_l,sup_f_acu.hec_acu_M_TEs_l,
               sup_f_acu.hecWacu_hec_F_TEs_l,sup_f_acu.hecWhec_hec_F_TEs_l, sup_f_acu.hec_hec_M_TEs_l,
               sup_f_hec.cinWhec_cin_F_TEs_l, sup_f_hec.cinWcin_cin_F_TEs_l, sup_f_hec.cin_cin_M_TEs_l,
               sup_f_hec.cinWhec_hec_F_TEs_l,sup_f_hec.cinWcin_hec_F_TEs_l, sup_f_hec.cin_hec_M_TEs_l)
fams_l <- lapply(fams_l, as.data.frame)
fams_l_up<- lapply(fams_l, subset, V2 > 0)
fams_l_dw<- lapply(fams_l, subset, V2 < 0)
fams_l_up <- bind_rows(fams_l_up, .id = "df_id")
fams_l_dw <- bind_rows(fams_l_dw, .id = "df_id")

class_l <- list(class_acu_hec_F_TEs_l, class_acu_hec_M_TEs_l, class_hec_cin_F_TEs_l, class_hec_cin_M_TEs_l, class_acu_cin_F_TEs_l,
                class_acu_cin_M_TEs_l, class_acu.hecWacu_acu_F_TEs_l,class_acu.hecWhec_acu_F_TEs_l,class_acu.hec_acu_M_TEs_l,
                class_acu.hecWacu_hec_F_TEs_l,class_acu.hecWhec_hec_F_TEs_l, class_acu.hec_hec_M_TEs_l,
                class_hec.cinWhec_cin_F_TEs_l, class_hec.cinWcin_cin_F_TEs_l, class_hec.cin_cin_M_TEs_l,
                class_hec.cinWhec_hec_F_TEs_l,class_hec.cinWcin_hec_F_TEs_l, class_hec.cin_hec_M_TEs_l)
class_l <- lapply(class_l, as.data.frame)
class_l_up<- lapply(class_l, subset, V2 > 0)
class_l_dw<- lapply(class_l, subset, V2 < 0)
class_l_up <- bind_rows(class_l_up, .id = "df_id")
class_l_dw <- bind_rows(class_l_dw, .id = "df_id")

for(i in 1:length(fams_l_up$df_id)){fams_l_up$contrast_l[i]<-level_order_l[as.numeric(fams_l_up$df_id[i])]}
for(i in 1:length(fams_l_dw$df_id)){fams_l_dw$contrast_l[i]<-level_order_l[as.numeric(fams_l_dw$df_id[i])]}
for(i in 1:length(class_l_up$df_id)){class_l_up$contrast_l[i]<-level_order_l[as.numeric(class_l_up$df_id[i])]}
for(i in 1:length(class_l_dw$df_id)){class_l_dw$contrast_l[i]<-level_order_l[as.numeric(class_l_dw$df_id[i])]}
for(i in 1:length(fams_l_up$df_id)){fams_l_up$sex_l[i]<-sex_l[as.numeric(fams_l_up$df_id[i])]}
for(i in 1:length(fams_l_dw$df_id)){fams_l_dw$sex_l[i]<-sex_l[as.numeric(fams_l_dw$df_id[i])]}
for(i in 1:length(class_l_up$df_id)){class_l_up$sex_l[i]<-sex_l[as.numeric(class_l_up$df_id[i])]}
for(i in 1:length(class_l_dw$df_id)){class_l_dw$sex_l[i]<-sex_l[as.numeric(class_l_dw$df_id[i])]}
fams_l_up_counts <- fams_l_up %>% group_by(contrast_l, V1) %>% summarise(counts=n(), sex_l=first(sex_l))
fams_l_dw_counts <- fams_l_dw %>% group_by(contrast_l, V1) %>% summarise(counts=n()*-1, sex_l=first(sex_l))
fams_l_counts <- rbind(fams_l_up_counts,fams_l_dw_counts)

class_l_up_counts <- class_l_up %>% group_by(contrast_l, V1) %>% summarise(counts=n(), sex_l=first(sex_l))
class_l_dw_counts <- class_l_dw %>% group_by(contrast_l, V1) %>% summarise(counts=n()*-1, sex_l=first(sex_l))
class_l_counts <- rbind(class_l_up_counts,class_l_dw_counts)

liver_plot_fams<- ggplot(fams_l_counts, aes(x=factor(contrast_l, levels = level_order_l), y=counts)) +
  geom_bar(stat='identity')+geom_col(color="black", linewidth=0.1,aes(fill=V1))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Liver TE superfamilies")+theme_bw()+
  #geom_bracket(xmin = 4, xmax = 11, label = "Hybrids", y.position = 180)+
  labs(fill="Superfamily")+geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-50,250), breaks = seq(-50, 300, 50), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~sex_l,scale="free_x")

liver_plot_class<- ggplot(class_l_counts, aes(x=factor(contrast_l, levels = level_order_l), y=counts)) +
  geom_bar(stat='identity')+geom_col(color="black", linewidth=0.1,aes(fill=V1))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Liver TE classes")+theme_bw()+
  #geom_bracket(xmin = 4, xmax = 11, label = "Hybrids", y.position = 180)+
  labs(fill="Order")+geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-50,250), breaks = seq(-50, 300, 50), labels = abs)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(~sex_l,scale="free_x")

#liver_plot_w<- ggplot(DETEs_l_melt_w, aes(x=factor(group, levels = level_order_l_w), y=value)) +
#  geom_bar(stat='identity', aes(fill=variable))+xlab("Species contrasts")+ylab("DETEs")+ggtitle("Liver W chrom. TE expression")+theme_bw()+
  #geom_bracket(xmin = 7, xmax = 8, label = "Hybrids", y.position = 1150)+
#  labs(fill="Regulation")+geom_hline(yintercept = 0)+ scale_y_continuous(limits = c(-200,650), breaks = seq(-200, 800, 200), labels = abs)+
 # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ovaries_plot
testes_plot
liver_plot

### Produce Venn diagrams ####
library(ggVennDiagram)
### acu-hec ovaries TEs and W ###
ov_acu.hec_m_up <- list("AxH(HW)H" = row.names(subset(contrast_hec.cin5_RUV_o_TEs, contrast_hec.cin5_RUV_o_TEs$log2FoldChange > 0)),
                        "AxH(HW)A"= row.names(subset(contrast_hec.cin6_RUV_o_TEs, contrast_hec.cin6_RUV_o_TEs$log2FoldChange > 0)),
                        "AxH(AW)H"=row.names(subset(contrast_hec.cin10_RUV_o_TEs, contrast_hec.cin10_RUV_o_TEs$log2FoldChange > 0)),
                        "AxH(AW)A"=row.names(subset(contrast_hec.cin11_RUV_o_TEs, contrast_hec.cin11_RUV_o_TEs$log2FoldChange > 0)))
o_venn_ah_up<- ggVennDiagram(ov_acu.hec_m_up, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "blue")

ov_acu.hec_m_dw <- list("AxH(HW)H" = row.names(subset(contrast_hec.cin5_RUV_o_TEs, contrast_hec.cin5_RUV_o_TEs$log2FoldChange < 0)),
                        "AxH(HW)A"= row.names(subset(contrast_hec.cin6_RUV_o_TEs, contrast_hec.cin6_RUV_o_TEs$log2FoldChange < 0)),
                        "AxH(AW)H"=row.names(subset(contrast_hec.cin10_RUV_o_TEs, contrast_hec.cin10_RUV_o_TEs$log2FoldChange < 0)),
                        "AxH(AW)A"=row.names(subset(contrast_hec.cin11_RUV_o_TEs, contrast_hec.cin11_RUV_o_TEs$log2FoldChange < 0)))
o_venn_ah_dw<-ggVennDiagram(ov_acu.hec_m_dw, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "red")


### hec-cin ovaries TEs and W
ov_hec.cin_m_up <- list("HxC(HW)H" = row.names(subset(contrast_hec.cin3_RUV_o_TEs, contrast_hec.cin3_RUV_o_TEs$log2FoldChange > 0)),
                        "HxC(HW)C"= row.names(subset(contrast_hec.cin4_RUV_o_TEs, contrast_hec.cin4_RUV_o_TEs$log2FoldChange > 0)),
                        "HxC(CW)H"=row.names(subset(contrast_hec.cin8_RUV_o_TEs, contrast_hec.cin8_RUV_o_TEs$log2FoldChange > 0)),
                        "HxC(CW)C"=row.names(subset(contrast_hec.cin9_RUV_o_TEs, contrast_hec.cin9_RUV_o_TEs$log2FoldChange > 0)))
o_venn_hc_up<-ggVennDiagram(ov_hec.cin_m_up, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "blue")

ov_hec.cin_m_dw <- list("HxC(HW)H" = row.names(subset(contrast_hec.cin3_RUV_o_TEs, contrast_hec.cin3_RUV_o_TEs$log2FoldChange < 0)),
                        "HxC(HW)C"= row.names(subset(contrast_hec.cin4_RUV_o_TEs, contrast_hec.cin4_RUV_o_TEs$log2FoldChange < 0)),
                        "HxC(CW)H"=row.names(subset(contrast_hec.cin8_RUV_o_TEs, contrast_hec.cin8_RUV_o_TEs$log2FoldChange < 0)),
                        "HxC(CW)C"=row.names(subset(contrast_hec.cin9_RUV_o_TEs, contrast_hec.cin9_RUV_o_TEs$log2FoldChange < 0)))
o_venn_hc_dw<-ggVennDiagram(ov_hec.cin_m_dw, label_alpha = 0.4,set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "red")

library(ggpubr)

o_venn_ah_all <- ggarrange(o_venn_ah_up,
                           o_venn_hc_up,
                           o_venn_ah_dw,
                           o_venn_hc_dw,
                           ncol = 2, nrow = 2, labels = c("A","B","C","D"))

ggsave("venn_ovaries_all_fams.png",plot=o_venn_ah_all, dpi = 'retina', height = 8, width = 12)

### Testes ####

### acu-hec testes ###
t_acu.hec_m_up <- list("AxH-H" = row.names(subset(contrast_hec.cin5_RUV_t_TEs, contrast_hec.cin5_RUV_t_TEs$log2FoldChange > 0)),
                       "AxH-A"= row.names(subset(contrast_hec.cin6_RUV_t_TEs, contrast_hec.cin6_RUV_t_TEs$log2FoldChange > 0)))
t_venn_ah_up<- ggVennDiagram(t_acu.hec_m_up, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "blue")

t_acu.hec_m_dw <- list("AxH-H" = row.names(subset(contrast_hec.cin5_RUV_t_TEs, contrast_hec.cin5_RUV_t_TEs$log2FoldChange < 0)),
                       "AxH-A"= row.names(subset(contrast_hec.cin6_RUV_t_TEs, contrast_hec.cin6_RUV_t_TEs$log2FoldChange < 0)))
t_venn_ah_dw<-ggVennDiagram(t_acu.hec_m_dw, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "red")

t_hec.cin_m_up <- list("HxC-H" = row.names(subset(contrast_hec.cin4_RUV_t_TEs, contrast_hec.cin4_RUV_t_TEs$log2FoldChange > 0)),
                       "HxC-C"= row.names(subset(contrast_hec.cin3_RUV_t_TEs, contrast_hec.cin3_RUV_t_TEs$log2FoldChange > 0)))
t_venn_hc_up<- ggVennDiagram(t_hec.cin_m_up, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "blue")

t_hec.cin_m_dw <- list("HxC-H" = row.names(subset(contrast_hec.cin4_RUV_t_TEs, contrast_hec.cin4_RUV_t_TEs$log2FoldChange < 0)),
                       "HxC-C"= row.names(subset(contrast_hec.cin3_RUV_t_TEs, contrast_hec.cin3_RUV_t_TEs$log2FoldChange < 0)))
t_venn_hc_dw<-ggVennDiagram(t_hec.cin_m_dw, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "red")

t_venn_ah_all <- ggarrange(t_venn_ah_up,t_venn_hc_up,t_venn_ah_dw,t_venn_hc_dw,
                           ncol = 2, nrow = 2, labels = c("A","B","C","D"))
t_venn_hc_all <- ggarrange(t_venn_hc_up,t_venn_hc_dw,
                           ncol = 1, nrow = 2, labels = c("A","B"))

ggsave("venn_testes_all_fams.png",plot=t_venn_ah_all, dpi = 'retina', height = 8, width = 12)

#### Liver #####

### Male liver ###
l_acu.hec_m_up <- list("AxH-H" = row.names(subset(contrast_hec.cin7_RUV_l_TEs, contrast_hec.cin7_RUV_l_TEs$log2FoldChange > 0)),
                       "AxH-A"= row.names(subset(contrast_hec.cin9_RUV_l_TEs, contrast_hec.cin9_RUV_l_TEs$log2FoldChange > 0)))
l_venn_ah_up<- ggVennDiagram(l_acu.hec_m_up, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "blue")

l_acu.hec_m_dw <- list("AxH-H" = row.names(subset(contrast_hec.cin7_RUV_l_TEs, contrast_hec.cin7_RUV_l_TEs$log2FoldChange < 0)),
                       "AxH-A"= row.names(subset(contrast_hec.cin9_RUV_l_TEs, contrast_hec.cin9_RUV_l_TEs$log2FoldChange < 0)))
l_venn_ah_dw<-ggVennDiagram(l_acu.hec_m_dw, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "red")

l_hec.cin_m_up <- list("HxC-H" = row.names(subset(contrast_hec.cin3_RUV_l_TEs, contrast_hec.cin3_RUV_l_TEs$log2FoldChange > 0)),
                       "HxC-C"= row.names(subset(contrast_hec.cin5_RUV_l_TEs, contrast_hec.cin5_RUV_l_TEs$log2FoldChange > 0)))
l_venn_hc_up<- ggVennDiagram(l_hec.cin_m_up, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "blue")

l_hec.cin_m_dw <- list("HxC-H" = row.names(subset(contrast_hec.cin3_RUV_l_TEs, contrast_hec.cin3_RUV_l_TEs$log2FoldChange < 0)),
                       "HxC-C"= row.names(subset(contrast_hec.cin5_RUV_l_TEs, contrast_hec.cin5_RUV_l_TEs$log2FoldChange < 0)))
l_venn_hc_dw<-ggVennDiagram(l_hec.cin_m_dw, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "red")

l_venn_ah_all_M <- ggarrange(l_venn_ah_up,l_venn_hc_up,l_venn_ah_dw,l_venn_hc_dw,
                             ncol = 2, nrow = 2, labels = c("A","B","C","D"))

ggsave("fams_venn_liver_all_M.png",plot=l_venn_ah_all_M, dpi = 'retina', height = 8, width = 12)

#### Female liver ####
### acu-hec liver TEs and W ###
l_acu.hec_m_up_F <- list("AxH(HW)H" = row.names(subset(contrast_hec.cin17_RUV_l_TEs, contrast_hec.cin17_RUV_l_TEs$log2FoldChange > 0)),
                         "AxH(HW)A"= row.names(subset(contrast_hec.cin18_RUV_l_TEs, contrast_hec.cin18_RUV_l_TEs$log2FoldChange > 0)),
                         "AxH(AW)H"=row.names(subset(contrast_hec.cin8_RUV_l_TEs, contrast_hec.cin8_RUV_l_TEs$log2FoldChange > 0)),
                         "AxH(AW)A"=row.names(subset(contrast_hec.cin10_RUV_l_TEs, contrast_hec.cin10_RUV_l_TEs$log2FoldChange > 0)))
l_venn_ah_up_F<- ggVennDiagram(l_acu.hec_m_up_F, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "blue")

l_acu.hec_m_dw_F <- list("AxH(HW)H" = row.names(subset(contrast_hec.cin17_RUV_l_TEs, contrast_hec.cin17_RUV_l_TEs$log2FoldChange < 0)),
                         "AxH(HW)A"= row.names(subset(contrast_hec.cin18_RUV_l_TEs, contrast_hec.cin18_RUV_l_TEs$log2FoldChange < 0)),
                         "AxH(AW)H"=row.names(subset(contrast_hec.cin8_RUV_l_TEs, contrast_hec.cin8_RUV_l_TEs$log2FoldChange < 0)),
                         "AxH(AW)A"=row.names(subset(contrast_hec.cin10_RUV_l_TEs, contrast_hec.cin10_RUV_l_TEs$log2FoldChange < 0)))
l_venn_ah_dw_F<-ggVennDiagram(l_acu.hec_m_dw_F, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "red")


### hec-cin ovaries TEs and W
l_hec.cin_m_up_F <- list("HxC(HW)H" = row.names(subset(contrast_hec.cin4_RUV_l_TEs, contrast_hec.cin4_RUV_l_TEs$log2FoldChange > 0)),
                         "HxC(HW)C"= row.names(subset(contrast_hec.cin6_RUV_l_TEs, contrast_hec.cin6_RUV_l_TEs$log2FoldChange > 0)),
                         "HxC(CW)H"=row.names(subset(contrast_hec.cin15_RUV_l_TEs, contrast_hec.cin15_RUV_l_TEs$log2FoldChange > 0)),
                         "HxC(CW)C"=row.names(subset(contrast_hec.cin16_RUV_l_TEs, contrast_hec.cin16_RUV_l_TEs$log2FoldChange > 0)))
l_venn_hc_up_F<-ggVennDiagram(l_hec.cin_m_up_F, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "blue")

l_hec.cin_m_dw_F <- list("HxC(HW)H" = row.names(subset(contrast_hec.cin4_RUV_l_TEs, contrast_hec.cin4_RUV_l_TEs$log2FoldChange < 0)),
                         "HxC(HW)C"= row.names(subset(contrast_hec.cin6_RUV_l_TEs, contrast_hec.cin6_RUV_l_TEs$log2FoldChange < 0)),
                         "HxC(CW)H"=row.names(subset(contrast_hec.cin15_RUV_l_TEs, contrast_hec.cin15_RUV_l_TEs$log2FoldChange < 0)),
                         "HxC(CW)C"=row.names(subset(contrast_hec.cin16_RUV_l_TEs, contrast_hec.cin16_RUV_l_TEs$log2FoldChange < 0)))
l_venn_hc_dw_F<-ggVennDiagram(l_hec.cin_m_dw_F, label_alpha = 0.4, set_size = 2.5)+ggplot2::scale_fill_gradient(low="white",high = "red")

ibrary(ggpubr)

l_venn_ah_all_F <- ggarrange(l_venn_ah_up_F,l_venn_hc_up_F,l_venn_ah_dw_F,l_venn_hc_dw_F,
                             ncol = 2, nrow = 2, labels = c("A","B","C","D"))
ggsave("fams_venn_liver_all_F.png",plot=l_venn_ah_all_F, dpi = 'retina', height = 8, width = 12)

ovaries_fig<- ggarrange(ovaries_plot,
                        ovaries_plot_class, 
                        ovaries_plot_fams,
                        ncol = 2, nrow = 2, labels = c("A","B","C"))
testes_fig<- ggarrange(testes_plot, testes_plot_class, testes_plot_fams,
                       ncol = 1, nrow = 3, labels = c("A","B","C"))
liver_fig<- ggarrange(liver_plot, 
                      liver_plot_class, liver_plot_fams,
                      ncol = 1, nrow = 3, labels = c("A","B","C"))

ggsave("ovaries_fig.png",plot=ovaries_fig, dpi = 'retina', height = 8, width = 12)
ggsave("testes_fig.png",plot=testes_fig, dpi = 'retina', height = 15, width = 8)
ggsave("liver_fig_fams.png",plot=liver_fig, dpi = 'retina', height = 15, width = 9)

all_o_list <- list(ov_acu.hec_m_up,ov_acu.hec_m_dw,ov_hec.cin_m_up,ov_hec.cin_m_dw) 
o_heatmap_d<- subset(counts(dds_o_RUV), str_detect(row.names(counts(dds_o_RUV)),"RND"))
o_heatmap_d<- subset(o_heatmap_d, row.names(o_heatmap_d) %in% unlist(all_o_list))
my_sample_col <- data.frame(Taxon = ovaries_meta$taxon)
row.names(my_sample_col)<-colnames(counts(dds_o_RUV))
ov_hm<-pheatmap(o_heatmap_d, scale = "row", annotation_col = my_sample_col, show_rownames = F)
dim(o_heatmap_d)

### Testes ###
all_t_list <- list(t_acu.hec_m_up,t_acu.hec_m_dw,t_hec.cin_m_up,t_hec.cin_m_dw)
t_heatmap_d<- subset(counts(dds_t_RUV), str_detect(row.names(counts(dds_t_RUV)),"RND"))
t_heatmap_d<- subset(t_heatmap_d, row.names(t_heatmap_d) %in% unlist(all_t_list))
my_sample_col_t <- data.frame(Taxon = testes_meta$taxon)
row.names(my_sample_col_t)<-colnames(counts(dds_t_RUV))
ts_hm<-pheatmap(t_heatmap_d, scale = "row", annotation_col = my_sample_col_t, show_rownames = F)
dim(t_heatmap_d)

### Liver ####
all_l_list <- list(l_acu.hec_m_up_F,l_acu.hec_m_dw_F,l_hec.cin_m_up_F,l_hec.cin_m_dw_F,
                   l_acu.hec_m_up,l_acu.hec_m_dw,l_hec.cin_m_up,l_hec.cin_m_dw)
l_heatmap_d<- subset(counts(dds_l_RUV), str_detect(row.names(counts(dds_l_RUV)),"RND"))
l_heatmap_d<- subset(l_heatmap_d, row.names(l_heatmap_d) %in% unlist(all_l_list))
my_sample_col_l <- data.frame(Taxon = liver_meta$taxon)
row.names(my_sample_col_l)<-colnames(counts(dds_l_RUV))
my_sample_col_l$Sex <- liver_meta$sex
lv_hm<-pheatmap(l_heatmap_d, scale = "row", annotation_col = my_sample_col_l, show_rownames = F)
dim(l_heatmap_d)

ggsave("liver_heatmap.png",plot=lv_hm, dpi = 'retina', height = 15, width = 18)
ggsave("testes_heatmap.png",plot=ts_hm, dpi = 'retina', height = 15, width = 18)
ggsave("ovaries_heatmap.png",plot=ov_hm, dpi = 'retina', height = 15, width = 18)      

### Boxplots per class ####
class_o_up_counts$State <- NA
class_o_up_counts$State[which(class_o_up_counts$contrast=="acu_hec")]<- "Non-hybrid"
class_o_up_counts$State[which(class_o_up_counts$contrast=="acu_cin")]<- "Non-hybrid"
class_o_up_counts$State[which(class_o_up_counts$contrast=="hec_cin")]<- "Non-hybrid"
class_o_up_counts$State[is.na(class_o_up_counts$State)]<- "Hybrid"
class_o_up_counts$State <- as.factor(class_o_up_counts$State)
class_o_up_counts<- subset(class_o_up_counts, class_o_up_counts$V1 != "DNA" & class_o_up_counts$V1 != "SINE" & class_o_up_counts$V1 != "Satellite")

class_boxplot_o<- ggplot(class_o_up_counts, aes(x=V1, y=counts))+geom_boxplot(aes(color=State))+
  stat_compare_means(aes(group = State), method = "t.test", label = "p.signif")+
  scale_y_continuous(limits = c(0,120), breaks = seq(0,120,20))+
  geom_point(aes(color=State),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6, jitter.height = 0), alpha = 0.5)+
  xlab("Order of DETEs")+ylab("DETEs per contrast")+ggtitle("Ovaries")+theme_bw()

class_t_up_counts$State <- NA
class_t_up_counts$State[which(class_t_up_counts$contrast=="acu_hec")]<- "Non-hybrid"
class_t_up_counts$State[which(class_t_up_counts$contrast=="acu_cin")]<- "Non-hybrid"
class_t_up_counts$State[which(class_t_up_counts$contrast=="hec_cin")]<- "Non-hybrid"
class_t_up_counts$State[is.na(class_t_up_counts$State)]<- "Hybrid"
class_t_up_counts$State <- as.factor(class_t_up_counts$State)
class_t_up_counts<- subset(class_t_up_counts, class_t_up_counts$V1 != "DNA" & class_t_up_counts$V1 != "SINE" & class_t_up_counts$V1 != "Satellite")

class_boxplot_t<-ggplot(class_t_up_counts, aes(x=V1, y=counts))+geom_boxplot(aes(color=State))+
  stat_compare_means(aes(group = State), method = "t.test", label = "p.signif")+
  scale_y_continuous(limits = c(0,120), breaks = seq(0,120,20))+
  geom_point(aes(color=State),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6, jitter.height = 0), alpha = 0.5)+
  xlab("Order of DETEs")+ylab("DETEs per contrast")+ggtitle("Testes")+theme_bw()

class_l_up_counts$State <- NA
class_l_up_counts$State[which(class_l_up_counts$contrast_l=="acu_hecM")]<- "Non-hybrid"
class_l_up_counts$State[which(class_l_up_counts$contrast_l=="acu_cinM")]<- "Non-hybrid"
class_l_up_counts$State[which(class_l_up_counts$contrast_l=="hec_cinM")]<- "Non-hybrid"
class_l_up_counts$State[which(class_l_up_counts$contrast_l=="acu_hecF")]<- "Non-hybrid"
class_l_up_counts$State[which(class_l_up_counts$contrast_l=="acu_cinF")]<- "Non-hybrid"
class_l_up_counts$State[which(class_l_up_counts$contrast_l=="hec_cinF")]<- "Non-hybrid"
class_l_up_counts$State[is.na(class_l_up_counts$State)]<- "Hybrid"
class_l_up_counts$State <- as.factor(class_l_up_counts$State)

class_l_up_counts<- subset(class_l_up_counts, class_l_up_counts$V1 != "DNA" & class_l_up_counts$V1 != "SINE" & class_l_up_counts$V1 != "Satellite")
missing_l <- data.frame(contrast_l="hec_cinF",V1="Unknown",counts=0,sex_l="F",State="Non-hybrid")
class_l_up_counts<- rbind(class_l_up_counts,missing_l) 

class_boxplot_l<-ggplot(class_l_up_counts, aes(x=V1, y=counts))+geom_boxplot(aes(color=State))+
  stat_compare_means(aes(group = State), method = "t.test", label = "p.signif")+
  scale_y_continuous(limits = c(0,120), breaks = seq(0,120,20))+
  geom_point(aes(color=State),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6, jitter.height = 0), alpha = 0.5)+
  xlab("Order of DETEs")+ylab("DETEs per contrast")+ggtitle("Liver")+facet_wrap(~sex_l)+theme_bw()


ggsave("order_bxplot_o.png",plot=class_boxplot_o, dpi = 'retina', height = 7, width = 7)
ggsave("order_bxplot_t.png",plot=class_boxplot_t, dpi = 'retina', height = 7, width = 7)
ggsave("order_bxplot_l.png",plot=class_boxplot_l, dpi = 'retina', height = 7, width = 10)

top_box <- ggarrange(class_boxplot_o,class_boxplot_t,ncol = 2, labels = "AUTO")
box_fig <- ggarrange(top_box,class_boxplot_l, ncol = 1,nrow = 2, labels =c("", "C"))
ggsave("boxplot_fig_fams.png",plot=box_fig, dpi = 'retina', height = 10, width = 10)

