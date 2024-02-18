library(tidyverse); library(DESeq2); library(ggplot2); library(edgeR); library(dynamicTreeCut)
library(circlize); library(dendextend); library(ggpubr); library(VennDiagram)
################################## Set Up ##################################

setwd("/Users/lindseyschwartz/Documents/SIFP/R Files/MGAL_new-genome/MolEcol2")
SampleInfo <- read.csv("SampleTable.csv")[,-2]

##read in the counts file and split protein from lncRNA
countsMat <- read.delim("MGAL_counts.txt")
counts_lncrna <- countsMat[grepl("NC", countsMat$Geneid),]
counts_NC <-  counts_lncrna[,-c(1:6)]
rownames(counts_NC) <- counts_lncrna$Geneid
keepNC <- rowSums(counts_lncrna >= 5) >= 24
NCslim <- counts_NC[keepNC,]

##read in ancestry coefficients
meanQ <- read.delim("miss30_pF_out.2.meanQ.txt")
meanQ$spp <- ifelse(meanQ$spp == 1, "GG", ifelse(meanQ$spp == 2, "EG", "EE"))
colnames(meanQ) <- c("genotype", "sample", "MG", "ME")
SampleInfo <- left_join(meanQ, SampleInfo, by = c("sample", "genotype"))

SampleInfo$genotype <- ifelse(SampleInfo$MG >= 0.95, "gallo.like", ifelse(SampleInfo$MG <= 0.05, "edulis.like", "hybrid"))
SampleInfo <- SampleInfo[order(SampleInfo$MG, decreasing = T),]

##create design matrix
coldat <- SampleInfo %>% remove_rownames() %>% column_to_rownames(var = "sample")
coldat <- coldat[order(rownames(coldat)),]

coldat$temp <- as.factor(coldat$temp) 
coldat$temp <- relevel(coldat$temp, ref = "C")
coldat$genotype <- as.factor(coldat$genotype)
coldat$genotype <- relevel(coldat$genotype, ref = "edulis.like")

############## Run DESeq ##############
#des.NC <- DESeqDataSetFromMatrix(countData = NCslim, colData = coldat, design = ~ genotype + temp + genotype:temp)
#des.NC <- DESeq(des.NC)
#saveRDS(des.NC, "desProtSlim.RDS")
des.NC <- readRDS("desProtSlim.RDS")
resultsNames(des.NC)
vst.NC <- vst(des.NC)
vst.NC.counts <- assay(vst.NC)

############## Protein Main Effects ##############

##main effect of temp in edulis
NC.tempEE <- lfcShrink(des.NC, coef = "temp_H_vs_C", type = "ashr")
tempEE.ordered <- NC.tempEE[order(NC.tempEE$padj),]
NCEE.genes <- which(tempEE.ordered$padj<0.05)
NCEE.sig <- as.data.frame(tempEE.ordered[NCEE.genes,])
NCEE.sig.up <- NCEE.sig[which(NCEE.sig$log2FoldChange > 0), ]
NCEE.sig.down <- NCEE.sig[which(NCEE.sig$log2FoldChange < 0), ]
NCEE.sig.genes <- rownames(NCEE.sig)
length(NCEE.sig.genes)

##main effect of temp in gallo
NC.tempGG <- lfcShrink(des.NC, contrast = list(c("temp_H_vs_C", "genotypegallo.like.tempH")), type = "ashr")
tempGG.ordered <- NC.tempGG[order(NC.tempGG$padj),]
NCGG.genes <- which(tempGG.ordered$padj<0.05)
NCGG.sig <- as.data.frame(tempGG.ordered[NCGG.genes,])
NCGG.sig.up <- NCGG.sig[which(NCGG.sig$log2FoldChange > 0), ]
NCGG.sig.down <- NCGG.sig[which(NCGG.sig$log2FoldChange < 0), ]
NCGG.sig.genes <- rownames(NCGG.sig)
length(NCGG.sig.genes)

##main effect of temp in hybrids
NC.tempEG <- lfcShrink(des.NC, contrast = list(c("temp_H_vs_C", "genotypehybrid.tempH")), type = "ashr")
tempEG.ordered <- NC.tempEG[order(NC.tempEG$padj),]
NCEG.genes <- which(tempEG.ordered$padj<0.05)
NCEG.sig <- as.data.frame(tempEG.ordered[NCEG.genes,])
NCEG.sig.up <- NCEG.sig[which(NCEG.sig$log2FoldChange > 0), ]
NCEG.sig.down <- NCEG.sig[which(NCEG.sig$log2FoldChange < 0), ]
NCEG.sig.genes <- rownames(NCEG.sig)
length(NCEG.sig.genes)

############## Protein Interactions ##############

##interaction terms for EE v GG
NC.intGG <- lfcShrink(des.NC, coef = "genotypegallo.like.tempH", type = "ashr")
NC.intGG.ordered <- NC.intGG[order(NC.intGG$padj),]
NC.intGG.genes <- which(NC.intGG$padj<0.05)
NC.intGG.sig <- as.data.frame(NC.intGG[NC.intGG.genes,])
NC.intGG.sig <- NC.intGG.sig[order(NC.intGG.sig$padj),]
NC.intGG.sig.genes <- rownames(NC.intGG.sig)
length(NC.intGG.sig.genes)

##interaction terms for EE v EG
NC.intEG <- lfcShrink(des.NC, coef = "genotypehybrid.tempH", type = "ashr")
NC.intEG.ordered <- NC.intEG[order(NC.intEG$padj),]
NC.intEG.genes <- which(NC.intEG$padj<0.05)
NC.intEG.sig <- as.data.frame(NC.intEG[NC.intEG.genes,])
NC.intEG.sig <- NC.intEG.sig[order(NC.intEG.sig$padj),]
NC.intEG.sig.genes <- rownames(NC.intEG.sig)
length(NC.intEG.sig.genes)

##interaction terms for EG v GG
NC.int.GGEG <- lfcShrink(des.NC, contrast = list(c("genotypehybrid.tempH", "genotypegallo.like.tempH")), type = "ashr")
NC.int.GGEG.ordered <- NC.int.GGEG[order(NC.int.GGEG$padj),]
NC.int.GGEG.genes <- which(NC.int.GGEG$padj<0.05)
NC.int.GGEG.sig <- as.data.frame(NC.int.GGEG[NC.int.GGEG.genes,])
NC.int.GGEG.sig <- NC.int.GGEG.sig[order(NC.int.GGEG.sig$padj),]
NC.int.GGEG.sig.genes <- rownames(NC.int.GGEG.sig)
length(NC.int.GGEG.sig.genes)


############## Plots ##############

edulis.labs <- c(expression(paste(italic("M. edulis"), " 15°")), expression(paste(italic("M. edulis"), " 23°")))
hybrid.labs <- c("Hybrids 15°", "Hybrids 23°")
gallo.labs <- c(expression(paste(italic("M. galloprovincialis"), " 15°")), expression(paste(italic("M. galloprovincialis"), " 23°")))
labs.plot <- c(edulis.labs, gallo.labs, hybrid.labs)

cols6 <- c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02", "#7570b3", "#7570b3")
cols3 <- c("#1b9e77", "#d95f02", "#7570b3")
cols2 <- c("#d95f02", "#1b9e77")

genotemp <- paste(coldat$genotype, coldat$temp, sep = "_")

##Venn diagrams

NCEE.up = which(NCEE.sig$log2FoldChange > 0); NCEE.down = which(NCEE.sig$log2FoldChange < 0)
NCGG.up = which(NCGG.sig$log2FoldChange > 0); NCGG.down = which(NCGG.sig$log2FoldChange < 0)
NCEG.up = which(NCEG.sig$log2FoldChange > 0); NCEG.down = which(NCEG.sig$log2FoldChange < 0)

NC.up.venn = venn.diagram(x = list(rownames(NCEE.sig[NCEE.up,]),
                                   rownames(NCGG.sig[NCGG.up,]), 
                                   rownames(NCEG.sig[NCEG.up,])), 
                          category.names = c(" ", " ", " "),filename = NULL,
                          col = cols3, disable.logging = T, 
                          fill = cols3, cex = 2.5, output = F)

NC.down.venn = venn.diagram(x = list(rownames(NCEE.sig[NCEE.down,]),
                                     rownames(NCGG.sig[NCGG.down,]), 
                                     rownames(NCEG.sig[NCEG.down,])), 
                            category.names = c(" ", " ", " "),filename = NULL,
                            col = cols3, disable.logging = T, 
                            fill = cols3, cex = 2.5, output = F)

venns = ggarrange(NC.up.venn, NC.down.venn, ncol = 2, 
                  labels = c("A", "B"), font.label = list(size = 18, color = "black"))

venns

##PCA

DESeq2::plotPCA(vst.NC, intgroup = c("genotype", "temp"), ntop = 500)
pca.bin <- DESeq2::plotPCA(vst.NC, intgroup = c("genotype", "temp"), returnData = T, ntop = 500)

bin.pca.plot <- ggplot(data = pca.bin, aes(PC1, PC2, color = genotemp, shape = genotemp)) + 
  geom_point(size = 2.5) + theme_classic(base_size = 18) + xlab("PC1 (13%)") + ylab("PC2 (4%)") +
  scale_color_manual(values = brew6, name = "Group", labels = labs.plot) + 
  scale_shape_manual(values = c(19, 17, 19, 17, 19, 17), name = "Group", labels = labs.plot) + 
  theme(legend.position = "bottom")

bin.pca.plot

##Interaction term clustering

int_all_NC <- union(NC.int.GGEG.sig.genes, union(NC.intGG.sig.genes, NC.intEG.sig.genes))
length(int_all_NC)
###excluding this, as there are not enough for a meaningful analysis

##LFC Plot

edulis.LFC <- as.data.frame(cbind(rownames(NC.tempEE), NC.tempEE$log2FoldChange))
gallo.LFC <- as.data.frame(cbind(rownames(NC.tempGG), NC.tempGG$log2FoldChange))
hybrid.LFC <- as.data.frame(cbind(rownames(NC.tempEG), NC.tempEG$log2FoldChange))

lfc.table <- cbind(edulis.LFC, gallo.LFC, hybrid.LFC)
lfc.table <- lfc.table[,c(1,2,4,6)]; colnames(lfc.table) = c("gene_id", "edulis", "gallo", "hybrid")
lfc.table$edulis <- as.numeric(lfc.table$edulis)
lfc.table$gallo <- as.numeric(lfc.table$gallo)
lfc.table$hybrid <- as.numeric(lfc.table$hybrid)

lfc.table$eg.ee <- lfc.table$hybrid - lfc.table$edulis
lfc.table$eg.gg <- lfc.table$hybrid - lfc.table$gallo

y.lab.italic <- expression(paste("Log"[2], " FC(hybrids) - Log"[2], " FC(", italic("M. galloprovincialis"), ")"))
x.lab.italic <- expression(paste("Log"[2], " FC(hybrids) - Log"[2], " FC(", italic("M. edulis"), ")"))


deUnion <- union(int_all_NC, union(NCEG.sig.genes, union(NCEE.sig.genes, NCGG.sig.genes)))
deUnion.indices <- which(lfc.table$gene_id %in% deUnion == T)
lfc.table.union <- lfc.table[deUnion.indices,]

me.dom <- which(abs(lfc.table.union$eg.ee)>0.32 & abs(lfc.table.union$eg.gg)<0.32)
mg.dom <- which(abs(lfc.table.union$eg.gg)>0.32 & abs(lfc.table.union$eg.ee)<0.32)
overdom <- which(lfc.table.union$eg.gg>0.32 & lfc.table.union$eg.ee>0.32)
underdom <- which(lfc.table.union$eg.gg<(-0.32) & lfc.table.union$eg.ee<(-0.32))

ggplot() + geom_point(data = lfc.table, aes(eg.ee, eg.gg), color = "gray") + 
  theme_classic(base_size = 15) + xlab(x.lab.italic) + ylab(y.lab.italic) + 
  geom_point(data = lfc.table.union[me.dom,], aes(eg.ee, eg.gg), color = "#1b9e77") + 
  geom_point(data = lfc.table.union[mg.dom,], aes(eg.ee, eg.gg), color = "#d95f02") + 
  geom_point(data = lfc.table.union[overdom,], aes(eg.ee, eg.gg), color = "firebrick2") +
  geom_point(data = lfc.table.union[underdom,], aes(eg.ee, eg.gg), color = "firebrick2") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  xlim(c(-3,3)) + ylim(c(-3,3))


