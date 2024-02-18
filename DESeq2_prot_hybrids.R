library(tidyverse); library(DESeq2); library(ggplot2); library(edgeR); library(dynamicTreeCut)
library(circlize); library(dendextend); library(ggpubr); library(VennDiagram)
################################## Set Up ##################################

setwd("/Users/lindseyschwartz/Documents/SIFP/R Files/MGAL_new-genome/MolEcol2")
SampleInfo <- read.csv("SampleTable.csv")[,-2]

##read in the counts file and split protein from lncRNA
countsMat <- read.delim("MGAL_counts.txt")
counts.prot <- countsMat[grepl("MGAL_10B", countsMat$Geneid),]
keep <- rowSums(counts.prot >= 5) >= 24
prot.slim <- counts.prot[keep,] %>% remove_rownames() %>%  column_to_rownames(var = "Geneid")
prot.slim.meta <- prot.slim[,1:5]
prot.slim <- prot.slim[,-c(1:5)]

##read in ancestry coefficients
meanQ <- read.delim("miss30_pF_out.2.meanQ.txt")
meanQ$spp <- ifelse(meanQ$spp == 1, "GG", ifelse(meanQ$spp == 2, "EG", "EE"))
colnames(meanQ) <- c("genotype", "sample", "MG", "ME")
SampleInfo <- left_join(meanQ, SampleInfo, by = c("sample", "genotype"))
SampleInfo <- SampleInfo[order(SampleInfo$MG, decreasing = T),]
hybrid <- SampleInfo[which(SampleInfo$genotype=="EG"),]

##create design matrix
coldat <- hybrid %>% remove_rownames() %>% column_to_rownames(var = "sample")
coldat$subtype <- c(rep("gallo", 11), rep("hybrid", 11), rep("edulis", 10))
coldat$temp <- as.factor(coldat$temp) 
coldat$temp <- relevel(coldat$temp, ref = "C")
coldat$subtype <- as.factor(coldat$subtype)
coldat$subtype <- relevel(coldat$subtype, ref = "edulis")

##ammend counts matrix
prot.hybrid <- prot.slim[,which(colnames(prot.slim) %in% rownames(coldat))]
prot.hybrid <- prot.hybrid[,order(match(colnames(prot.hybrid), rownames(coldat)))]

############## Run DESeq ##############
#des.prot <- DESeqDataSetFromMatrix(countData = prot.hybrid, colData = coldat, design = ~ subtype + temp + subtype:temp)
#des.prot <- DESeq(des.prot)
#saveRDS(des.prot, "desProtHybrid.RDS")
des.prot <- readRDS("desProtHybrid.RDS")
resultsNames(des.prot)
vst.prot <- vst(des.prot)
vst.prot.counts <- assay(vst.prot)

############## Protein Main Effects ##############

##main effect of temp in edulis
prot.tempEE <- lfcShrink(des.prot, coef = "temp_H_vs_C", type = "ashr")
tempEE.ordered <- prot.tempEE[order(prot.tempEE$padj),]
protEE.genes <- which(tempEE.ordered$padj<0.05)
protEE.sig <- as.data.frame(tempEE.ordered[protEE.genes,])
protEE.sig.up <- protEE.sig[which(protEE.sig$log2FoldChange > 0), ]
protEE.sig.down <- protEE.sig[which(protEE.sig$log2FoldChange < 0), ]
protEE.sig.genes <- rownames(protEE.sig)
length(protEE.sig.genes)

##main effect of temp in gallo
prot.tempGG <- lfcShrink(des.prot, contrast = list(c("temp_H_vs_C", "subtypegallo.tempH")), type = "ashr")
tempGG.ordered <- prot.tempGG[order(prot.tempGG$padj),]
protGG.genes <- which(tempGG.ordered$padj<0.05)
protGG.sig <- as.data.frame(tempGG.ordered[protGG.genes,])
protGG.sig.up <- protGG.sig[which(protGG.sig$log2FoldChange > 0), ]
protGG.sig.down <- protGG.sig[which(protGG.sig$log2FoldChange < 0), ]
protGG.sig.genes <- rownames(protGG.sig)
length(protGG.sig.genes)

##main effect of temp in hybrids
prot.tempEG <- lfcShrink(des.prot, contrast = list(c("temp_H_vs_C", "subtypehybrid.tempH")), type = "ashr")
tempEG.ordered <- prot.tempEG[order(prot.tempEG$padj),]
protEG.genes <- which(tempEG.ordered$padj<0.05)
protEG.sig <- as.data.frame(tempEG.ordered[protEG.genes,])
protEG.sig.up <- protEG.sig[which(protEG.sig$log2FoldChange > 0), ]
protEG.sig.down <- protEG.sig[which(protEG.sig$log2FoldChange < 0), ]
protEG.sig.genes <- rownames(protEG.sig)
length(protEG.sig.genes)

############## Protein Interactions ##############

##interaction terms for EE v GG
prot.intGG <- lfcShrink(des.prot, coef = "subtypegallo.tempH", type = "ashr")
prot.intGG.ordered <- prot.intGG[order(prot.intGG$padj),]
prot.intGG.genes <- which(prot.intGG$padj<0.05)
prot.intGG.sig <- as.data.frame(prot.intGG[prot.intGG.genes,])
prot.intGG.sig <- prot.intGG.sig[order(prot.intGG.sig$padj),]
prot.intGG.sig.genes <- rownames(prot.intGG.sig)
length(prot.intGG.sig.genes)

##interaction terms for EE v EG
prot.intEG <- lfcShrink(des.prot, coef = "subtypehybrid.tempH", type = "ashr")
prot.intEG.ordered <- prot.intEG[order(prot.intEG$padj),]
prot.intEG.genes <- which(prot.intEG$padj<0.05)
prot.intEG.sig <- as.data.frame(prot.intEG[prot.intEG.genes,])
prot.intEG.sig <- prot.intEG.sig[order(prot.intEG.sig$padj),]
prot.intEG.sig.genes <- rownames(prot.intEG.sig)
length(prot.intEG.sig.genes)

##interaction terms for EG v GG
prot.int.GGEG <- lfcShrink(des.prot, contrast = list(c("subtypehybrid.tempH", "subtypegallo.tempH")), type = "ashr")
prot.int.GGEG.ordered <- prot.int.GGEG[order(prot.int.GGEG$padj),]
prot.int.GGEG.genes <- which(prot.int.GGEG$padj<0.05)
prot.int.GGEG.sig <- as.data.frame(prot.int.GGEG[prot.int.GGEG.genes,])
prot.int.GGEG.sig <- prot.int.GGEG.sig[order(prot.int.GGEG.sig$padj),]
prot.int.GGEG.sig.genes <- rownames(prot.int.GGEG.sig)
length(prot.int.GGEG.sig.genes)


############## Plots ##############
edulis.labs <- c(expression(paste(italic("M. edulis"), " 15°")), expression(paste(italic("M. edulis"), " 23°")))
hybrid.labs <- c("Hybrids 15°", "Hybrids 23°")
gallo.labs <- c(expression(paste(italic("M. galloprovincialis"), " 15°")), expression(paste(italic("M. galloprovincialis"), " 23°")))
labs.plot <- c(edulis.labs, gallo.labs, hybrid.labs)

brew6 <- c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02", "#7570b3", "#7570b3")
genotemp <- paste(coldat$subtype, coldat$temp, sep = "_")

##PCA

DESeq2::plotPCA(vst.prot, intgroup = c("subtype", "temp"), ntop = 500)
pca.bin <- DESeq2::plotPCA(vst.prot, intgroup = c("subtype", "temp"), returnData = T, ntop = 500)

bin.pca.plot <- ggplot(data = pca.bin, aes(PC1, PC2, color = genotemp, shape = genotemp)) + 
  geom_point(size = 2.5) + theme_classic(base_size = 18) + xlab("PC1 (9%)") + ylab("PC2 67%)") +
  scale_color_manual(values = brew6, name = "Group", labels = labs.plot) + 
  scale_shape_manual(values = c(19, 17, 19, 17, 19, 17), name = "Group", labels = labs.plot) + 
  theme(legend.position = "bottom")

bin.pca.plot

ancestry <- coldat %>% rownames_to_column("name") %>% select(name, MG, temp, subtype)
pca.continuous <- left_join(pca.bin, ancestry, by = c("name", "temp", "subtype"))

pca.cont.plot <- ggplot(data = pca.continuous, aes(PC1, PC2, color = MG, shape = temp)) + 
  geom_point(size = 2.5) + theme_classic(base_size = 18) + xlab("PC1 (8%)") + ylab("PC2 (7%)") +
  scale_color_viridis_c() + theme(legend.position = "bottom") + 
  scale_shape_manual(values = c(19, 17, 19, 17, 19, 17), name = "Group")
  
pca.cont.plot

ggarrange(bin.pca.plot, pca.cont.plot, nrow = 1)

##LFC Plot
int_all_prot <- union(prot.int.GGEG.sig.genes, union(prot.intGG.sig.genes, prot.intEG.sig.genes))

edulis.LFC <- as.data.frame(cbind(rownames(prot.tempEE), prot.tempEE$log2FoldChange))
gallo.LFC <- as.data.frame(cbind(rownames(prot.tempGG), prot.tempGG$log2FoldChange))
hybrid.LFC <- as.data.frame(cbind(rownames(prot.tempEG), prot.tempEG$log2FoldChange))

lfc.table <- cbind(edulis.LFC, gallo.LFC, hybrid.LFC)
lfc.table <- lfc.table[,c(1,2,4,6)]; colnames(lfc.table) = c("gene_id", "edulis", "gallo", "hybrid")
lfc.table$edulis <- as.numeric(lfc.table$edulis)
lfc.table$gallo <- as.numeric(lfc.table$gallo)
lfc.table$hybrid <- as.numeric(lfc.table$hybrid)

lfc.table$eg.ee <- lfc.table$hybrid - lfc.table$edulis
lfc.table$eg.gg <- lfc.table$hybrid - lfc.table$gallo

y.lab.italic <- expression(paste("Log"[2], " FC(true hybrids) - Log"[2], " FC(", italic("M. galloprovincialis-like"), ")"))
x.lab.italic <- expression(paste("Log"[2], " FC(true hybrids) - Log"[2], " FC(", italic("M. edulis-like"), ")"))


deUnion <- union(int_all_prot, union(protEG.sig.genes, union(protEE.sig.genes, protGG.sig.genes)))
deUnion.indices <- which(lfc.table$gene_id %in% deUnion == T)
lfc.table.union <- lfc.table[deUnion.indices,]

me.dom <- which(abs(lfc.table.union$eg.ee)>0.32 & abs(lfc.table.union$eg.gg)<0.32)
mg.dom <- which(abs(lfc.table.union$eg.gg)>0.32 & abs(lfc.table.union$eg.ee)<0.32)
overdom <- which(lfc.table.union$eg.gg>0.32 & lfc.table.union$eg.ee>0.32)
underdom <- which(lfc.table.union$eg.gg<(-0.32) & lfc.table.union$eg.ee<(-0.32))

ggplot() + geom_point(data = lfc.table, aes(eg.ee, eg.gg), color = "gray") + 
  theme_classic(base_size = 15) + xlab(x.lab.italic) + ylab(y.lab.italic) + 
  geom_point(data = lfc.table.union[me.dom,], aes(eg.ee, eg.gg), color = "#d95f02") + 
  geom_point(data = lfc.table.union[mg.dom,], aes(eg.ee, eg.gg), color = "#1b9e77") + 
  geom_point(data = lfc.table.union[overdom,], aes(eg.ee, eg.gg), color = "#0f4c5c") +
  geom_point(data = lfc.table.union[underdom,], aes(eg.ee, eg.gg), color = "#0f4c5c") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  xlim(c(-10,10)) + ylim(c(-10,10))








