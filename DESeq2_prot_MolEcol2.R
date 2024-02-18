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
#des.prot <- DESeqDataSetFromMatrix(countData = prot.slim, colData = coldat, design = ~ genotype + temp + genotype:temp)
#des.prot <- DESeq(des.prot)
#saveRDS(des.prot, "desProtSlim.RDS")
des.prot <- readRDS("desProtSlim.RDS")
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
prot.tempGG <- lfcShrink(des.prot, contrast = list(c("temp_H_vs_C", "genotypegallo.like.tempH")), type = "ashr")
tempGG.ordered <- prot.tempGG[order(prot.tempGG$padj),]
protGG.genes <- which(tempGG.ordered$padj<0.05)
protGG.sig <- as.data.frame(tempGG.ordered[protGG.genes,])
protGG.sig.up <- protGG.sig[which(protGG.sig$log2FoldChange > 0), ]
protGG.sig.down <- protGG.sig[which(protGG.sig$log2FoldChange < 0), ]
protGG.sig.genes <- rownames(protGG.sig)
length(protGG.sig.genes)

##main effect of temp in hybrids
prot.tempEG <- lfcShrink(des.prot, contrast = list(c("temp_H_vs_C", "genotypehybrid.tempH")), type = "ashr")
tempEG.ordered <- prot.tempEG[order(prot.tempEG$padj),]
protEG.genes <- which(tempEG.ordered$padj<0.05)
protEG.sig <- as.data.frame(tempEG.ordered[protEG.genes,])
protEG.sig.up <- protEG.sig[which(protEG.sig$log2FoldChange > 0), ]
protEG.sig.down <- protEG.sig[which(protEG.sig$log2FoldChange < 0), ]
protEG.sig.genes <- rownames(protEG.sig)
length(protEG.sig.genes)

############## Protein Interactions ##############

##interaction terms for EE v GG
prot.intGG <- lfcShrink(des.prot, coef = "genotypegallo.like.tempH", type = "ashr")
prot.intGG.ordered <- prot.intGG[order(prot.intGG$padj),]
prot.intGG.genes <- which(prot.intGG$padj<0.05)
prot.intGG.sig <- as.data.frame(prot.intGG[prot.intGG.genes,])
prot.intGG.sig <- prot.intGG.sig[order(prot.intGG.sig$padj),]
prot.intGG.sig.genes <- rownames(prot.intGG.sig)
length(prot.intGG.sig.genes)

##interaction terms for EE v EG
prot.intEG <- lfcShrink(des.prot, coef = "genotypehybrid.tempH", type = "ashr")
prot.intEG.ordered <- prot.intEG[order(prot.intEG$padj),]
prot.intEG.genes <- which(prot.intEG$padj<0.05)
prot.intEG.sig <- as.data.frame(prot.intEG[prot.intEG.genes,])
prot.intEG.sig <- prot.intEG.sig[order(prot.intEG.sig$padj),]
prot.intEG.sig.genes <- rownames(prot.intEG.sig)
length(prot.intEG.sig.genes)

##interaction terms for EG v GG
prot.int.GGEG <- lfcShrink(des.prot, contrast = list(c("genotypehybrid.tempH", "genotypegallo.like.tempH")), type = "ashr")
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

cols6 <- c("#1b9e77", "#1b9e77", "#d95f02", "#d95f02", "#7570b3", "#7570b3")
cols3 <- c("#1b9e77", "#d95f02", "#7570b3")
cols2 <- c("#d95f02", "#1b9e77")

genotemp <- paste(coldat$genotype, coldat$temp, sep = "_")

##Venn diagrams

protEE.up = which(protEE.sig$log2FoldChange > 0); protEE.down = which(protEE.sig$log2FoldChange < 0)
protGG.up = which(protGG.sig$log2FoldChange > 0); protGG.down = which(protGG.sig$log2FoldChange < 0)
protEG.up = which(protEG.sig$log2FoldChange > 0); protEG.down = which(protEG.sig$log2FoldChange < 0)

prot.up.venn = venn.diagram(x = list(rownames(protEE.sig[protEE.up,]),
                                     rownames(protGG.sig[protGG.up,]), 
                                     rownames(protEG.sig[protEG.up,])), 
                            category.names = c(" ", " ", " "),filename = NULL,
                            col = cols3, disable.logging = T, 
                            fill = cols3, cex = 2.5, output = F)

prot.down.venn = venn.diagram(x = list(rownames(protEE.sig[protEE.down,]),
                                       rownames(protGG.sig[protGG.down,]), 
                                       rownames(protEG.sig[protEG.down,])), 
                              category.names = c(" ", " ", " "),filename = NULL,
                              col = cols3, disable.logging = T, 
                              fill = cols3, cex = 2.5, output = F)

venns = ggarrange(prot.up.venn, prot.down.venn, ncol = 2, 
                  labels = c("A", "B"), font.label = list(size = 18, color = "black"))

venns

##PCA

DESeq2::plotPCA(vst.prot, intgroup = c("genotype", "temp"), ntop = nrow(counts.prot))
pcaProt <- DESeq2::plotPCA(vst.prot, intgroup = c("genotype", "temp"), 
                           returnData = T, ntop = nrow(counts.prot))

pc1 <- round((attr(pcaProt, "percentVar")[1])*100, 2)
pc2 <- round((attr(pcaProt, "percentVar")[2])*100, 2)
pcaPlot <- ggplot(data = pcaProt, aes(PC1, PC2, color = genotemp, shape = genotemp)) + 
  geom_point(size = 2.5) + theme_classic(base_size = 18) + 
  xlab(paste0("PC1 (", pc1, "%)")) + ylab(paste0("PC2 (", pc2, "%)")) +
  scale_color_manual(values = cols6, name = "Group", labels = labs.plot) + 
  scale_shape_manual(values = c(19, 17, 19, 17, 19, 17), name = "Group", labels = labs.plot) + 
  theme(legend.position = "bottom")

pcaPlot

#Top 500 only

pcaTop500 <- DESeq2::plotPCA(vst.prot, intgroup = c("genotype", "temp"), returnData = T, ntop = 500)

pc1t500 <- round((attr(pcaTop500, "percentVar")[1])*100, 2)
pc2t500 <- round((attr(pcaTop500, "percentVar")[2])*100, 2)
pcaPlot500 <- ggplot(data = pcaTop500, aes(PC1, PC2, color = genotemp, shape = genotemp)) + 
  geom_point(size = 2.5) + theme_classic(base_size = 18) + 
  xlab(paste0("PC1 (", pc1t500, "%)")) + ylab(paste0("PC2 (", pc2t500, "%)")) +
  scale_color_manual(values = cols6, name = "Group", labels = labs.plot) + 
  scale_shape_manual(values = c(19, 17, 19, 17, 19, 17), name = "Group", labels = labs.plot) + 
  theme(legend.position = "bottom")

pcaPlot500

ggarrange(pcaPlot, pcaPlot500, ncol = 2, common.legend = T, legend = "bottom")

##Interaction term clustering

int_all_prot <- union(prot.int.GGEG.sig.genes, union(prot.intGG.sig.genes, prot.intEG.sig.genes))
coldat_int <- cbind(coldat[,c(1,4)], genotemp)
colnames(coldat_int) <- c("genotype", "temp", "gentemp")
coldat_int <- as.data.frame(coldat_int %>% rownames_to_column("sample"))

intMat.prot <- t(scale(t(vst.prot.counts[int_all_prot,])))
clust.interact <- hclust(as.dist(1-cor(t(intMat.prot), method="spearman")), method="complete")

int_all_tree <- as.dendrogram(clust.interact, method="complete")
plot(int_all_tree, leaflab = "none", main = "Gene Clustering", ylab = "Height")
k4 <- cutree(int_all_tree, k = 4)
cols4 <- c("#E69F00", "#009E73", "#0072B2", "#000000", "#CC79A7", "#D55D00")[k4]
k5 <- cutree(int_all_tree, k = 5)
cols5 <- c("#E69F00", "#009E73", "#0072B2", "#000000", "#CC79A7", "#D55D00")[k5]
k6 <- cutree(int_all_tree, k = 6)
cols6 <- c("#E69F00", "#009E73", "#0072B2", "#000000",  "#D55D00", "#CC79A7")[k6]
colored_bars(cols4, int_all_tree, sort_by_labels_order = T, 
             y_shift=-0.1, rowLabels = c("k = 4"),cex.rowLabels=0.7)
colored_bars(cols5, int_all_tree, sort_by_labels_order = T, 
             y_shift=-0.2, rowLabels = c("k = 5"),cex.rowLabels=0.7)
colored_bars(cols6, int_all_tree, sort_by_labels_order = T, 
             y_shift=-0.3, rowLabels = c("k = 6"),cex.rowLabels=0.7)

hclust.clusterIDs.prot <- as.data.frame(cutree(int_all_tree, k = 5)) %>% rownames_to_column()
colnames(hclust.clusterIDs.prot) <- c("gene_id", "module")

clusters.prot <- intMat.prot %>% as.data.frame() %>% rownames_to_column("gene_id") %>% gather(-gene_id, key = "sample", value = "expression") %>% 
  left_join(coldat_int) %>%  left_join(hclust.clusterIDs.prot, by = "gene_id", copy = T) %>% group_by(genotype, temp, module) %>% 
  summarize(n_samples = n(), median_expression = median(expression), mean_expression = mean(expression), se_expression = sd(expression)/sqrt(n_samples)) %>% 
  ggplot(aes(x = temp, y = mean_expression, group = genotype, color = genotype)) + geom_point() + 
  geom_errorbar(aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression), width = 0.15) +
  facet_grid(~module, scales = "free") + geom_line() + theme_classic(base_size = 12) + 
  scale_color_manual(values = cols3, labels = c(expression(italic("M. edulis")), "Hybrids", expression(italic("M. galloprovincialis")))) +
  theme(text = element_text(size=22), legend.text.align = 0) + labs(x = "Temperature", y = "Mean Expression", color = "Genotype") +
  theme(strip.text.y = element_text(angle = 0)) + theme(legend.position = "bottom")

clusters.prot

##LFC Plot

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

y.lab.italic <- expression(paste("Log"[2], " FC(hybrids) - Log"[2], " FC(", italic("M. galloprovincialis"), ")"))
x.lab.italic <- expression(paste("Log"[2], " FC(hybrids) - Log"[2], " FC(", italic("M. edulis"), ")"))


deUnion <- union(int_all_prot, union(protEG.sig.genes, union(protEE.sig.genes, protGG.sig.genes)))
deUnion.indices <- which(lfc.table$gene_id %in% deUnion == T)
lfc.table.union <- lfc.table[deUnion.indices,]

me.dom <- which(abs(lfc.table.union$eg.ee)>0.32 & abs(lfc.table.union$eg.gg)<0.32)
mg.dom <- which(abs(lfc.table.union$eg.gg)>0.32 & abs(lfc.table.union$eg.ee)<0.32)
overdom <- which(lfc.table.union$eg.gg>0.32 & lfc.table.union$eg.ee>0.32)
underdom <- which(lfc.table.union$eg.gg<(-0.32) & lfc.table.union$eg.ee<(-0.32))

length(which(lfc.table.union$eg.gg>0.32 & lfc.table.union$eg.ee<(-0.32)))
length(which(lfc.table.union$eg.gg<(-0.32) & lfc.table.union$eg.ee>0.32))


ggplot() + geom_point(data = lfc.table.union, aes(eg.ee, eg.gg), color = "gray") + 
  theme_classic(base_size = 15) + xlab(x.lab.italic) + ylab(y.lab.italic) + 
  geom_point(data = lfc.table.union[me.dom,], aes(eg.ee, eg.gg), color = "#d95f02") + 
  geom_point(data = lfc.table.union[mg.dom,], aes(eg.ee, eg.gg), color = "#1b9e77") + 
  geom_point(data = lfc.table.union[overdom,], aes(eg.ee, eg.gg), color = "#0f4c5c") +
  geom_point(data = lfc.table.union[underdom,], aes(eg.ee, eg.gg), color = "#0f4c5c") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  xlim(c(-4,4)) + ylim(c(-4,4))

ggplot() + geom_point(data = lfc.table.union, aes(eg.ee, eg.gg), color = "gray") + 
  theme_classic(base_size = 15) + xlab(x.lab.italic) + ylab(y.lab.italic) + 
  geom_point(data = lfc.table.union[me.dom,], aes(eg.ee, eg.gg), color = "#d95f02") + 
  geom_point(data = lfc.table.union[mg.dom,], aes(eg.ee, eg.gg), color = "#1b9e77") + 
  geom_point(data = lfc.table.union[overdom,], aes(eg.ee, eg.gg), color = "#0f4c5c") +
  geom_point(data = lfc.table.union[underdom,], aes(eg.ee, eg.gg), color = "#0f4c5c") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) 

### Population Structure Plots

fastS30 = t((SampleInfo[,3:4]))
barplot(fastS30, col = cols2, main = "70% missing data")

##by site
PQrows <- which(grepl("PQ", SampleInfo$sample))
TWrows <- which(grepl("TW", SampleInfo$sample))
WBrows <- which(grepl("WB", SampleInfo$sample))
SDrows <- which(grepl("SD", SampleInfo$sample))

par(mfrow=c(1,4))
barplot(fastS30[,PQrows], col = cols2)
barplot(fastS30[,TWrows], col = cols2)
barplot(fastS30[,WBrows], col = cols2)
barplot(fastS30[,SDrows], col = cols2)
dev.off()

##by gebotype group
MGrows <- which(SampleInfo$genotype=="gallo.like")
EGrows <- which(SampleInfo$genotype=="hybrid")
MErows <- which(SampleInfo$genotype=="edulis.like")

par(mfrow=c(1,3))
barplot(fastS30[,MGrows], col = cols2)
barplot(fastS30[,EGrows], col = cols2)
barplot(fastS30[,MErows], col = cols2)
dev.off()

########################### Mapping Data ########################### 
mapping <- read.csv("mapping_percentages.csv", header = F)
colnames(mapping) <- c("ID", "genotype", "map")
model <- aov(map~genotype, data = mapping)
summary(model)
TukeyHSD(model, conf.level=.95)






