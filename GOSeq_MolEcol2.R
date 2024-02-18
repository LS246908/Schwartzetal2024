################################## Set Up ##################################

library(goseq)
setwd("/Users/lindseyschwartz/Documents/SIFP/R files/MGAL_new-genome/GOSeq_MolEcol2/")

prot.genes <- rownames(prot.slim.meta)
length.prot <- prot.slim.meta$Length
names(length.prot) <- rownames(prot.slim.meta)
length.prot <- length.prot[order(match(length.prot, prot.genes))]
annot <- read.delim("goTerms_reformat_MGAL10B.txt")
hits.list <- str_split(annot$GO.IDs, ";")
names(hits.list) <- annot$SeqName

############################### M. edulis 15 v 23 ###############################
EE_CvH.up = rownames(protEE.sig.up)
EE_CvH.down = rownames(protEE.sig.down)

## Upregulated
edulis.up = ifelse(prot.genes %in% EE_CvH.up, 1, 0)
names(edulis.up) = prot.genes
pwf.ee.up = nullp(DEgenes = edulis.up, bias.data = length.prot)
go.ee.up = goseq(pwf = pwf.ee.up, gene2cat = hits.list, method = "Wallenius")
go.ee.up.over = which(go.ee.up$over_represented_pvalue < 0.01)
ee.up.overRep = go.ee.up[go.ee.up.over,]
write.csv(ee.up.overRep[,c(1,2,6,7)], "eeCvH_up_01.csv", row.names = F)

## Downregulated
edulis.down = ifelse(prot.genes %in% EE_CvH.down, 1, 0)
names(edulis.down) = prot.genes
pwf.ee.down = nullp(DEgenes = edulis.down, bias.data = length.prot)
go.ee.down = goseq(pwf = pwf.ee.down, gene2cat = hits.list, method = "Wallenius")
go.ee.down.over = which(go.ee.down$over_represented_pvalue < 0.01)
ee.down.overRep = go.ee.down[go.ee.down.over,]
write.csv(ee.down.overRep[,c(1,2,6,7)], "eeCvH_down_01.csv", row.names = F)

############################### M. gallo 15 v 23 ###############################
GG_CvH.up = rownames(protGG.sig.up)
GG_CvH.down = rownames(protGG.sig.down)

## Upregulated
gallo.up = ifelse(prot.genes %in% GG_CvH.up, 1, 0)
names(gallo.up) = prot.genes
pwf.gg.up = nullp(DEgenes = gallo.up, bias.data = length.prot)
go.gg.up = goseq(pwf = pwf.gg.up, gene2cat = hits.list, method = "Wallenius")
go.gg.up.over = which(go.gg.up$over_represented_pvalue < 0.01)
gg.up.overRep = go.gg.up[go.gg.up.over,]
write.csv(gg.up.overRep[,c(1,2,6,7)], "ggCvH_up_01.csv", row.names = F)

## Downregulated
gallo.down = ifelse(prot.genes %in% GG_CvH.down, 1, 0)
names(gallo.down) = prot.genes
pwf.gg.down = nullp(DEgenes = gallo.down, bias.data = length.prot)
go.gg.down = goseq(pwf = pwf.gg.down, gene2cat = hits.list, method = "Wallenius")
go.gg.down.over = which(go.gg.down$over_represented_pvalue < 0.01)
gg.down.overRep = go.gg.down[go.gg.down.over,]
write.csv(gg.down.overRep[,c(1,2,6,7)], "ggCvH_down_01.csv", row.names = F)

############################### Hybrids 15 v 23 ###############################
EG_CvH.up = rownames(protEG.sig.up)
EG_CvH.down = rownames(protEG.sig.down)

## Upregulated
hybrid.up = ifelse(prot.genes %in% EG_CvH.up, 1, 0)
names(hybrid.up) = prot.genes
pwf.eg.up = nullp(DEgenes = hybrid.up, bias.data = length.prot)
go.eg.up = goseq(pwf = pwf.eg.up, gene2cat = hits.list, method = "Wallenius")
go.eg.up.over = which(go.eg.up$over_represented_pvalue < 0.01)
eg.up.overRep = go.eg.up[go.eg.up.over,]
write.csv(eg.up.overRep[,c(1,2,6,7)], "egCvH_up_01.csv", row.names = F)

## Downregulated
hybrid.down = ifelse(prot.genes %in% EG_CvH.down, 1, 0)
names(hybrid.down) = prot.genes
pwf.eg.down = nullp(DEgenes = hybrid.down, bias.data = length.prot)
go.eg.down = goseq(pwf = pwf.eg.down, gene2cat = hits.list, method = "Wallenius")
go.eg.down.over = which(go.eg.down$over_represented_pvalue < 0.01)
eg.down.overRep = go.eg.down[go.eg.down.over,]
write.csv(eg.down.overRep[,c(1,2,6,7)], "egCvH_down_01.csv", row.names = F)

############################### Interaction Term ###############################

cluster1 = hclust.clusterIDs.prot[which(hclust.clusterIDs.prot$module==1),1]
cluster2 = hclust.clusterIDs.prot[which(hclust.clusterIDs.prot$module==2),1]
cluster3 = hclust.clusterIDs.prot[which(hclust.clusterIDs.prot$module==3),1]
cluster4 = hclust.clusterIDs.prot[which(hclust.clusterIDs.prot$module==4),1]
cluster5 = hclust.clusterIDs.prot[which(hclust.clusterIDs.prot$module==5),1]

#### Cluster 1 ####

c1.genes = ifelse(prot.genes %in% cluster1, 1, 0)
names(c1.genes) = prot.genes
pwf.c1 = nullp(DEgenes = c1.genes, bias.data = length.prot)
go.c1 = goseq(pwf = pwf.c1, gene2cat = hits.list, method = "Wallenius")
go.c1.over = which(go.c1$over_represented_pvalue < 0.01)
c1.over = go.c1[go.c1.over,]
write.table(c1.over[,c(1,2,6,7)], "cluster1_01.csv", sep = ",", row.names = F)

#### Cluster 2 ####

c2.genes = ifelse(prot.genes %in% cluster2, 1, 0)
names(c2.genes) = prot.genes
pwf.c2 = nullp(DEgenes = c2.genes, bias.data = length.prot)
go.c2 = goseq(pwf = pwf.c2, gene2cat = hits.list, method = "Wallenius")
go.c2.over = which(go.c2$over_represented_pvalue < 0.01)
c2.over = go.c2[go.c2.over,]
write.table(c2.over[,c(1,2,6,7)], "cluster2_01.csv", sep = ",", row.names = F)

#### Cluster 3 ####

c3.genes = ifelse(prot.genes %in% cluster3, 1, 0)
names(c3.genes) = prot.genes
pwf.c3 = nullp(DEgenes = c3.genes, bias.data = length.prot)
go.c3 = goseq(pwf = pwf.c3, gene2cat = hits.list, method = "Wallenius")
go.c3.over = which(go.c3$over_represented_pvalue < 0.01)
c3.over = go.c3[go.c3.over,]
write.table(c3.over[,c(1,2,6,7)], "cluster3_01.csv", sep = ",", row.names = F)

#### Cluster 4 ####

c4.genes = ifelse(prot.genes %in% cluster4, 1, 0)
names(c4.genes) = prot.genes
pwf.c4 = nullp(DEgenes = c4.genes, bias.data = length.prot)
go.c4 = goseq(pwf = pwf.c4, gene2cat = hits.list, method = "Wallenius")
go.c4.over = which(go.c4$over_represented_pvalue < 0.01)
c4.over = go.c4[go.c4.over,]
write.table(c4.over[,c(1,2,6,7)], "cluster4_01.csv", sep = ",", row.names = F)

#### Cluster 5 ####

c5.genes = ifelse(prot.genes %in% cluster5, 1, 0)
names(c5.genes) = prot.genes
pwf.c5 = nullp(DEgenes = c5.genes, bias.data = length.prot)
go.c5 = goseq(pwf = pwf.c5, gene2cat = hits.list, method = "Wallenius")
go.c5.over = which(go.c5$over_represented_pvalue < 0.01)
c5.over = go.c5[go.c5.over,]
write.table(c5.over[,c(1,2,6,7)], "cluster5_01.csv", sep = ",", row.names = F)

