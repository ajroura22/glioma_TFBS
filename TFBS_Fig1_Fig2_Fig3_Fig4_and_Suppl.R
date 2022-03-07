## -

## Purpose of script: TF prediction in ATAC-seq and its integration with transcriptomic datasets (cell line RNAseq, TCGA) and histone markers (Stepniak et al., Nat. Comms.2021)
##
## Author: Adria-Jaume Roura
##
## Date Created: 2020-22-10
##

## libraries used
library(preprocessCore)
library(biomaRt)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(DESeq2)
library(ReactomePA)
library(clusterProfiler)
library(GenomicFeatures)
library(stringr)
library(BiocParallel)
library(gplots)
library(ggpmisc)
library(tidygenomics)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RColorBrewer)
library(VennDiagram)
library(reshape2)
library(dplyr)
library(viridis)
library(scales)
library(tidyverse)
library(network)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(ELMER)
library(wesanderson)
library(data.table)
library(survival)
library(survminer)
library(rtracklayer)
library(limma)
library(Rsubread)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tracktables)
register(MulticoreParam(8))

setwd("/media/adria/5c3592e6-f482-4fb8-8e84-3c77a7776bd3/BMO_analysis/BMO_output_cell_lines/test")
wd <- "/media/adria/5c3592e6-f482-4fb8-8e84-3c77a7776bd3/BMO_analysis/BMO_output_cell_lines/test/"

#mart hg38
mart <- readRDS(file = paste(wd,"mart.rds", sep = ""))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


############################### TCGA dataset preparation ###############################
load("/media/adria/5c3592e6-f482-4fb8-8e84-3c77a7776bd3/SMARCA_TCGA_chinchu/Glioma_TCGA_Merged.RData")

#select glioma samples where we have RNAseq data
sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)
data<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$HTseq) #obtaining raw counts from TCGA
nm<-sapply(1:length(fd), function(i) names(glioma.tcga[[i]]$mRNAseq[[1]]))
colnames(data)<-unique(substr(unlist(nm),1,12))

# Adjust Gene expression data
##changing RPKM<1 into 1(log0 is not calculatable)
#GENE data filtering for median>1
sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)
data2<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$FPKM) #obtaining the FPKMs
nm<-sapply(1:length(fd), function(i) names(glioma.tcga[[i]]$mRNAseq[[1]]))
colnames(data2)<-unique(substr(unlist(nm),1,12))
a = (data2 < 1)
#create back-up copy
tuned_RPKMs = data2
#changing RPM<1 into 1 (avoiding log 0)
tuned_RPKMs[a] = 1
logs = log2(tuned_RPKMs)

ctrl<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$isControl)
CTRL<-which(ctrl=="YES")
dip<-sapply(1:length(fd), function(i) as.vector(unlist(glioma.tcga[[i]]$clinical))[26]) #taking the grade
dip<-str_replace(dip, "NULL", "CTRL")
dip<- as.factor(dip)
G2<-which(dip=="G2")
G3<-which(dip=="G3")
G4<-which(dip=="G4")
HGG<-which(dip == "G3" | dip == "G4")
all_grades <- which(dip == "G2" | dip == "G3" | dip == "G4")
G2G3 <- which(dip == "G2" | dip == "G3")
non<-setdiff(c(1:677),c(CTRL,G2,G3,G4))

#select only WHO glioma II and IV
a <-which((dip=="G2") | (dip =="G4"))
data <- data[,a]
dip <- dip[a]
dip <- as.factor(dip)
condition <- dip
condition <- as.vector(condition); condition <- as.factor(condition)
condition <- factor(condition, levels = c("G4", "G2"))

#prefilter based on mean (this cutoff can be different)
keep <- rowSums(data) >= 5
data <- data[keep,]
data_matrix <- as.matrix(data)

run_DESeq2 <- function(group_1, group_2){
  print(paste("Running DESeq2 and comparing",group_1, "vs", group_2))
  coldata <- data.frame(row.names=colnames(data), condition)
  print(paste("Number of samples analyzed", dim(coldata)[1], sep = ":"))
  dds <- DESeqDataSetFromMatrix(countData=data_matrix, 
                               colData=coldata, 
                               design=~condition)
  dds <- DESeq(dds)
  levels <- c()
  if (group_1=='G4' | group_2=='G2'){
    levels <- c(levels, 'G4', 'G2')
  }
  vst <- vst(dds, blind = FALSE)
  results <- results(dds, contrast=c("condition", group_1, group_2), independentFiltering = TRUE)
  results_ordered <- data.frame(results[order(results$padj),])
  print(paste("Differentially UP-regulated mRNA with q<0.05: ",sum((results_ordered$padj < 0.05 & results_ordered$log2FoldChange > 0), na.rm=TRUE)))
  print(paste("Differentially DOWN-regulated mRNA with q<0.05: ",sum((results_ordered$padj < 0.05 & results_ordered$log2FoldChange < 0), na.rm=TRUE)))
  write.table(results_ordered, paste0(wd,'DESeq2_', levels[1],'_', 'vs_',levels[2],'.csv'), sep='\t', quote = F)
  return(list(dds, results_ordered, vst))
}

deseq2_comp <- run_DESeq2("G4", "G2")

#sample-to-sample distance - Unpublished Figure
sampleDists <- as.matrix(dist(t(assay(deseq2_comp[[3]]))))
pdf(paste(wd, "TCGA_sample_to_sample_distance.pdf"), height = 10, width = 10)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",col=colorpanel(100, "black", "white"),margin=c(10, 10), main="Sample Distance Matrix", Rowv = TRUE, cexRow = 0.3, cexCol = 0.3)
dev.off()

#Principal component analysis (PCA) in WHO grade IV and WHO grade II from TCGA - Supplementary Figure 2A
pdf(file = paste(wd,"TCGA_PCA_WHO_glioma_II_and_IV.pdf"), width = 10, height = 8)
comp1_PCA <- plotPCA(deseq2_comp[[3]], intgroup = c("condition")) +
  scale_color_manual(values = c("firebrick3", "#1B9E77"), labels=c("GIV", "GII")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14,face = "plain", colour = "black"), 
        axis.title.x = element_text(size = 14,face = "bold"),
        axis.title.y = element_text(size = 14,face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.2)) +
        labs(color='')
comp1_PCA 
dev.off()

#transcriptomic differences between glioma grades (DESeq2) 
results <- deseq2_comp[[2]]
results$ensembl_gene_id <- rownames(results)
biomart_genes <- getBM(rownames(results),  filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "gene_biotype"))
results <- merge(results, biomart_genes, by="ensembl_gene_id")
sum(results$padj < 0.05, na.rm=TRUE)

#volcano plot and DEG between WHO grade IV and WHO grade II from TCGA - Supplementary Figure 2B
topT <- results
pdf(file = paste(wd,"TCGA_volcano_WHO_glioma_II_and_IV.pdf"), width = 10, height = 8)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.2, cex.lab=1.2)
with(topT, plot(log2FoldChange, -log10(padj), col = "gray18",pch=20, main="",xlim=c(-10,10),ylim=c(0,300),cex=1.0, xlab="", ylab=""))
with(subset(topT, padj<=0.05 & abs(log2FoldChange)>=1), points(log2FoldChange, -log10(padj), pch=20, col = ifelse(log2FoldChange >= 0, "firebrick3", "#1B9E77"), cex=0.5))
#grid(nx = NULL, ny = NULL, col = "gray95", lty = 1,lwd = par("lwd"), equilogs = TRUE)
abline(v=0, col="black", lty=3, lwd=1.5)
abline(h = 1.122018, col="black", lty=3, lwd=1.5)
mtext(side=1, line=2, "Log2 Fold Change (GIV/GII)", col="black", font=2,cex=1.2)
mtext(side=2, line=3, "Significance (-Log10 adj. p-value)", col="black", font=2,cex=1.2)
dev.off()

#selecting signifcant genes (up- and down-regulated) with strict thresholds
rownames(results) <- make.names(results[,1], unique = TRUE)
results_pc <- results[results$gene_biotype == "protein_coding",]
results_pc_up <- subset(results_pc, padj < 0.01 & log2FoldChange > 2.5)
results_pc_up <- as.vector(as.numeric(results_pc_up[,9])) #only entrez genes
results_pc_up <- results_pc_up[!is.na(results_pc_up)] #removing NAs (filtered out by DESeq2, outliers or low-expressed features)
results_pc_down <- subset(results_pc, padj < 0.01 & log2FoldChange < -2.5)
results_pc_down <- as.vector(as.numeric(results_pc_down[,9])) 
results_pc_down <- results_pc_down[!is.na(results_pc_down)] 
DEG_all <- c(results_pc_up, results_pc_down)
x <- length(results_pc_up)
y <- length(results_pc_down)
DEG_df <- data.frame(Entrez=DEG_all,group = c(replicate(x, "upregulated"), replicate(y, "downregulated")))
colnames(DEG_df) <- c("entrezgene", "group")

#pathway enrichment analysis - Supplementary Figure 2E and 2F
compKEGG <- compareCluster(entrezgene~group, data=DEG_df,fun = "enrichKEGG", pvalueCutoff=0.05, pAdjustMethod = "BH", organism = "human")
pdf(file = paste(wd,"TCGA_KEGG_strict_cutoff_WHO_glioma_II_and_IV.pdf"), width = 10, height = 10)
KEGG_plot <- dotplot(compKEGG, showCategory = 30, title = "GIV vs GII glioma transcriptomic differences (KEGG)") +
  theme(plot.title = element_text(face="bold"))
KEGG_plot
dev.off()

compGO <- compareCluster(entrezgene~group, data=DEG_df,fun = "enrichGO", pvalueCutoff=0.05, pAdjustMethod = "BH", OrgDb='org.Hs.eg.db',ont="BP")
pdf(file = paste(wd,"TCGA_GO_BP_strict_cutoff_WHO_glioma_II_and_IV.pdf"), width = 12, height = 10)
GO_plot <- dotplot(compGO, showCategory = 30, title = "GIV vs GII glioma transcriptomic differences (GO)") +
  theme(plot.title = element_text(face="bold"))
GO_plot
dev.off()


############################### Randomized patient selection (G2/G4) in TCGA ###############################
data_matrix
data_matrix_G2 <- data_matrix[,coldata == 'G2']
data_matrix_G4 <- data_matrix[,coldata == 'G4']

List = list()
for(a in 1:200){
  print(a)
  #for (i in 1:length(data_matrix)){
  G2 <-data_matrix_G2[,sample(ncol(data_matrix_G2), 20)]
  G4 <-data_matrix_G4[,sample(ncol(data_matrix_G4), 20)]
  pvalues <- sapply(1:nrow(data_matrix), function(i) t.test(as.numeric(G2[i,]),as.numeric(G4[i,]), exact=FALSE)$p.value)
  qvalues <- p.adjust(pvalues, method = "fdr")
  List[[length(List)+1]] = qvalues 
  df_final <- do.call(cbind.data.frame, List)
}

colnames(df_final) <- NULL
df_final2 <- df_final 
df_final2$medians <- rowMeans(as.matrix(df_final2))
rownames(df_final2) <- rownames(data_matrix)
df_final2 <- df_final2[df_final2$medians < 0.05,]
df_final2 <- df_final2[order(df_final2$medians, decreasing = FALSE),]
colnames(df_final2) <- NULL

#select GIV or GII genes based on DESeq2's direcctionality (log2FoldChange)
G4_genes <- merge(results, df_final2, by=0)
G4_genes <- G4_genes[G4_genes$log2FoldChange > 0,]; G4_genes <- G4_genes[,c(2,4,8,9,10,11,212)];colnames(G4_genes)[7] <- "mean"
G4_genes <- G4_genes[order(G4_genes$mean, decreasing = FALSE),];G4_genes <- G4_genes[G4_genes$gene_biotype == "protein_coding",]
G4_genes_later <- G4_genes$ensembl_gene_id;G4_genes <- G4_genes$entrezgene_id;G4_genes <- G4_genes[!is.na(G4_genes)];G4_genes_later <- G4_genes_later[!is.na(G4_genes_later)]
G2_genes <- merge(results, df_final2, by=0)
G2_genes <- G2_genes[G2_genes$log2FoldChange < 0,];G2_genes <- G2_genes[,c(2,4,8,9,10,11,212)];colnames(G2_genes)[7] <- "mean"
G2_genes <- G2_genes[order(G2_genes$mean, decreasing = FALSE),];G2_genes <- G2_genes[G2_genes$gene_biotype == "protein_coding",]
G2_genes_later <- G2_genes$ensembl_gene_id;G2_genes <- G2_genes$entrezgene_id;G2_genes <- G2_genes[!is.na(G2_genes)];G2_genes_later <- G2_genes_later[!is.na(G2_genes_later)]

#load pre-computed G4_genes and G2_genes objects
G4_genes <- readRDS(file = paste(wd,"G4_TCGA_significant_genes.rds", sep = ""))
G2_genes <- readRDS(file = paste(wd,"G2_TCGA_significant_genes.rds", sep = ""))


############################### TFBS prediction on ATAC-seq data (BMO) ###############################
BMO_LN18 <- read.table("ATACseq_LN18_BMO_bound.bed",header = FALSE, sep = "\t") #LN18 GBM cell line
BMO_LN229 <- read.table("ATACseq_LN229_BMO_bound.bed",header = FALSE, sep = "\t") #LN229 GBM cell line
BMO_GBM1 <- read.table("ATACseq_GBM1_BMO_bound.bed",header = FALSE, sep = "\t") #GBM specimen 1
BMO_GBM2 <- read.table("ATACseq_GBM2_BMO_bound.bed",header = FALSE, sep = "\t") #GBM specimen 2

change_colnames <- function(x) {  
  colnames(x) <- c("chr", "start", "end", "TF", "BMO_score", "strand")
  x$chr <- str_remove_all(string = x$chr, pattern = "chr")
  x$strand <- NULL
  return(x)
}

all_BMO <- lapply(list(BMO_LN18, BMO_LN229,BMO_GBM1, BMO_GBM2), change_colnames)

#select predictions that are found to be intersected in both cell lines or in both GBM specimens
all_BMO[[1]]$compare <- paste(all_BMO[[1]]$chr, all_BMO[[1]]$start, all_BMO[[1]]$end, all_BMO[[1]]$TF, sep = ":" )
all_BMO[[2]]$compare <- paste(all_BMO[[2]]$chr, all_BMO[[2]]$start, all_BMO[[2]]$end, all_BMO[[2]]$TF, sep = ":" )
all_BMO[[3]]$compare <- paste(all_BMO[[3]]$chr, all_BMO[[3]]$start, all_BMO[[3]]$end, all_BMO[[3]]$TF, sep = ":" )
all_BMO[[4]]$compare <- paste(all_BMO[[4]]$chr, all_BMO[[4]]$start, all_BMO[[4]]$end, all_BMO[[4]]$TF, sep = ":" )
BMO_LN18_LN229 <- merge(all_BMO[[1]], all_BMO[[2]], by= "compare")
BMO_LN18_LN229 <- BMO_LN18_LN229[,c(2,3,4,5,6,11)]
colnames(BMO_LN18_LN229) <- c("chr", "start", "end", "TF", "BMO_score_LN18", "BMO_score_LN229")
BMO_GBM1_GBM2 <- merge(all_BMO[[3]], all_BMO[[4]], by= "compare")
BMO_GBM1_GBM2 <- BMO_GBM1_GBM2[,c(2,3,4,5,6,11)]
colnames(BMO_GBM1_GBM2) <- c("chr", "start", "end", "TF", "BMO_score_GBM1", "BMO_score_GBM2")
BMO_GBM1_GBM2$chr <-  paste("chr", BMO_GBM1_GBM2$chr, sep = "") 
write.table(BMO_GBM1_GBM2, "common_TF_GBM1_GBM2_coordinates.bed",quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

#intersect TFBS calls with putative glioma enhancers (H3K27ac - Stepniak et al., 2021)
enhancers <- read.table("common_putative_enhancers.bed",header = FALSE, sep = "\t")
colnames(enhancers) <- c("chr", "start", "end"); enhancers[,4] <-NULL
enhancers_GRanges <- makeGRangesFromDataFrame(df = enhancers, seqnames.field = c("chr"), 
                                              keep.extra.columns = FALSE, 
                                              start.field = "start",
                                              end.field = "end")

peakAnno <- annotatePeak(enhancers_GRanges, tssRegion=c(-500, 500),TxDb=txdb, annoDb="org.Hs.eg.db")
COLS <- brewer.pal(n =9, name = 'Paired')

enhancers$chr <- str_replace(enhancers$chr, pattern = "chr", replacement = "")
BMO_LN18_LN229_enhancers <- genome_intersect(BMO_LN18_LN229, enhancers, by=c("chr", "start", "end"), mode = "both") #intersecting based on genomic coordinates

#most frequently TFBS in enhancers
BMO_LN18_LN229_enhancers_table <- table(BMO_LN18_LN229_enhancers$TF)
BMO_LN18_LN229_enhancers_table <- as.data.frame(BMO_LN18_LN229_enhancers_table)
colnames(BMO_LN18_LN229_enhancers_table) <- c("TF", "TF_occurrence_experiment")
BMO_LN18_LN229_enhancers_table$TF <- str_remove(BMO_LN18_LN229_enhancers_table$TF, "_TFBS")
BMO_LN18_LN229_enhancers_table <- BMO_LN18_LN229_enhancers_table[order(BMO_LN18_LN229_enhancers_table$TF_occurrence_experiment, decreasing = TRUE),]

#ENHANCERS: hypergeometric test (x = TFBS_in_enhancers, m = TFBS, n = all_TFBS, k = nr of H3K27ac enhancers)
#what is the probability that my TFBS in enhancers and in TFBS in all the genome (overlap) is any better than obtained by chance alone
total_population <- BMO_LN18_LN229
total_population$TF <- gsub(total_population$TF, pattern = "_TFBS", replacement = "")
hyper_results = list()
for (i in BMO_LN18_LN229_enhancers_table$TF){
  x <- as.numeric(BMO_LN18_LN229_enhancers_table[BMO_LN18_LN229_enhancers_table$TF == i,][2]) 
  m <- dim(total_population[total_population$TF == i,])[1] #nr of times a specific TFBS is found
  n <- 145123 #total number of TFBS that are common in glioblastoma cell lines
  k <- 7571 #total number of TFBS in enhancers
  if (x != 0){
    print(x)
    hyper_results[i] <- unlist(phyper(x - 1, m, n-m, k, lower.tail = FALSE, log.p = FALSE))
  }
}

total_TFBS <- as.data.frame(table(total_population$TF))
colnames(total_TFBS) <- c("TF", "total_TFBS_genome")
hyper_results <- as.data.frame(sort(unlist(hyper_results)))
hyper_results$TF <- rownames(hyper_results)
hyper_results <- merge(hyper_results, BMO_LN18_TF_enhancers_table, by = "TF" )
hyper_results <- merge(hyper_results, total_TFBS, by = "TF" )
hyper_results <- hyper_results[order(hyper_results$`sort(unlist(hyper_results))`, decreasing = FALSE),]
colnames(hyper_results) <- c("TF", "hypergeo_pvalue", "TFBS_occurrence_enhancers", "TFBS_occurrence_genome")
hyper_results$hyperge_adj_pvalue <- p.adjust(hyper_results$hypergeo_pvalue, method = "BH")
hyper_results <- hyper_results[,c(1,3,4,2,5)]

#print final results - Table 1 and Supplementary Table 1
write.table(hyper_results, "enhancers_TFBS_hypergeometic_calculation.bed",quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

#focusing on JUN's binding sites
JUN_enhancers <- as.data.frame(BMO_LN18_LN229_enhancers[grep("JUN_0_A_TFBS", BMO_LN18_LN229_enhancers$TF),])
JUN_enhancers$chr <- paste("chr", JUN_enhancers$chr, sep = "") 
JUN_enhancers <- makeGRangesFromDataFrame(df = JUN_enhancers, seqnames.field = c("chr"), 
                                                 keep.extra.columns = FALSE, 
                                                 start.field = "start",
                                                 end.field = "end")

#TFBS intersection in cell lines and GBM specimens - Figure 1A, upper plot
venn.diagram(list(all_BMO[[1]]$compare, all_BMO[[2]]$compare), category.names = c("TFBS LN18", "TFBS LN229"), output=TRUE, width = 1200 , 
             resolution = 1200, col=c("khaki2","deepskyblue"), filename = "common_TFBS_LN18_LN229.tiff", height = 900, compression = "lzw",
             lwd = 0.7, fill = c(alpha("khaki2",0.5),alpha("deepskyblue", 0.5)), cex = 0.3, fontface = "bold", fontfamily = "arial",
             cat.pos = c(10, 180), cat.cex = 0.3, main.cex = 0.5, cat.dist = c(0.015, 0.015),cat.fontface = "bold", cat.fontfamily = "arial",
             cat.default.pos = "outer", cat.col = c("khaki2", "deepskyblue"), area.vector = 0.5)

venn.diagram(list(all_BMO[[3]]$compare, all_BMO[[4]]$compare), category.names = c("TFBS GBM1", "TFBS GBM2"), output=TRUE, width = 1200, 
             resolution = 1200, col=c("brown1", "darkorchid1"), filename = "common_TFBS_GBM1_GBM2.tiff", height = 900, compression = "lzw",
             lwd = 0.7, fill = c(alpha("brown1",0.5),alpha("darkorchid1", 0.5)), cex = 0.3, fontface = "bold", fontfamily = "arial",
             cat.pos = c(10, 180), cat.cex = 0.3, main.cex = 0.5, cat.dist = c(0.015, 0.015), cat.fontface = "bold", cat.fontfamily = "arial",
             cat.default.pos = "outer", cat.col = c("brown1", "darkorchid1"), area.vector = 0.5)

#link TFBS calls with genes
all_genes <- getBM(mart = mart, attributes = c("chromosome_name","transcription_start_site","hgnc_symbol", "ensembl_gene_id", "gene_biotype", "entrezgene_id"))
all_genes <- all_genes[- grep("CHR*", all_genes$chromosome_name),] 
all_genes$start <- all_genes$transcription_start_site - 1500
all_genes$end <- all_genes$transcription_start_site + 1500
colnames(all_genes)[1] <- "chr"
#all_genes <- all_genes[all_genes$gene_biotype == "protein_coding",] #considering only protein coding
BMO_LN18_LN229_pc_genes <- genome_intersect(BMO_LN18_LN229, all_genes, by=c("chr", "start", "end"), mode = "both") #intersecting based on genomic coordinates
BMO_LN18_LN229_pc_genes <- BMO_LN18_LN229_pc_genes[!duplicated(BMO_LN18_LN229_pc_genes[,c('TF', 'hgnc_symbol')]),] #just considering 1 TF promoter, but we have TF multi-seeding 1 promoter as well
BMO_LN18_LN229_pc_genes <- BMO_LN18_LN229_pc_genes[!BMO_LN18_LN229_pc_genes$hgnc_symbol == "",] #inexistent locus only
BMO_LN18_LN229_pc_genes_tb <- table(BMO_LN18_LN229_pc_genes$TF)
BMO_LN18_LN229_pc_genes_tb <- as.data.frame(BMO_LN18_LN229_pc_genes_tb)
colnames(BMO_LN18_LN229_pc_genes_tb) <- c("TF", "TF_occurrence_experiment")
BMO_LN18_LN229_pc_genes_tb$TF <- str_remove(BMO_LN18_LN229_pc_genes_tb$TF, "_TFBS")
BMO_LN18_LN229_pc_genes_tb <- BMO_LN18_LN229_pc_genes_tb[order(BMO_LN18_LN229_pc_genes_tb$TF_occurrence_experiment, decreasing = TRUE),]

#plot TF models, most common 50 - Figure 1C
paletteLength = length(unique(BMO_LN18_LN229_pc_genes_tb$TF))
library(viridis)
library(scales)
lots_colors <- viridis_pal(option = "B")(50)

BMO_LN18_TF_plot <- BMO_LN18_LN229_pc_genes_tb
BMO_LN18_TF_plot$TF <- as.factor(BMO_LN18_TF_plot$TF)
BMO_LN18_TF_plot$TF <- factor(BMO_LN18_TF_plot$TF, levels = BMO_LN18_TF_plot$TF[1:50])
pdf(file = paste(wd,"TF_models_50_hocomoco.pdf"), width = 10, height = 8)
p <- ggplot(data=BMO_LN18_TF_plot[c(1:50),], aes(x=reorder(TF, TF_occurrence_experiment) , y=TF_occurrence_experiment, width=0.9, fill = TF)) +
  geom_bar(stat="identity", position="stack") +
  xlab("") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("TF binding site instances in LN18/LN229") +
  theme(legend.position = "none") +
  scale_fill_manual(values = lots_colors) 
p + theme_bw() + theme(legend.position = "none", axis.text.y =element_text(hjust = 1, size = 9, face ="bold", colour = "black"), axis.text.x =element_text(face ="bold", colour = "black", size = 12), 
                       axis.ticks = element_line(colour = "black", size = 1), axis.title.x=element_text(face ="bold", colour = "black", size = 12),
                       panel.border = element_rect(colour = "black", fill=NA, size=1.5))
dev.off()


############################### GIV and GII genes, obtain their TSS and promoter and intersect with TFBS predictions ###############################
get_TSS_TFBS <- function(genes, grade){ #dataframe with entrez id's
  genes <- getBM(genes$entrezgene_id, 
                      filters = "entrezgene_id", 
                      mart = mart, 
                      attributes = c("chromosome_name","transcription_start_site","hgnc_symbol", "ensembl_gene_id", "gene_biotype", "entrezgene_id"))
  genes <- genes[!duplicated(genes[,c('hgnc_symbol')]),] 
  genes$start <- genes$transcription_start_site - 1500 #3kB window
  genes$end <- genes$transcription_start_site + 1500
  colnames(genes)[1] <- "chr"
  genes <- genome_intersect(genes, BMO_LN18_LN229, by=c("chr", "start", "end"))
  genes <- as.data.frame(genes)
  genes  <- genes[!duplicated(genes[,c('ensembl_gene_id', 'TF')]),] #remove TF multiseed, 1 TFBS - 1 specific promoter
  genes$TF <- str_remove(genes$TF, "_TFBS")
  genes <- merge(genes, BMO_LN18_LN229_pc_genes_tb, by = "TF") #add TF occurrence
  genes$BMO_score_mean <- apply(genes[,c(8,9)],1,mean) #order by the BMO score
  genes  <- genes[order(genes$BMO_score_mean, decreasing = TRUE),]
  FREQ <- table(genes$TF)
  FREQ <- as.data.frame(FREQ)
  colnames(FREQ)[1] <- "TF"
  genes <- merge(genes, FREQ, by="TF")
  if (grade == "G2"){
    print(grade)
    genes$Freq_norm <- (genes$Freq * 100) / dim(G2_genes)[1]  #normalize by geneset number (DEGs in GII / GIV)
    genes <- genes[order(genes$Freq_norm, decreasing = TRUE),]
    return(genes)
  } else {
    print(grade)
    genes$Freq_norm <- (genes$Freq * 100) / dim(G4_genes)[1]  #normalize by geneset number (DEGs in GIV / GII)
    genes <- genes[order(genes$Freq_norm, decreasing = TRUE),]
    return(genes)
}}

#obtain TFBS predictions in promoters of overexpressed genes only (GIV and GII genes)
G4_genes_TSS_BMO <- get_TSS_TFBS(G4_genes, "G4")
G2_genes_TSS_BMO <- get_TSS_TFBS(G2_genes, "G2")

#intersect with GBM cell lines LN18 and LN229 RNA-sequencing; calculate RPKMs values
CL_expression <- read.table("SAR_LN18_LN229_WT.featureCounts",header = TRUE, sep = "\t")
raw_to_RPKM <- function(raw_data){
  exon_length <- raw_data$Length
  raw_data <- raw_data[c(1,7,8,9,10)]
  colnames(raw_data) <- c("ensembl_gene_id", "LN18_1", "LN18_2", "LN229_1", "LN229_2")
  dgList <- DGEList(counts=as.matrix(raw_data[,-1]), genes=raw_data$ensembl_gene_id)
  RPKM <- rpkm(dgList, normalized.lib.sizes=TRUE, gene.length = exon_length)
  RPKM <- as.data.frame(RPKM, row.names = raw_data$ensembl_gene_id)
  RPKM$ensembl_gene_id <- rownames(RPKM)
  return(RPKM)
}

RPKM <- raw_to_RPKM(CL_expression)

G2_genes_TSS_BMO <- merge(G2_genes_TSS_BMO, RPKM, by="ensembl_gene_id")
G2_genes_TSS_BMO <- G2_genes_TSS_BMO[order(G2_genes_TSS_BMO$Freq_norm, decreasing = TRUE),]
G4_genes_TSS_BMO <- merge(G4_genes_TSS_BMO, RPKM, by="ensembl_gene_id")
G4_genes_TSS_BMO <- G4_genes_TSS_BMO[order(G4_genes_TSS_BMO$Freq_norm, decreasing = TRUE),]

#check different TFBS frequencies between GIV and GII
G2_genes_TSS_BMO_TF  <- G2_genes_TSS_BMO[!duplicated(G2_genes_TSS_BMO[,c('TF', 'Freq_norm')]),]
G2_genes_TSS_BMO_TF <- G2_genes_TSS_BMO_TF[,c(2,15)]
G4_genes_TSS_BMO_TF  <- G4_genes_TSS_BMO[!duplicated(G4_genes_TSS_BMO[,c('TF', 'Freq_norm')]),] 
G4_genes_TSS_BMO_TF <- G4_genes_TSS_BMO_TF[,c(2,15)]
final <- merge(G4_genes_TSS_BMO_TF,G2_genes_TSS_BMO_TF, by = "TF", all=T)
final_G4 <- final[is.na(final$Freq_norm.y),]
final_G2 <- final[is.na(final$Freq_norm.x),]

#A) Grade-specific TFBS - Figure 1D, right barplot
final_G4$condition <- "G4 upregulation"
final_G2$condition <- "G2 upregulation"
zz <- rbind(final_G2, final_G4)
zz$Frequency <- c(zz$Freq_norm.y[c(1:16)], zz$Freq_norm.x[c(17:82)]) #transform df to a new column
zz$TF <- as.factor(zz$TF)
zz <- zz[order(zz[,4], zz[,5]),]
zz$TF <- factor(zz$TF, levels = zz$TF)

plot_grade_specific <- ggplot(data=zz, aes(x=TF, y=Frequency, fill=condition,width=0.7)) +
  geom_bar(stat="identity", position="stack") +
  xlab("") +
  ggtitle("Grade-specific TFBS") +
  coord_flip() +
  ylab("Normalized TF binding site frequency") +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#1B9E77", "firebrick3")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(title = element_text(size = 10, face = "bold", hjust = 1), axis.text.y =element_text(hjust = 1, size = 9, face ="bold", colour = "black"), axis.text.x =element_text(face ="bold", colour = "black", size = 12), 
        axis.ticks = element_line(colour = "black", size = 1), axis.title.x=element_text(face ="bold", colour = "black", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5), panel.background = element_blank()) 

#B) Generic TFBS - found in both GIV and GII gliomas with similar frequencies - Figure 1D, left barplot
final <- final[complete.cases(final),] #remove NAs
final <- melt(final, id.vars = c("TF"))
final$condition <- "G4 upregulation"
final$condition[589:1178] <- "G2 upregulation"

#TF ordering in stacked bar plot
new_order <- final %>%
  group_by(TF) %>%
  summarise(var1=sum(value))

new_order <- new_order[order(new_order$var1, decreasing = FALSE),]; new_order <- as.data.frame(new_order)
final$TF <- factor(final$TF,levels = new_order$TF)
zzz <- final
#select only top TF
a <- as.list(tail(new_order, n=100)[1]); a<-unlist(a)
zzz <- zzz[zzz$TF %in% a,]

plot_generic <- ggplot(data=zzz, aes(x=TF, y=value, fill=condition,width=0.8)) +
  geom_bar(stat="identity", position="stack") +
  xlab("") +
  ggtitle("Generic TFBS") +
  coord_flip() +
  ylab("Normalized TF binding site frequency") +
  theme(legend.position="none") +
  scale_fill_manual(values = c("#1B9E77", "firebrick3")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(title = element_text(size = 10, face = "bold", hjust = 1),axis.text.y =element_text(hjust = 1, size = 9, face ="bold", colour = "black"), axis.text.x =element_text(face ="bold", colour = "black", size = 12), 
        axis.ticks = element_line(colour = "black", size = 1), axis.title.x=element_text(face ="bold", colour = "black", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5), panel.background = element_blank())

#plot TFBS found in overexpessed genes of both grades and unique - Figure 1D
pdf(file = paste(wd,"Figure_1D.pdf"), width = 10, height = 12)
ggarrange(plot_generic,plot_grade_specific)
dev.off()


############################### Gene-TF association of grade-specific TFBS ###############################
genes_pathways_G2 <- G2_genes_TSS_BMO[G2_genes_TSS_BMO[,2] %in% final_G2[,1],] #selecting 
genes_pathways_G4 <- G4_genes_TSS_BMO[G4_genes_TSS_BMO[,2] %in% final_G4[,1],]
network_G2_TF <- genes_pathways_G2[,c(2,5)]
network_G2_TF <- tibble::as_tibble(network_G2_TF)
network_G4_TF <- genes_pathways_G4[,c(2,5)]
network_G4_TF <- tibble::as_tibble(network_G4_TF)

gene_TF_network <- function(genes_TFBS){ #link TFBS with their genes in a network
  
  sources <- genes_TFBS %>%
    dplyr::distinct(TF) %>%
    dplyr::rename(label = TF)
  
  destinations <- genes_TFBS %>%
    dplyr::distinct(hgnc_symbol) %>%
    dplyr::rename(label = hgnc_symbol)
  
  nodes <- full_join(sources, destinations, by = "label")
  nodes <- nodes %>% tibble::rowid_to_column("id")
  
  per_route <- genes_TFBS %>%  
    group_by(TF, hgnc_symbol) %>%
    summarise(weight = n()) %>% 
    ungroup()
  
  edges <- per_route %>% 
    left_join(nodes, by = c("TF" = "label")) %>% 
    rename(from = id)
  
  edges <- edges %>% 
    left_join(nodes, by = c("hgnc_symbol" = "label")) %>% 
    rename(to = id)
  
  edges <- select(edges, from, to, weight)
  
  routes_network <- network(edges, vertex.attr = nodes, matrix.type = "edgelist", ignore.eval = FALSE, directed = TRUE)
  network.vertex.names(routes_network) <- nodes$label
  #plot(routes_network, vertex.cex = 2, displaylabels = TRUE)
  
  #change color by group
  routes_network %v% "Family" = ifelse(network.vertex.names(routes_network) %in% genes_TFBS$TF, "TF", 
                                       ifelse(network.vertex.names(routes_network) %in% genes_TFBS$hgnc_symbol, "gene", "other"))
  
  # modify color
  routes_network %v% "color" = ifelse(routes_network %v% "Family" == "TF", "violet", 
                                      ifelse(routes_network %v% "Family" == "gene", "gold1","gray60"))
  # rescale edge size
  network::set.edge.attribute(routes_network, "weights", ifelse(routes_network %e% "weights" <= 1, 0.25, 
                                                                ifelse(routes_network %e% "weights" <= 3, .5, 1)))
  # define line type
  network::set.edge.attribute(routes_network, "lty", ifelse(routes_network %e% "weights" == 0.25, 3, 
                                                            ifelse(routes_network %e% "weights" == .5, 2, 1)))
  
  ggnet2(routes_network, 
         color = "color", 
         label = TRUE,arrow.size = 7, arrow.gap = 0.01,arrow.type = "open",
         label.size = 3,
         alpha = 0.95,
         size = "Family",
         edge.alpha = 1, edge.label.fill = "black", edge.label.size = 15) +
    guides(color = FALSE, size = FALSE)
}

#networks - Supplementary Figure 7
pdf(paste(wd, "Supplementary_Figure7.pdf"), height = 12, width = 12)
gene_TF_network(network_G2_TF)
gene_TF_network(network_G4_TF)
dev.off()

############################### GIV specific TFBS and theirs genes ###############################
final_G4 <- readRDS(file = paste(wd,"final_G4.rds", sep = "")) #TFBS identified only in genes overexpressed in GBM
G4_genes_TSS_BMO <- readRDS(file = paste(wd,"G4_genes_TSS_BMO.rds", sep = "")) #df with all information, previous mart object

genes_pathways <- G4_genes_TSS_BMO[G4_genes_TSS_BMO$TF %in% final_G4$TF,]
genes_pathways <-getBM(values = unique(genes_pathways$ensembl_gene_id),  filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"))
REACTOME <- enrichPathway(genes_pathways$entrezgene_id, pvalueCutoff=0.05, qvalueCutoff = 0.05,pAdjustMethod = "BH" , organism = "human")

#plotting pathways related with G4 genes - Figure 3C
pdf(paste(wd, "Figure3C.pdf", height = 8, width = 12))
REACTOME_plot <- dotplot(REACTOME, showCategory = 12)
REACTOME_plot + theme(panel.grid.major = element_line(colour = "azure3"),panel.border = element_rect(colour = "black", fill=NA, size=1.5),
                      panel.spacing = unit(1, "lines"), strip.background = element_rect(colour = "black", fill = "gray20"), 
                      axis.text.y = element_text(hjust=1, size = 12,color = "black", face = "bold"),axis.text.x=element_text(angle=0, size = 12,color = "black", face = "bold"), 
                      legend.position = "top", legend.justification = "top", legend.title = element_text(face = "bold"), 
                      strip.text.x = element_text(colour = "white", face = "bold", size = 12)) 
dev.off()

############################### Gene set enrichment with selected GIV and GII genes ###############################
#gse input
results <- deseq2_comp[[2]]
results$ensembl_gene_id <- rownames(results)
biomart_genes <- getBM(rownames(results),  filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "gene_biotype"))
results <- merge(results, biomart_genes, by="ensembl_gene_id")
sum(results$padj < 0.05, na.rm=TRUE)
deseq2_results <- results

original_gene_list <- deseq2_results[deseq2_results$padj < 0.05,]
original_gene_list <- deseq2_results$log2FoldChange
names(original_gene_list) <- deseq2_results$ensembl_gene_id
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(gene_list)] #remove duplicates
#gene_list <- gene_list[names(gene_list) %in% genes_pathways$ensembl_gene_id]
write.table(gene_list, paste(wd, "genes_pathways_GSEA_FC.bed", sep = ""),quote = FALSE, sep = "\t",row.names = TRUE, col.names = FALSE)

#GO
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 5, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

quant10 <- as.numeric(quantile(gse@result$enrichmentScore, probs = c(0.05, 0.95))[1])
quant90 <- as.numeric(quantile(gse@result$enrichmentScore, probs = c(0.05, 0.95))[2])
gse@result <- gse@result[(gse@result$enrichmentScore < quant10 | gse@result$enrichmentScore > quant90),]

#plotting GSE pathways related with G4 genes - Supplementary Figure 2D
pdf(paste(wd, "Supplementary_Figure2D.pdf", sep = ""), height = 10, width = 12)
GO_gsea<- ridgeplot(gse, fill = "enrichmentScore") + labs(x = "enrichment distribution") + scale_fill_gradient(low = "royalblue1", high = "red1") + 
  theme(legend.text=element_text(size=12, face = "bold"),panel.grid.major = element_line(colour = "azure3"),panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        panel.spacing = unit(1, "lines"), strip.background = element_rect(colour = "black", fill = "gray20"), 
        axis.text.y = element_text(hjust=1, size = 12,color = "black", face = "bold"),axis.text.x=element_text(angle=0, size = 12,color = "black", face = "bold"), 
        legend.position = "right", legend.justification = "right", legend.title = element_text(size = 12,face = "bold"), axis.title.x =element_text(size = 14,color = "black", face = "bold"), 
        strip.text.x = element_text(colour = "white", face = "bold", size = 12)) + geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.8)
GO_gsea
dev.off()

#change format
deseq2_results <- results
original_gene_list <- deseq2_results[deseq2_results$padj < 0.05,]
original_gene_list <- deseq2_results$log2FoldChange
names(original_gene_list) <- deseq2_results$entrezgene_id
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list = gene_list[!duplicated(gene_list)] #remove duplicates

#KEGG (nPerm parameter avoids the error!)
gse2 <- gseKEGG(geneList = gene_list,
                organism = "hsa",
                minGSSize = 5,
                nPerm = 10000, 
                maxGSSize = 800,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                keyType = "ncbi-geneid",
                verbose = TRUE)

require(DOSE)
library(viridis)
dotplot(gse2, showCategory=10, split=".sign") + facet_grid(.~.sign)

#calculate quantiles and apply thresholds to select pathways
quant25 <- as.numeric(quantile(gse2@result$enrichmentScore)[2])
quant75 <- as.numeric(quantile(gse2@result$enrichmentScore)[4])
gse2@result <- gse2@result[(gse2@result$enrichmentScore < quant25 | gse2@result$enrichmentScore > quant75),]

#plotting GSE pathways related with G4 genes - Supplementary Figure 2C
pdf(paste(wd, "Supplementary_Figure2C.pdf", sep = ""), height = 10, width = 10)
KEGG_gsea<- ridgeplot(gse2, fill = "enrichmentScore") + labs(x = "enrichment distribution") + scale_fill_gradient(low = "royalblue1", high = "red1") + 
  theme(legend.text=element_text(size=12, face = "bold"),panel.grid.major = element_line(colour = "azure3"),panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        panel.spacing = unit(1, "lines"), strip.background = element_rect(colour = "black", fill = "gray20"), 
        axis.text.y = element_text(hjust=1, size = 12,color = "black", face = "bold"),axis.text.x=element_text(angle=0, size = 12,color = "black", face = "bold"), 
        legend.position = "right", legend.justification = "right", legend.title = element_text(size = 12,face = "bold"), axis.title.x =element_text(size = 14,color = "black", face = "bold"), 
        strip.text.x = element_text(colour = "white", face = "bold", size = 12)) + geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.8)
KEGG_gsea
dev.off()

############################### TCGA mRNA expression of genes coding for TF which TFBS are enriched in G4 genes ###############################
#directly substract TF names: some names do not concorde between TF from hocomoco db and the protein coding gene, done it manually in 657
TF_expression <- levels(zz$TF)
TF_expression <- rev(TF_expression)
TF_expression <- str_remove(TF_expression, pattern = "_._.")

#manually modify gene ids for some TF, different names and alliase: eliminated are DUXA, HSFY1, SOX10 due to pre-filtering
TF_expression <- c("JUN", "SCRT2", "PITX3", "ESRRA", "ZN784", "PBX2", "FOXH1", "ZFP28", "NRF1", "EN2", "TEF", "SOX8", "ZNF146", "NR6A1", "MAFF", 
                   "HNF1B", "GATA2", "FUBP1", "BARX1", "PITX2", "LBX2", "HOXD11", "HOXC6", "HOXC11", "HSF2", "HBP1", "GATA1", "FOXO4", "FOXB1", "DDIT3",
                   "ZNF502", "ZFP82", "TWIST1", "POU4F2", "MSX1", "MEOX2", "HOXD9", "HNF1B", "BARX2", "BARHL1", "ZNF394", "VAX2", "SOX21",
                   "SOX11", "POU5F1", "PDX1", "OVOL1", "OLIG3", "OLIG1", "NKX6-2", "MEIS1", "MAFG", "HOXC10", "HOXB3", "HOXA2", "HOXA1", "HMX1", 
                   "GATA6", "FOXL1", "MECOM", "DMBX1", "DLX4", "CUX2", "BHLHE22")

#JUN TARGETS
gene_list <- c("IFRD1", "VIM", "FOSL2", "RAB36", "PTN", "SPATA1", "TMEM43", "SLFN12", "TRIB1", "RIN1", "GRP3", "SIAH2", "UPP1", "FAM111B", "S100A2", "S100A10")

#G4 TF genes
anno <- as.data.frame(colData(deseq2_comp[[3]])[, c("condition")])
rownames(anno) <- colnames(assay(deseq2_comp[[3]]))
colnames(anno) <- c("condition")

#three different genesets to plot
DEG <- deseq2_results[deseq2_results$hgnc_symbol %in% TF_expression, ] #G4 specific TFBS
DEG <- deseq2_results[deseq2_results$hgnc_symbol %in% genes_pathways$hgnc_symbol, ] #genes that might be regulated by G4 specific TFBS
DEG <- deseq2_results[deseq2_results$hgnc_symbol %in% gene_list, ] #c-Jun targets
mat <- assay(deseq2_comp[[3]])[DEG$ensembl_gene_id, ]
a <- getBM(rownames(mat),  filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id", "external_gene_name"))
rownames(mat) <- make.names(a$external_gene_name, unique = TRUE)

my_colour = list(condition = c("G4" = "firebrick3", "G2" = "#1B9E77"))
#selecting extra clinical data from Verhaak et at., Cancer Cell (2010)
TCGA_clinical <- read.table("TableS1.PatientData.20151020.v3.csv",header = TRUE, sep = "\t")
TCGA_clinical <- TCGA_clinical[,c(1,13,14,15,16,22)]
rownames(TCGA_clinical) <- TCGA_clinical[,1]
TCGA_clinical[,1] <- NULL
#select patients in mat, so focusing on G2 and G4
TCGA_clinical <- TCGA_clinical[rownames(TCGA_clinical) %in% colnames(mat),]
mat <- mat[,colnames(mat) %in% rownames(TCGA_clinical)]
#reorder
TCGA_clinical <- TCGA_clinical[order(row.names(TCGA_clinical)),]
mat <- mat[,rownames(TCGA_clinical)]
colnames(TCGA_clinical) <- c("Histology", "Grade", "Age", "Gender", "IDH_status")
TCGA_clinical$Grade <- gsub(x = TCGA_clinical$Grade, pattern = "G4", replacement = "GIV")
TCGA_clinical$Grade <- gsub(x = TCGA_clinical$Grade, pattern = "G2", replacement = "GII")
TCGA_clinical$Histology <- factor(TCGA_clinical$Histology, levels = c("oligodendroglioma", "oligoastrocytoma", "astrocytoma", "glioblastoma"))

#palete colours
f1 = colorRamp2(seq(min(mat), max(mat), length = 3), c("blue", "#EEEEEE", "red"))

#annotations
ann <- TCGA_clinical
colours <- list('Histology' = c('glioblastoma' = 'darkorchid4', 'astrocytoma' = 'orange3', 'oligoastrocytoma' = 'goldenrod1', 'oligodendroglioma' = 'lightgoldenrod1'),
                'Grade' = c('GII' = '#1B9E77', 'GIV' = 'firebrick3'),
                'Age' = colorRamp2(c(0, 90), c("white", "#007700")),
                'Gender' =  c('female' = 'dodgerblue4', 'male' = 'darkorange'),
                'IDH_status' = c('WT' = 'seashell2', 'Mutant' = 'gray7'))

colAnn <- HeatmapAnnotation(df = TCGA_clinical,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

#simple_anno_size = unit(1, "cm")

#a) plotting TF coding genes expression using TCGA data - Figure 3A
scaled_mat = t(scale(t(mat)))
ht <- Heatmap(scaled_mat, name = "z-score",
              column_title = "TCGA samples", 
              column_title_side = "bottom",
              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              row_title = "", 
              row_title_gp = gpar(fontsize = 20, fontface = "bold"),
              rect_gp = gpar(col = "white", lwd = 0),
              cluster_rows = TRUE, cluster_columns = TRUE, show_column_dend = TRUE, show_row_dend = TRUE,
              clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
              row_names_gp = gpar(fontsize = 8), show_column_names = FALSE, row_km = 2, column_km = 2,
              top_annotation = colAnn, border = TRUE)

tiff(paste(wd, "grade_specific_TFBS_TCGA_expression.tff", sep = ""), units="in", width=12, height=10, res=500, compression = 'lzw')
draw(ht)
dev.off()

#b) plotting c-Jun targets (promoter) using TCGA data - Figure 2D
scaled_mat = t(scale(t(mat)))
ht <- Heatmap(scaled_mat, name = "z-score", 
              column_title = "", 
              column_title_side = "bottom",
              column_title_gp = gpar(fontsize = 15, fontface = "bold"),
              row_title = "",
              row_title_gp = gpar(fontsize = 20, fontface = "bold"),
              rect_gp = gpar(col = "white", lwd = 0),
              cluster_rows = TRUE, cluster_columns = TRUE, show_column_dend = FALSE, show_row_dend = TRUE,
              clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
              row_names_gp = gpar(fontsize = 8), show_column_names = FALSE,
              top_annotation = colAnn, border = TRUE, heatmap_height = unit(10, "cm"))

tiff(paste(wd, "TCGA_cJun_targets_expression.tff", sep = ""), units="in", width=6, height=10, res=500, compression = 'lzw')
draw(ht)
dev.off()

#c) plotting all genes targeted by GIV-specific TFs using TCGA data - Supplementary Figure 5A
scaled_mat = t(scale(t(mat)))
ht <- Heatmap(scaled_mat, name = "z-score", 
              column_title = "", 
              column_title_side = "bottom",
              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              row_title = "",
              row_title_gp = gpar(fontsize = 20, fontface = "bold"),
              rect_gp = gpar(col = "white", lwd = 0),
              cluster_rows = TRUE, cluster_columns = TRUE, show_column_dend = FALSE, show_row_dend = TRUE,
              clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
              row_names_gp = gpar(fontsize = 6), show_column_names = FALSE,
              top_annotation = colAnn, border = TRUE) + 
  rowAnnotation(link = anno_mark(at = which(rownames(scaled_mat) %in% gene_list), 
                                 labels = rownames(scaled_mat)[which(rownames(scaled_mat) %in% gene_list)], 
                                 labels_gp = gpar(fontsize = 10), padding = unit(2, "mm")))    

tiff(paste(wd, "grade_specific_TFBS_genes_TCGA_expression.tff", sep = ""), units="in", width=12, height=10, res=500, compression = 'lzw')
draw(ht)
dev.off()

############################### Gene expression correlation in the TCGA ###############################
vst_data <- assay(deseq2_comp[[3]])

#JUN correlation in vst object, preparing df
JUN <- G4_genes_TSS_BMO[G4_genes_TSS_BMO$TF == "JUN_0_A",][,c(1,2)] 
rownames(JUN) <- JUN$ensembl_gene_id
JUN$TF <- "ENSG00000177606" #ensembl id
test <- getBM(JUN$ensembl_gene_id, filters = "ensembl_gene_id", mart = mart, attributes = c("hgnc_symbol"))
JUN <- cbind(JUN,test)

#plotting Pearsons correlation and its significance between JUN gene and the genes it targets using TCGA data - Supplementary Figure 5A
names_final <- "JUN"
pdf(paste(wd, "Supplementary_Figure5A.pdf", sep = ""), height = 10, width = 12)
par(mfrow=c(1,1))
JUN_ordered <- JUN[c(3,2,9,16,13,6,15,12,7,11,1,8,10,4,14,5),] #re-order
for (i in 1:length(rownames(JUN_ordered))){
  vst_data_1 <- rbind(vst_data[(rownames(vst_data) %in% JUN_ordered$ensembl_gene_id[i]),], vst_data[(rownames(vst_data) %in% JUN_ordered$TF[i]),])
  res1 <- cor.mtest(t(vst_data_1), conf.level=0.95) # by default is pearson
  M = cor(t(vst_data_1))
  print(c(JUN$hgnc_symbol[i], M[1,2]))
  colnames(M) <- c(names_final[1], JUN$hgnc_symbol[i])
  rownames(M) <- c(names_final[1], JUN$hgnc_symbol[i])
  corrplot(M, p.mat = res1$p, method = "pie", insig = "label_sig", sig.level = c(.0001, .001, .005), pch.cex = 2.2, tl.cex = 1.7, 
           pch.col = "black",na.label = "NA", tl.srt = 0, tl.col = "black", cl.pos = "r", cl.cex = 1.5,cl.align.text = "c", type = "upper",
           cl.ratio = 0.3,cl.length = 5,col = brewer.pal(n = 8, name = "RdYlBu"),diag = TRUE)
}
dev.off()

correlation_mRNA <- function(genes, grade){
  AP1_complex_genes <- genes
  AP1_complex_genes <- as.data.frame(AP1_complex_genes)
  colnames(AP1_complex_genes) <- "hgnc_symbol"
  AP1_complex_genes <- getBM(AP1_complex_genes$hgnc_symbol, filters = "hgnc_symbol", mart = mart, attributes = c("hgnc_symbol", "ensembl_gene_id"))
  AP1_complex <- AP1_complex_genes
  
  #split the dataset depending on the grade
  if (grade=='G2'){
    vst_data <- assay(deseq2_comp[[3]])
    grade <- dip == "G2"
    vst_data <- vst_data[,grade]
    
  } else {
    vst_data <- assay(deseq2_comp[[3]])
    grade <- dip == "G4"
    vst_data <- vst_data[,grade]
  }
  AP1_complex <- vst_data[(rownames(vst_data) %in% AP1_complex$ensembl_gene_id),]
  AP1_complex <- AP1_complex[c(7,3,4,5,1,6,2),]
  res1 <- cor.mtest(t(AP1_complex), conf.level=0.95) 
  M = cor(t(AP1_complex))
  colnames(M) <- rownames(M) <- c("JUN", "JUND", "FOS", "JUNB", "FOSL2", "FOSL1", "FOSB")
  rownames(M) <- c("JUN", "JUND", "FOS", "JUNB", "FOSL2", "FOSL1", "FOSB")
  return(corrplot(M, p.mat = res1$p, method = "pie", insig = "label_sig", sig.level = c(.0001, .001, .005), pch.cex = 2.2, tl.cex = 1.7, 
                    pch.col = "black",na.label = "NA", tl.srt = 0, tl.col = "black", cl.pos = "r", cl.cex = 1.5,cl.align.text = "c", type = "upper",
                    cl.ratio = 0.3,cl.length = 5,col = brewer.pal(n = 8, name = "RdYlBu"),diag = TRUE))
  
}

#plotting Pearsons correlation and its significance between JUN gene and AP1 partners - Supplementary Figure 6B, first row of data
pdf(paste(wd, "Supplementary_Figure6B.pdf", sep = ""), height = 10, width = 12)
correlation_mRNA((c("FOS", "FOSL1", "FOSL2", "FOSB", "JUNB", "JUND", "JUN")), "G4")
correlation_mRNA((c("FOS", "FOSL1", "FOSL2", "FOSB", "JUNB", "JUND", "JUN")), "G2")
dev.off()

#AP-1 complex together with JUN
genes <- getBM(values = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"),  filters = "hgnc_symbol", mart = mart, attributes = c("ensembl_gene_id","hgnc_symbol" )) 

# keep in mind there are differences between the TF and the gene codigin e.g. EN / HME2
compare_TCGA_expression <- function(gene){
  genes <- getBM(values = gene,  filters = "hgnc_symbol", mart = mart, attributes = c("ensembl_gene_id","hgnc_symbol")) 
  DFRAME_0 = data.frame()
  for (i in 1:length(genes$hgnc_symbol)){
    expression <- logs[genes$ensembl_gene_id[i],]
    expression <- as.data.frame(expression)
    expression$group <- "group"
    expression$group[non] <- "non"
    expression$group[CTRL] <- "CTRL"
    expression$group[G2] <- "GII"
    expression$group[G3] <- "GIII"
    expression$group[G4] <- "GIV"
    expression <- expression[!expression$group == "non",]
    expression <- expression[!expression$group == "CTRL",]
    expression <- expression[expression$group == "GII" | expression$group == "GIII" | expression$group == "GIV",]
    expression$variable <- genes$hgnc_symbol[i]
    DFRAME_0 <- rbind(DFRAME_0,expression)
  }
  
  colores <- c("#1B9E77","chocolate3", "firebrick3")
  names(colores) = c("GII", "GIII", "GIV")
  my_comparisons <- list(c("GII", "GIV"),c("GII", "GIII"),c("GIII", "GIV") )
  
  DFRAME_0$variable <- factor(DFRAME_0$variable, levels = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"))
  
  p <- ggboxplot(DFRAME_0, x = "group", y = "expression", facet.by = "variable", bxp.errorbar = TRUE,size = 0.6, ggtheme = theme_bw(),
                 color = "group", palette = colores,add="jitter",
                 order = c("GII", "GIII", "GIV"),
                 ylab = "log2(FPKM)", xlab = "",
                 title = "",
                 legend = "right",
                 outline = TRUE,
                 outlier.shape =NA, notch = TRUE) +
    stat_compare_means( comparisons = my_comparisons, method = "wilcox.test",size = 9, vjust = 0.4,hide.ns = FALSE,
                        symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))) +
    font("ylab", size = 18, color = "black", face = "bold") +
    scale_y_continuous(limits = c(0,13)) +
    theme(strip.text = element_text(size = 14, face ="bold", colour = "black"),axis.text.y =element_text(size = 14, face ="bold", colour = "black"), axis.text.x =element_text(face ="bold", colour = "black", size = 12), 
          axis.ticks = element_line(colour = "black", size = 1), axis.title.x=element_text(face ="bold", colour = "black", size = 14),
          panel.background = element_blank()) + geom_hline(yintercept = 5, linetype = "dashed", color = "black", size=0.5)
  return(p)
}

#plotting JUN gene and AP1 partners expression across glioma grades using TCGA data - Supplementary Figure 6A
pdf(paste(wd, "Supplementary_Figure6A.pdf", sep = ""), onefile = T, height = 5, width = 4)
compare_TCGA_expression("JUN")
compare_TCGA_expression("JUND")
compare_TCGA_expression("FOS")
compare_TCGA_expression("JUNB")
compare_TCGA_expression("FOSL2")
compare_TCGA_expression("FOSL1")
compare_TCGA_expression("FOSB")
dev.off()


############################### Glioma LN18 and LN229 cell line mRNA expression of GIV TFBS ###############################
genes <- getBM(values = TF_expression,  filters = "hgnc_symbol", mart = mart, attributes = c("ensembl_gene_id","hgnc_symbol" )) 
genes <- genes[!duplicated(genes[,2], fromLast=T),] #there are 2 ensembl_gene_id per hugo symbol
DFRAME_0 = data.frame()
for (i in 1:length(genes$hgnc_symbol)){
  print(i)
  print(genes$ensembl_gene_id[i])
  expression <- logs[genes$ensembl_gene_id[i],]
  expression <- as.data.frame(expression)
  expression$group <- "group"
  expression$group[non] <- "non"
  expression$group[CTRL] <- "CTRL"
  expression$group[G2] <- "G2"
  expression$group[G3] <- "G3"
  expression$group[G4] <- "G4"
  expression <- expression[!expression$group == "non",]
  expression <- expression[!expression$group == "CTRL",]
  expression <- expression[expression$group == "G2" | expression$group == "G4",]
  expression$variable <- genes$hgnc_symbol[i]
  DFRAME_0 <- rbind(DFRAME_0,expression)
}

DFRAME_CL <- merge(RPKM, genes,by = "ensembl_gene_id")

#reformatting of the df
f<-data.frame(c(DFRAME_CL[,2], DFRAME_CL[,3], DFRAME_CL[,4], DFRAME_CL[,5]))
f$group <- c("CL")
f$variable <- rep(DFRAME_CL$hgnc_symbol, 4)
colnames(f)[1] <- "expression"
f$expression <- log2(f$expression)

#low values and negative to 0, non expressed TF (approximation)
f$expression <- ifelse(f$expression <0, 0,  f$expression)

#merge TCGA and cell line
DFRAME_0_f <- rbind(DFRAME_0, f)
colores <- c("#1B9E77", "firebrick3", "black")
names(colores) = c("G2", "G4", "CL")
my_comparisons <- list(c("G4", "G2"))
colnames(DFRAME_0_f) <- c("expression", "group", "TF")

#order the TF list based on the heatmap from Figure 3A
new_levels <- c("HOXD11","HOXD9","HOXC10","HOXC11","HOXC6","HOXB3","HOXA2","HOXA1","MEOX2","MECOM","TWIST1",
                  "FOXL1","BARX1","DMBX1","GATA6","MAFF","JUN","DDIT3","ZNF502","PDX1","OLIG3","POU5F1","MEIS1",
                  "EN2","PITX2","FOXB1","MSX1","LBX2","ZNF394","HBP1","ESRRA", "PITX3","FOXH1","OVOL1","GATA1",
                  "GATA2","TEF","CUX2","FOXO4","HNF1B", "BARX2","NKX6-2","BHLHE22","MAFG","BARHL1","SCRT2","SOX11",
                  "HSF2","OLIG1","SOX8", "HMX1", "VAX2", "SOX21", "FUBP1","ZNF146","NR6A1","PBX2","NRF1","DLX4","POU4F2","ZFP82","ZFP28")
                
DFRAME_0_f <- DFRAME_0_f[order(factor(DFRAME_0_f$TF, levels = new_levels)),]
DFRAME_0_f$TF <- factor(DFRAME_0_f$TF, levels = new_levels)

#compute wilcoxon and correct (multiple testing) p-values 
anno_df = compare_means(expression ~ group, group.by = "TF", data = DFRAME_0_f, method = "wilcox.test", p.adjust.method = "BH") %>%
  mutate(y_pos = 9) %>%
  filter(group1 == "G2" & group2 == "G4")

#geom_signif function allows us to use adjusted p-values instead of p-values
#plotting genes coding for TF using TCGA data and clinical information - Figure 3A
pdf(paste(wd, "Figure3A.pdf", sep = ""), onefile = T, height = 5, width = 4)
p <- ggboxplot(DFRAME_0_f, x = "group", y = "expression", facet.by = "TF", size = 0.7, ggtheme = theme_bw(),
               color = "group", palette = colores,add="jitter",
               order = c("G2", "G4", "CL"),
               ylab = "log2(FPKM)", xlab = "",
               title = "",
               legend = "right",
               outline = TRUE,
               outlier.shape = NA, notch = FALSE,bxp.errorbar = TRUE) +
  ggsignif::geom_signif(
    data=anno_df, 
    aes(xmin=group1, xmax=group2, annotations=p.adj, y_position=y_pos), 
    manual=TRUE, map_signif_level = TRUE, tip.length = c(0.1, 0.1)) + 
  font("ylab", size = 14, color = "black", face = "bold") +
  font("xy.text", size = 12, color = "black", face = "bold") +
  font("title", size = 12, color = "black", face = "bold.italic") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
p
dev.off()


############################### Transcription factor families (HOCOMOCO DB) ###############################
TF_families <- function(grade){
  TF_families <- read.table("/media/adria/5c3592e6-f482-4fb8-8e84-3c77a7776bd3/BMO_analysis/HUMAN_mono_motifs.tsv",header = TRUE, sep = "\t")
  TF_families$Model <- str_remove(TF_families$Model, "HUMAN.H11MO.")
  TF_families$Model <- str_replace(TF_families$Model, "\\.", "_")
  colnames(TF_families)[1] <- "TF"
  if (grade == "II"){
    TF_families_merged <- merge(final_G2, TF_families, by="TF")
  }
  if (grade == "IV"){
    TF_families_merged <- merge(final_G4, TF_families, by="TF")
  }
  TFs <- as.data.frame(sort(table(TF_families_merged$TF.family)))
  TFs <- TFs[TFs$Freq >0,]
  TFs <- TFs[TFs$Var1 != "",]
  paletteLength = length(unique(TFs$Var1))
  lots_colors <- inferno(paletteLength)
  TFs$Var1 <- factor(TFs$Var1,levels = rev(TFs$Var1) )
  table_percent <- TFs %>%
    mutate(Var1 = Var1,
           cumulative = cumsum(Freq),
           midpoint = cumulative - Freq / 2,
           labels = paste0(round((Freq/ sum(Freq)) * 100, 1), "%"))
  return(list(table_percent, lots_colors))
}

TF_families_G2 <- TF_families("II")
TF_families_G4 <- TF_families("IV")

TF_families_G2[[1]]$group <- "G2"
TF_families_G4[[1]]$group <- "G4"

paletteLength = length(unique(TF_families_G4[[1]]$Var1))
lots_colors <- viridis(paletteLength)
#PIEplot 
G4_fam <- ggplot(TF_families_G4[[1]], aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar(theta = "y", start=0,direction = -1) +
  theme_void() +
  labs(fill = "") +
  guides(fill = guide_legend(reverse=FALSE)) +
  scale_fill_manual(values = lots_colors) +
  geom_text(aes(x = 1.2, y = midpoint , label = labels), color="white",
            fontface = "bold", size = 6) +
  theme(legend.text = element_text(size=12, face = "plain"))

paletteLength = length(unique(TF_families_G2[[1]]$Var1))
lots_colors <- viridis(paletteLength)
G2_fam <- ggplot(TF_families_G2[[1]], aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar(theta = "y", start=0,direction = -1) +
  theme_void() +
  labs(fill = "") +
  guides(fill = guide_legend(reverse=FALSE)) +
  scale_fill_manual(values = lots_colors) +
  geom_text(aes(x = 1.2, y = midpoint , label = labels), color="white",
            fontface = "bold", size = 6) +
  theme(legend.text = element_text(size=14, face = "plain"))

#to plot grade-specific TF families using HOCOMOCO's database - Figure 1E and 1F
pdf(paste(wd, "Figure1E_Figure1F.pdf", sep = ""), onefile = T, height = 10, width = 10)
ggarrange(G4_fam, G2_fam, ncol = 1,nrow = 2, align = "v")
dev.off()

############################### cJUN protein activation (pS73) ###############################
#reading again TCGA data considering only G2 and G4 patients
sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)
data<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$FPKM) #obtaining normalized counts
nm<-sapply(1:length(fd), function(i) names(glioma.tcga[[i]]$mRNAseq[[1]]))
colnames(data)<-unique(substr(unlist(nm),1,12))

ctrl<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$isControl)
CTRL<-which(ctrl=="YES")
dip<-sapply(1:length(fd), function(i) as.vector(unlist(glioma.tcga[[i]]$clinical))[26]) #taking the grade
dip<-str_replace(dip, "NULL", "CTRL")
dip<- as.factor(dip)

#select only G2 and G4
a <-which((dip=="G2") | (dip =="G4"))
data <- data[,a]
dip <- dip[a]
dip <- as.factor(dip)
condition <- dip
condition <- as.vector(condition); condition <- as.factor(condition)

#prefilter based on mean
keep <- rowSums(data) >= 5
data <- data[keep,]
data_matrix <- as.matrix(data)

protein_expression <- read.table("GBMLGG.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt",header = TRUE, sep = "\t")
c_JUN <- protein_expression[protein_expression$Sample.REF == "c-Jun_pS73",]
names(c_JUN) <- substring(names(c_JUN),1,12) #make names shorted and consisten with gliona.tcga object
c_JUN$Sample.REF <- NULL
colnames(c_JUN) <- gsub("\\.", "-", colnames(c_JUN))
c_JUN_t <- as.data.frame(t(as.matrix(c_JUN)))
rownames(c_JUN_t) <- gsub("\\.", "-", rownames(c_JUN_t))

#discard some patients
remove <- c("TCGA-TQ-A8XE-1", "TCGA-TQ-A7RV-1", "TCGA-TQ-A7RK-1","TCGA-19-0957-1","TCGA-06-0190-1", "TCGA-FG-5963-1", "TCGA-06-0171-1",
            "TCGA-FG-A4MT-1","TCGA-06-0125-1","TCGA-19-1389-1")

c_JUN_t <- subset(c_JUN_t, !(rownames(c_JUN_t)) %in% remove)

#obtaining glioma's grade from clinical object
grade_info = list()
counter=0
for (id in rownames(c_JUN_t)){
  counter = counter + 1
  grade_info[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[26] 
}

grade_info[sapply(grade_info, is.null)] <- NA 

grade_info <- unlist(grade_info)
c_JUN_t$grade <- grade_info
colnames(c_JUN_t) <- c("protein_expression", "group")
c_JUN_t <- c_JUN_t[!c_JUN_t$group == "-",]
c_JUN_t$protein_expression <- as.numeric(as.character(c_JUN_t$protein_expression))
c_JUN_t <- na.omit(c_JUN_t)

myplots <- list()
for (i in 1:nrow(JUN)){
  print(i)
  JUN_df <- as.data.frame(data_matrix[JUN$ensembl_gene_id[i],]) #data2 matrix contains FKPM values from the TCGA
  colnames(JUN_df) <- "RPKM"
  JUN_FINAL <- merge(c_JUN_t, JUN_df, by = 0) #by rowNames
  JUN_FINAL$group <- as.factor(JUN_FINAL$group)
  plot <- ggscatter(size = 4, JUN_FINAL, y = 'RPKM', x = 'protein_expression',
                       add = "reg.line", conf.int = FALSE, shape = 19, 
                       cor.coef = TRUE, cor.method = "pearson", 
                       ylab = paste(JUN$hgnc_symbol[i], "mRNA (FPKM)", sep = " "), 
                       color = "group", xlab = "c-Jun_pS73 expression (RPPA)", 
                       palette = c("G4"="firebrick3", "G2"="#1B9E77"), 
                       point = TRUE, cor.coef.size=6)
  stat_cor
  plot <- ggpar(plot, font.xtickslab = c(11, "black"), font.ytickslab = c(11, "black"), font.x = c(12, "bold"), font.y = c(12, "bold"))
  myplots[[i]] <- plot
}

#plotting c-Jun phosphorylation correlation with mRNA expression of c-Jun targets using the TCGA - Figure 3E
pdf(paste(wd, "Figure3E.pdf", sep = ""), onefile = T, height = 5, width = 7)
myplots
dev.off()


###############################  GEPIA2 expression in normal and cancer samples ###############################
GEPIA_JUN <- read.table("GEPIA2_JUN_expression_normal_cancers.txt", sep = "\t", header = F); colnames(GEPIA_JUN) <- c("cancer_type", "tumor", "normal")

#order based on differences between normal and cancer
GEPIA_JUN2 <- GEPIA_JUN; GEPIA_JUN <- melt(GEPIA_JUN, id.vars = c("cancer_type"));GEPIA_JUN2$dif <- GEPIA_JUN2[,3] - GEPIA_JUN2[,2]
GEPIA_JUN2 <- GEPIA_JUN2[,c(1,4)]; GEPIA_JUN <- merge(GEPIA_JUN, GEPIA_JUN2, by="cancer_type")

p <- ggplot(GEPIA_JUN, aes(x=reorder(cancer_type, -dif), y=value, width=0.8, fill = variable)) +
  geom_bar(stat="identity", position = "dodge", width = 0.6) +
  xlab("") +
  ylab("JUN's expression (meadian TPM)") +
  scale_y_continuous(breaks = c(25,50,75,100,125,150,175,200), expand = c(0,0)) +
  scale_fill_manual(values = c("brown2", "cornflowerblue"))
#geom_vline(xintercept = seq(from = 1.5, to = 31, by = 1), linetype="dotdash", color = "black", size=0.8)

#to plot GEPIA2 mRNA JUN expression frequency in PAN cancer - Figure 2A; statistics applied externally
pdf(paste(wd, "Figure2A.pdf", sep = ""), onefile = T, height = 8, width = 14)
p + theme(axis.text.y =element_text(hjust = 1, size = 14, face ="bold", colour = "black"), axis.text.x =element_text(hjust = 1,face ="bold", colour = "black", size = 14, angle = 45), 
          axis.ticks = element_line(colour = "black", size = 1), axis.title.x=element_text(face ="bold", colour = "black", size = 14),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5), axis.title.y =element_text(size = 15, face ="bold", colour = "black"),
          legend.text = element_text(size=12, face = "bold"),panel.background = element_blank())
dev.off()


###############################  Circular plot to show TF-genes in G4 and G2 ###############################
chrOrder <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
              "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

bed <- read.csv("JUN_targets_genes_coordinates_hugo.bed", header = TRUE, sep = "\t") #hg38 coordinates
bed <- bed[,c(2,3,4,1)]
colnames(bed) <- c("chr", "start", "end", "gene_name")
bed[,4] <- as.character(bed[,4])
bed$chr <- factor(bed$chr, chrOrder, ordered=TRUE)
bed <- bed[do.call(order, bed[, c("chr", "start")]), ]

#reading MACS2 peaks
GBM1_TF_peaks <- read.csv("BR20170224_peaks.broadPeak", header = FALSE, sep = "\t")
GBM2_TF_peaks <- read.csv("BR20170331_peaks.broadPeak", header = FALSE, sep = "\t")
LN18_TF_peaks <- read.csv("ATACseq_LN18_peaks.broadPeak",header = FALSE, sep = "\t")
LN229_TF_peaks <- read.csv("ATACseq_LN229_peaks.broadPeak",header = FALSE, sep = "\t")

peaks_ready_circos <- function(peak_file){
  peaks <- peak_file
  peaks <- peaks[c(1,2,3)]
  
  if (colnames(peaks)[1] == "V1"){
    colnames(peaks) <- c("chr", "start", "end")
  }
  return(peaks)
}

GBM1_TF_peaks <- peaks_ready_circos(GBM1_TF_peaks)
GBM2_TF_peaks <- peaks_ready_circos(GBM2_TF_peaks)
LN18_TF_peaks <- peaks_ready_circos(LN18_TF_peaks)
LN229_TF_peaks <- peaks_ready_circos(LN229_TF_peaks)

G4_genes2 <- getBM(G4_genes,  filters = "entrezgene_id", mart = mart, attributes = c("ensembl_gene_id", "entrezgene_id"))
G4_genes2 <- G4_genes2[!duplicated(G4_genes2[ , c("entrezgene_id")]),] #duplicated entrezgenes 
G4_genes2$entrezgene_id <- NULL
G2_genes2 <- getBM(G2_genes,  filters = "entrezgene_id", mart = mart, attributes = c("ensembl_gene_id", "entrezgene_id"))
G2_genes2 <- G2_genes2[!duplicated(G2_genes2[ , c("entrezgene_id")]),] #duplicated entrezgenes 
G2_genes2$entrezgene_id <- NULL

#select TFBS in promoters (common per LN18 and LN229)
LN18_TF_BMO_common_promoters <- as.data.frame(BMO_LN18_LN229_pc_genes)
LN18_TF_BMO_common_promoters <- LN18_TF_BMO_common_promoters[,c(1,10,11)]
LN18_TF_BMO_common_promoters$chr <- paste("chr", LN18_TF_BMO_common_promoters$chr, sep = "")

#select TFBS OUTSIDE promoters (common per LN18 and LN229), 'anti' option used
LN18_TF_BMO_common_outside_promoters <- genome_intersect(BMO_LN18_LN229, all_genes, by=c("chr", "start", "end"), mode = "anti")
LN18_TF_BMO_common_outside_promoters <- as.data.frame(LN18_TF_BMO_common_outside_promoters)
LN18_TF_BMO_common_outside_promoters <- LN18_TF_BMO_common_outside_promoters[,c(1,2,3)]
LN18_TF_BMO_common_outside_promoters$chr <- paste("chr", LN18_TF_BMO_common_outside_promoters$chr, sep = "")
#select TFBS in over expressed G4 genes / G2 genes (UNIQUE!)
LN18_TF_BMO_common_promoters_G4_TFs <- BMO_LN18_TF[BMO_LN18_TF$TF %in% paste(final_G4$TF, "_TFBS",sep = ""),]
LN18_TF_BMO_common_promoters_G4_TFs <- as.data.frame(LN18_TF_BMO_common_promoters_G4_TFs)
#select high-expressed in G4
LN18_TF_BMO_common_promoters_G4_TFs <- LN18_TF_BMO_common_promoters_G4_TFs[LN18_TF_BMO_common_promoters_G4_TFs$ensembl_gene_id %in% G4_genes2$ensembl_gene_id,]
LN18_TF_BMO_common_promoters_G4_TFs <- LN18_TF_BMO_common_promoters_G4_TFs[,c(1,10,11)]
LN18_TF_BMO_common_promoters_G4_TFs$chr <- paste("chr", LN18_TF_BMO_common_promoters_G4_TFs$chr, sep = "")
LN18_TF_BMO_common_promoters_G4_TFs$chr <- factor(LN18_TF_BMO_common_promoters_G4_TFs$chr, chrOrder, ordered=TRUE)
LN18_TF_BMO_common_promoters_G4_TFs <- LN18_TF_BMO_common_promoters_G4_TFs[do.call(order, LN18_TF_BMO_common_promoters_G4_TFs[, c("chr", "start", "end")]), ]
LN18_TF_BMO_common_promoters_G4_TFs <- LN18_TF_BMO_common_promoters_G4_TFs[!duplicated(LN18_TF_BMO_common_promoters_G4_TFs$start),]
LN18_TF_BMO_common_promoters_G2_TFs <- BMO_LN18_TF[BMO_LN18_TF$TF %in% paste(final_G2$TF, "_TFBS",sep = ""),]
LN18_TF_BMO_common_promoters_G2_TFs <- as.data.frame(LN18_TF_BMO_common_promoters_G2_TFs)
#select high-expressed in G2
LN18_TF_BMO_common_promoters_G2_TFs <- LN18_TF_BMO_common_promoters_G2_TFs[LN18_TF_BMO_common_promoters_G2_TFs$ensembl_gene_id %in% G2_genes2$ensembl_gene_id,]
LN18_TF_BMO_common_promoters_G2_TFs <- LN18_TF_BMO_common_promoters_G2_TFs[,c(1,10,11)]
LN18_TF_BMO_common_promoters_G2_TFs$chr <- paste("chr", LN18_TF_BMO_common_promoters_G2_TFs$chr, sep = "")
LN18_TF_BMO_common_promoters_G2_TFs$chr <- factor(LN18_TF_BMO_common_promoters_G2_TFs$chr, chrOrder, ordered=TRUE)
LN18_TF_BMO_common_promoters_G2_TFs <- LN18_TF_BMO_common_promoters_G2_TFs[do.call(order, LN18_TF_BMO_common_promoters_G2_TFs[, c("chr", "start")]), ]
LN18_TF_BMO_common_promoters_G2_TFs <- LN18_TF_BMO_common_promoters_G2_TFs[!duplicated(LN18_TF_BMO_common_promoters_G2_TFs$start),]

#prepare legend for the circos
lgd1 = Legend(labels = c("LN18", "LN229", "GBM1", "GBM2"), labels_gp = gpar(fontsize = 10), border = "black", legend_gp = gpar(fill = c("lightgoldenrod1", "deepskyblue", "brown1", "mediumorchid1")), title = "ATACseq peaks density")
lgd2 = Legend(labels = c("TFBS on promoters (n=101962)", "TFBS outside promoters (n=23463)"), labels_gp = gpar(fontsize = 10),border = "black",legend_gp = gpar(fill = c("aquamarine3", "orange3")), title = "Common TF binding sites")
lgd3 = Legend(labels = c("TFBS on GII- genes promoters (n=24)", "TFBS on GIV- genes promoters (n=240)"), labels_gp = gpar(fontsize = 10),border = "black",legend_gp = gpar(fill = c("#1B9E77", "firebrick3")), title = "Grade-specific TF binding sites")
lgd4 = Legend(at = c("c-Jun's binding"), type = "lines", 
              legend_gp = gpar(col = "firebrick1", lwd = 2))
pd = packLegend(lgd1, lgd2, lgd3, lgd4, direction = "horizontal")

#a) Circos plotting ATACseq peaks in both GBM cell lines and GBM specimens, as well as TFBS - Figure 2C
pdf(file = paste(wd,"Figure2C.pdf"), width = 12, height = 10)
#parameters and initialize the plot
par(mar = c(1, 1, 1, 1)) 
circos.par("start.degree" = 90)
circos.par("track.height" = 0.1)
circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))
circos.initializeWithIdeogram(species='hg38', chromosome.index = paste0("chr", c(1:22, "X", "Y")))
set.seed(123)
circos.initializeWithIdeogram(plotType = NULL)
draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.10, bg.border = NA)

circos_colors <- c(adjustcolor("lightgoldenrod1", alpha.f = 0.8), adjustcolor("deepskyblue", alpha.f = 0.3))
circos_colors2 <- c(adjustcolor("aquamarine3", alpha.f = 0.4), adjustcolor("peru", alpha.f = 0.8))
circos_colors3 <- c(adjustcolor("brown1", alpha.f = 0.4), adjustcolor("mediumorchid1", alpha.f = 0.8))
circos_colors_test <- c(adjustcolor("lightgoldenrod1", alpha.f = 0.8), adjustcolor("deepskyblue", alpha.f = 0.3),
                        adjustcolor("brown1", alpha.f = 0.4), adjustcolor("mediumorchid1", alpha.f = 0.8))

#plot the broad peaks from ATACseq (GBM cell lines)
circos.genomicDensity(list(LN18_TF_peaks, LN229_TF_peaks), baseline = 0, col = circos_colors, track.height = 0.08,border = TRUE, lwd = 0.8, type = "l")

#plot the broad peaks from ATACseq (GBM tumours)
circos.genomicDensity(list(GBM1_TF_peaks, GBM2_TF_peaks), baseline = 0, col = circos_colors3, track.height = 0.08,border = TRUE, lwd = 0.8, type = "l")

#plot the TFBS in promoter and outside regions that are common in both cell lines
circos.genomicRainfall(list(LN18_TF_BMO_common_promoters, LN18_TF_BMO_common_outside_promoters), pch = 16, cex = 0.2, col = circos_colors2)

#plot the TF specific in G2 and G4 genes
circos.genomicDensity(list(LN18_TF_BMO_common_promoters_G2_TFs, LN18_TF_BMO_common_promoters_G4_TFs), 
                      baseline = 0, col = c("#1B9E77", "firebrick3"), track.height = 0.08,border = FALSE, lwd = 0.8, type = "l")

#links
bed_JUN <- bed
bed_JUN$chr <- rep("chr1", 16)
bed_JUN$start <- as.integer(rep("58776845", 16))
bed_JUN$end <- as.integer(rep("58784048", 16))
bed_JUN$gene_name <- "cJUN"

circos.genomicLink(bed_JUN, bed, col = adjustcolor("firebrick3", alpha.f = 1))

#labbels
circos.genomicLabels(bed, labels.column = 4, side = "inside",cex = 0.7,connection_height = convert_height(4, "mm"), 
                     niceFacing = TRUE, line_lty = "solid", padding = 1)
dev.off()

#Stepniak et al. glioma enhancers (H3K27ac)
H3K27ac_GBM_5_samples <- read.csv("H3K27ac_GBM_peaks_in_5_samples.narrowPeak",header = FALSE, sep = "\t")
H3K27ac_DA_5_samples <- read.csv("H3K27ac_DA_peaks_in_5_samples.narrowPeak",header = FALSE, sep = "\t")
#all enhancers
enhancers <- read.table("common_putative_enhancers.bed",header = FALSE, sep = "\t")
colnames(enhancers) <- c("chr", "start", "end"); enhancers[,4] <-NULL

#enhancers intersected with BMO TFBS
enhancers_TF <- as.data.frame(BMO_LN18_LN229_enhancers)
enhancers_TF <- enhancers_TF[,c(1,5,6)]
enhancers_TF$chr <- paste("chr",enhancers_TF$chr, sep = "")
enhancers_TF <- makeGRangesFromDataFrame(df = enhancers_TF, seqnames.field = c("chr"), 
                                         keep.extra.columns = FALSE, 
                                         start.field = "start")

#which peaks in promoters? Discard Symfonia peaks within promoters
peaks <- as.data.frame(enhancers_peakAnno@anno)
peaks_promoters <- subset(peaks, annotation == "Promoter")
peaks_no_promoters <- genome_intersect(peaks, peaks_promoters, by=c("seqnames", "start", "end"), mode="anti")
peaks_no_promoters <- makeGRangesFromDataFrame(df = as.data.frame(peaks_no_promoters), seqnames.field = c("seqnames"), 
                                               keep.extra.columns = FALSE, 
                                               start.field = "start")
enhancers_peakAnno <- annotatePeak(peaks_no_promoters, tssRegion=c(-500, 500),TxDb=txdb, annoDb="org.Hs.eg.db")

#discarding peaks thay are within gene promoters
enhancers_peakAnno <- annotatePeak(enhancers_TF, tssRegion=c(-500, 500),TxDb=txdb, annoDb="org.Hs.eg.db")

#plotting H3K27ac peaks distribution - Figure 4B
pdf(paste(wd, "Figure4B.pdf", sep = ""), height = 3, width = 8)
plot_anno <- plotAnnoBar(enhancers_peakAnno, col = COLS)
plot_anno + theme(axis.text.y = element_text(hjust=0.5, size = 12,color = "black",face = "bold",angle = 90),
                  axis.text.x=element_text(size = 12,color = "black",face = "bold"),axis.title =element_text(size = 12,color = "black",face = "bold"),
                  legend.title = element_text(face = "bold", size= 0), legend.text = element_text(face = "bold", size= 12))
dev.off()

enhancers_peakAnno <- as.data.frame(enhancers_peakAnno)
colnames(enhancers_peakAnno)[1] <- "chr"
enhancers_peakAnno$annotation <- gsub(enhancers_peakAnno$annotation,pattern = "Intron.*", replacement = "Intron")
enhancers_peakAnno$annotation <- gsub(enhancers_peakAnno$annotation,pattern = "Exon.*", replacement = "Exon")
col_anno_gene = structure(COLS, names = unique(enhancers_peakAnno$annotation))

#enhancers intersected with BMO JUN TFBS
JUN_in_enhancers <- as.data.frame(BMO_LN18_LN229_enhancers[ grep("JUN_0_A*", BMO_LN18_LN229_enhancers$TF),])
JUN_in_enhancers$chr <- paste("chr", JUN_in_enhancers$chr, sep = "") 

#prepare legend
lgd1 = Legend(labels = c("GBM", "DA"), labels_gp = gpar(fontsize = 10), border = "black", legend_gp = gpar(fill = c("brown1", "gold")), title = "ChIPseq H3K27ac peaks density")
lgd2 = Legend(labels = c("H3K27ac peaks (n=10673)"), labels_gp = gpar(fontsize = 10),border = "black",legend_gp = gpar(fill = c("darkmagenta")), title = "Putative glioma enhancers")
lgd3 = Legend(labels = c("TFBS on glioma enhancers (n=7571)", "cJun's binding sites on glioma enhancers (n=94)"), labels_gp = gpar(fontsize = 10),border = "black",legend_gp = gpar(fill = c("cyan3", "coral3")), title = "TF binding sites")
lgd4 = Legend(at = c("c-Jun's binding"), type = "lines", 
              legend_gp = gpar(col = "firebrick1", lwd = 2))

#lgd3 = Legend(title = "Gene anno", at = names(col_anno_gene), legend_gp = gpar(fill = col_anno_gene))
pd = packLegend(lgd1, lgd2, lgd3, lgd4, direction = "horizontal")

#b) Circos plotting H3K27ac peaks (Stepniak et al.,) in both GBMs and DA patients - Figure 4A
pdf(file = paste(wd,"Figure4A.pdf"), width = 13, height = 13)
par(mar = c(1, 1, 1, 1)) 
circos.par("start.degree" = 90)
circos.par("track.height" = 0.1)
circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))
circos.initializeWithIdeogram(species='hg38', chromosome.index = paste0("chr", c(1:22, "X", "Y")))
set.seed(123)
circos.initializeWithIdeogram(plotType = NULL)
draw(pd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.10, bg.border = NA)

#plot top H3K27ac in GBM and DA
circos.genomicDensity(list(H3K27ac_GBM_5_samples, H3K27ac_DA_5_samples), baseline = 0, col = c("brown1", "orange"), track.height = 0.08,border = TRUE, lwd = 0.8, type = "l")

#all enhancers
circos.genomicRainfall(enhancers, pch = 16, cex = 0.3, col = "darkmagenta")

#JUN and FOS related bindings in contrast with total TFBS in enhancers
circos.genomicRainfall(list(enhancers_peakAnno, JUN_in_enhancers[,c(1,5,6)]), pch = c(16,17), cex = c(0.3, 0.4), col = c("cyan3", "coral3"))

bed_JUN <- data.frame(chr = "chr1", start = "58776845", end = "58784048")
bed_JUN <- as.data.frame(lapply(bed_JUN, rep, 94))

circos.genomicLink(bed_JUN, JUN_in_enhancers[,c(1,5,6)], col = adjustcolor("firebrick3", alpha.f = 0.4))
dev.off()


############################### Computing ATACseq signal within TSS regions ###############################
files <- list(LN18 = "ATACseq_LN18_peaks.broadPeak", 
              LN229 = "ATACseq_LN229_peaks.broadPeak",
              GBM1 = "BR20170224_peaks.broadPeak",
              GBM2 = "BR20170331_peaks.broadPeak")

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

#ATACseq broad peaks - Figure 1B
pdf(file = paste(wd,"Figure1B.pdf"), width = 12, height = 10)
peakHeatmap(files, TxDb=txdb, upstream=3000, downstream=3000, color=c("lightgoldenrod1", "deepskyblue", "brown1", "darkorchid1"))
dev.off()


########################## Enhancers with JUNs TFBS and HiC #############################
HiC_targets <- read.table("enhancers_hg38_TSS_Won_et_al_HiC.txt",header = TRUE)
colnames(HiC_targets) <- c("chr", "start", "end", "ensembl_transcript_id", "ensembl_gene_id", "-log10(qval)", "confirmed_both_ways")

#focusing on enhancers harboring a JUN motif
JUN_in_enhancers2 <- JUN_in_enhancers
JUN_in_enhancers2 <- JUN_in_enhancers2[,c(1,5,6,2,3,4)]
HiC_targets_JUN_enhancers <- genome_intersect(HiC_targets, JUN_in_enhancers2, mode = "both")
HiC_targets_JUN_enhancers <- as.data.frame(HiC_targets_JUN_enhancers)
HiC_targets_JUN_enhancers <- subset(HiC_targets_JUN_enhancers, ensembl_gene_id != "-")
HiC_targets_JUN_enhancers <- HiC_targets_JUN_enhancers[!duplicated(HiC_targets_JUN_enhancers$ensembl_gene_id),]
HiC_targets_JUN_enhancers_genes <- getBM(HiC_targets_JUN_enhancers$ensembl_gene_id, filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"))
HiC_targets_JUN_enhancers_genes <- HiC_targets_JUN_enhancers_genes[complete.cases(HiC_targets_JUN_enhancers_genes$hgnc_symbol), ]
HiC_targets_JUN_enhancers_genes <- merge(HiC_targets_JUN_enhancers, HiC_targets_JUN_enhancers_genes, by = "ensembl_gene_id")
write.table(HiC_targets_JUN_enhancers_genes, "enhancers_common_TF_LN18_LN229_JUN_HiC_genes.bed",quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

#discard enhancers in promoters
enhancers_top_TF_CL <- read.table("common_putative_enhancers_common_TF_LN18_LN229.bed", header = FALSE)
enhancers_top_TF_tumour <- read.table("common_putative_enhancers_common_TF_GBM1_GBM2.bed", header = FALSE)
enhancers_top_TF <- cbind(enhancers_top_TF_CL, enhancers_top_TF_tumour[,c(4)])

colnames(enhancers_top_TF) <- c("seqnames", "Start", "End", "number_TFBS_CL", "number_TFBS_tumour")
enhancers_top_TF <- makeGRangesFromDataFrame(df = enhancers_top_TF, seqnames.field = c("seqnames"), 
                                             keep.extra.columns = TRUE, 
                                             start.field = "Start")
enhancers_top_TF <- annotatePeak(enhancers_top_TF, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
enhancers_top_TF <- as.data.frame(enhancers_top_TF@anno)
enhancers_top_TF <- enhancers_top_TF[- grep("Promoter*", enhancers_top_TF$annotation),]
enhancers_top_TF <- enhancers_top_TF[,c(1,2,3,6,7)]
enhancers_top_TF <- enhancers_top_TF[order(enhancers_top_TF$number_TFBS_CL, decreasing = TRUE),]
#write.table(enhancers_top_TF, "common_putative_enhancers_common_TF_LN18_LN229_GBM1_GBM2_no_promoters_fix.bed",quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

enhancers_top_TF2 <- enhancers_top_TF
enhancers_top_TF2 <- enhancers_top_TF2[enhancers_top_TF2$number_TFBS_CL >=30,] #selecting top enhancers with more motifs
enhancers_top_TF2$region <- paste(enhancers_top_TF2$seqnames, enhancers_top_TF2$start, enhancers_top_TF2$end, sep = "-")

#plot total number of predicted TFBS in glioma enhancers, using n=>30 as an arbitrary threshold - Supplementary Figure 8B
pdf(file = paste(wd,"Supplementary_Figure8B.pdf"), width = 10, height = 8)
q <- ggplot(data=enhancers_top_TF2, aes(x=reorder(region, -number_TFBS_CL) , y=number_TFBS_CL, width=0.8, fill=region)) +
  geom_bar(stat="identity", position="stack") +
  xlab("") +
  ylab("TF binding sites in glioma enhancers") +
  theme(legend.position="none") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept=30, linetype = "dashed", color = "black", size=1) +
  theme(title = element_text(size = 12, face = "bold"),axis.text.y =element_text(size = 12, face ="bold", colour = "black"), axis.text.x =element_text(hjust = 1,face ="bold", angle = 90,colour = "black", size = 12), 
        axis.ticks = element_line(colour = "black", size = 1), axis.title.x=element_text(face ="bold", colour = "black", size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5), plot.margin = margin(1, 1, 1, 1, "cm"),panel.background = element_blank())
q
dev.off()

######################## H3K27ac signal differences between glioma grades (Stepniak et al.,) ################################
peaks_H3K27ac <- dir(wd, pattern = "*peaks.narrowPeak", full.names = TRUE)
myPeaks <- lapply(peaks_H3K27ac, ChIPQC:::GetGRanges, simple = TRUE)

names(myPeaks) <- c("DA01", "DA02", "DA03", "DA04", "DA05", "DA06", "DA07", "GB01", "GB02", "GB03", "GB04", "GB05", "GB06", "GB07", "GB08","GB09", "GB10")
Group <- factor(c(rep("DA", 7), rep("GBM", 10)))
myGRangesList <- GRangesList(myPeaks)   
reduced <- GenomicRanges::reduce(unlist(myGRangesList))
consensusIDs <- paste0("consensus_", seq(1, length(reduced)))
mcols(reduced) <- do.call(cbind, lapply(myGRangesList, function(x) (reduced %over% x) + 0))
reducedConsensus <- reduced
mcols(reducedConsensus) <- cbind(as.data.frame(mcols(reducedConsensus)), consensusIDs)
consensusIDs <- paste0("consensus_", seq(1, length(reducedConsensus)))
consensusToCount <- reducedConsensus 

#discard peaks within ensembl blacklist
blklist <- import.bed("GRCh38_unified_blacklist.bed")
consensusToCount <- consensusToCount[!consensusToCount %over% blklist & !seqnames(consensusToCount) %in% "chrM"]
consensusToCount

#peak overlaps to get an overall view of correspondance between peak calls (not published)
myPlot <- as.data.frame(elementMetadata(consensusToCount)) %>% dplyr::select(-consensusIDs) %>% 
  as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = Group)) + geom_point(size = 5) + scale_color_manual(values = c("firebrick3","#1B9E77")) +
  theme(axis.text = element_text(size = 14,color = "black",face = "bold"), 
        legend.title = element_text(size = 0), 
        axis.title = element_text(size = 14,color = "black",face = "bold"),
        legend.text = element_text(size=12))
myPlot 

#identifying changes of ChiP-seq signal within peaks will allow us to better understand changes accurately
occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% rowSums
table(occurrences) %>% rev %>% cumsum #unique and common peaks
consensusToCount <- consensusToCount[occurrences >= 2, ] #a peaks identified at least on 2 samples

#set regions to count paired reads leading to peaks
bamsToCount <- dir(wd, full.names = TRUE, pattern = "*uniq.bam$")
regionsToCount <- data.frame(GeneID = paste("ID", seqnames(consensusToCount), start(consensusToCount), end(consensusToCount), sep = "_"), Chr = seqnames(consensusToCount), 
                             Start = start(consensusToCount), End = end(consensusToCount), Strand = strand(consensusToCount))

#counting reads within peaks 
fcResults <- featureCounts(bamsToCount, annot.ext = regionsToCount, isPairedEnd = FALSE, countMultiMappingReads = FALSE,nthreads = 8)
myCounts <- fcResults$counts
colnames(myCounts) <- c("DA01", "DA02", "DA03", "DA04", "DA05", "DA06", "DA07", "GB01", "GB02", "GB03", "GB04", "GB05", "GB06", "GB07", "GB08","GB09", "GB10")

#DESeq2 for differential ChIP-seq
metaData <- data.frame(Group, row.names = colnames(myCounts))
chipDDS <- DESeqDataSetFromMatrix(myCounts, metaData, design = ~Group, rowRanges = consensusToCount)
#atacDDS_test <- estimateSizeFactors(atacDDS) #to estimate size factors (library size) for bamCoverage
chipDDS <- DESeq(chipDDS)
chip_vst <- vst(chipDDS)
z <- plotPCA(chip_vst, intgroup = "Group", ntop = nrow(chip_vst))
nudge <- position_nudge(y = 1)
#not published
z  + scale_color_manual(values = c("firebrick3","#1B9E77")) +
  theme(legend.title = element_text(size = 0), 
        axis.title = element_text(size = 14,color = "black",face = "bold"),
        axis.text = element_text(size = 14,color = "black",face = "bold"),
        legend.text = element_text(size=12))

results <- results(chipDDS, c("Group", "GBM", "DA"), format = "GRanges")
results <- results[order(results$padj)]
results <- results[(!is.na(results) & results$padj < 0.05),]
sum(results$padj < 0.05, na.rm=TRUE)
#OPTIONAL results <- results[(!is.na(results) & results$pvalue < 0.01 & results$log2FoldChange < 0),]
results <- as.data.frame(results)
JUN_in_enhancers3 <- JUN_in_enhancers2
colnames(JUN_in_enhancers3)[1] <- "seqnames"

#intersect dif H3K27ac regions with super-enhancers
H3K27ac_signal_enhancers <- genome_intersect(results, enhancers_top_TF3, by=c("seqnames", "start", "end"), mode = "both") 
H3K27ac_signal_enhancers <- as.data.frame(H3K27ac_signal_enhancers)
H3K27ac_signal_enhancers <- subset(H3K27ac_signal_enhancers, number_TFBS > 0)

#intersect diffferent H3K27ac regions with c-JUN enhancers - Supplementary Table 3
H3K27ac_signal_JUN <- genome_intersect(results, JUN_in_enhancers3, by=c("seqnames", "start", "end"), mode = "both") 
H3K27ac_signal_JUN <- as.data.frame(H3K27ac_signal_JUN)
H3K27ac_signal_JUN <- H3K27ac_signal_JUN[,c(1,13,14,2,3,4,5,6,7,8,9,10,11,12)]
write.table(H3K27ac_signal_JUN, "Supplementary_Table3.bed",quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)



####### END #######





