library(GenomicAlignments)
library(GenomicFeatures)
library(Rsamtools)
library("BiocParallel")
library("DESeq2")
library("AnnotationDbi")
library("Mus.musculus")
library(ggrepel)
library(reshape2)
library(fgsea)
library(msigdbr)


setwd("~/MGN/BAM_mergeReference/")
set.seed(seed = 10044) # Set seed for reproducible results in permutation tests

# Load aligned BAM files --------------------------------------------------

filenames <- c("M1-IOL-2dpi_S9Aligned.sortedByCoord.out.bam",
               "M2-IOL-2dpi_S10Aligned.sortedByCoord.out.bam",
               "M3-IOL-2-dpi_S11Aligned.sortedByCoord.out.bam",
               "M4-IOL-2dpi_S12Aligned.sortedByCoord.out.bam",
               "M5-HI-control-2dpi_S13Aligned.sortedByCoord.out.bam",
               "M6-HI-control-2dpi_S14Aligned.sortedByCoord.out.bam",
               "M7-IOL-5dpi_S15Aligned.sortedByCoord.out.bam",
               "M8-IOL-5dpi_S16Aligned.sortedByCoord.out.bam",
               "M9-IOL-5dpi_S17Aligned.sortedByCoord.out.bam",
               "M10-IOL-5dpi_S18Aligned.sortedByCoord.out.bam",
               "M11-HI-control-5dpi_S19Aligned.sortedByCoord.out.bam",
               "M12-HI-control-5dpi_S20Aligned.sortedByCoord.out.bam"
)
file.exists(filenames)

bamfiles <- BamFileList(filenames, yieldSize=2000000)

gtffile <- file.path("~/MGN/Genomes/mm10/gencode.vM25.annotation.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf",metadata = "symbol")

ebg <- transcriptsBy(txdb, by="gene")


#  Add nonstructural polyprotein coding sequence for annotation of the genome --------
ebg$CHIKVgp1 = GRanges(seqnames = "AM258994.1",
                       ranges = IRanges(start = 77,end = 7501),
                       strand = "*",
                       tx_id = 956309,
                       tx_name = "CHIKVgp1")

# Add structural polyprotein coding sequence for annotation of combined genome --------
ebg$CHIKVgp2 = GRanges(seqnames = "AM258994.1",
                       ranges = IRanges(start = 7567,end = 11313),
                       strand = "*",
                       tx_id = 956308,
                       tx_name = "CHIKVgp2")

register(MulticoreParam(18))

se <- summarizeOverlaps(features=ebg,
                        reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE)

#save(se, file = "~/MGN/ProcessedData/se.transcripts.merge.genomes.Rdata")
# Checkpoint - load SE data -----------------------------------------------
load(file = "~/MGN/ProcessedData/se.transcripts.merge.genomes.Rdata")

se
colnames(se)
se$treatment = factor(x = c(rep("IOL",4),rep("HI",2),rep("IOL",4), rep("HI",2)))
se$days = factor(x = c(rep("Day2",6),rep("Day5",6)))
se$combined = factor(x = c(rep("IOL_Day2",4),rep("HI",2),rep("IOL_Day5",4),rep("HI",2)))
se$treatment_day = factor(x = c(rep("IOL_Day2",4),rep("HI_Day2",2),rep("IOL_Day5",4),rep("HI_Day5",2)))

# To check CHIKV genomes -------------------------------------------------
df = as.data.frame(t(assay(se)[grep(rownames(se), pattern = "CHIK", value = T),]))
total.reads = colSums(assay(se))
df$CHIKVgp1.norm = df$CHIKVgp1/total.reads*1000000
df$CHIKVgp2.norm = df$CHIKVgp2/total.reads*1000000

df$treatment = se$treatment
df$days = se$days
df$combined  = paste0(df$treatment,"_",df$days)
df$combined[grep(df$combined, pattern = "HI")] = "HI" # Combine controls for viral reads count

m = melt(df)
m = m[m$variable %in% c("CHIKVgp1.norm","CHIKVgp2.norm"),]
m$variable = as.character(m$variable)
m$variable[m$variable == "CHIKVgp1.norm"] = "Non-structural protein"
m$variable[m$variable == "CHIKVgp2.norm"] = "Structural protein"

med = lapply(unique(as.character(m$combined)), function(x){
  m1 = median(m$value[m$variable == "Non-structural protein" & m$combined == x])
  m2 = median(m$value[m$variable == "Structural protein" & m$combined == x])
  c(m1,m2)
})
med = as.data.frame(do.call(rbind,med))
rownames(med) = unique(as.character(m$combined))
colnames(med) = c("Non-structural protein", "Structural protein")
med$delta = med$`Structural protein` - med$`Non-structural protein`

ggplot(m, aes(x = combined, y = value, fill = variable)) +
  geom_boxplot() +
  geom_point(pch=21, position = position_dodge(width = 0.75)) +
  labs(y = "Expression (CPM)", x = "Treatment", fill = "Genomic region",caption = paste0("dDay2 : ", round(med$delta[1], digits = 3), "- dDay5 : ", round(med$delta[3], digits = 3))) +
  scale_fill_manual(values = c("blue","firebrick2")) +
  theme_classic()

wilcox.test(m$value[m$combined == "HI" & m$variable == "Non-structural protein"],
            m$value[m$combined == "HI" & m$variable == "Structural protein"])

wilcox.test(m$value[m$combined == "IOL_Day2" & m$variable == "Non-structural protein"],
            m$value[m$combined == "IOL_Day2" & m$variable == "Structural protein"])

wilcox.test(m$value[m$combined == "IOL_Day5" & m$variable == "Non-structural protein"],
            m$value[m$combined == "IOL_Day5" & m$variable == "Structural protein"])


# Check read distribution per gene ----------------------------------------
hist(log10(rowSums(assay(se))), breaks = 100, main = "Distribution of counts per gene")
abline(v = log10(10), lty=2, col = "red")

# Run DEseq2 and check PCA distributions ----------------------------------
# DEseq2
dds <- DESeqDataSet(se, design = ~ combined)
dds$combined <- relevel(dds$combined, ref = "HI")
vsd <- vst(dds)

pca.dat <- plotPCA(vsd, "treatment_day", returnData = T)
ggplot(pca.dat, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("lightgrey","darkgrey","darkblue","darkred"))+
  theme_classic()

dds <- DESeq(dds)
res_d2 <- results(dds, contrast=c("combined","IOL_Day2","HI")) # Specify numerator first and then denominator for FC
res_d5 <- results(dds, contrast=c("combined","IOL_Day5","HI"))

summary(res_d2)
summary(res_d5)

# Add gene symbols --------------------------------------------------------
res_d2$symbol <- mapIds(Mus.musculus,
                        keys=unlist(lapply(strsplit(row.names(res_d2),"[.]"), function(x) x[1])),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

res_d5$symbol <- mapIds(Mus.musculus,
                        keys=unlist(lapply(strsplit(row.names(res_d5),"[.]"), function(x) x[1])),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

res_d2$padj = p.adjust(res_d2$pvalue, method = "BH")
res_d5$padj = p.adjust(res_d5$pvalue, method = "BH")

# Genes removed due to outlier by Cook distance ---------------------------
# Day 2:
out.genes.symbol = as.data.frame(res_d2)[as.data.frame(res_d2)$baseMean !=0 & is.na(as.data.frame(res_d2)$padj),]$symbol
out.genes.ensembl = rownames(as.data.frame(res_d2))[as.data.frame(res_d2)$baseMean !=0 & is.na(as.data.frame(res_d2)$padj)]
names(out.genes.ensembl) = out.genes.symbol
out.genes.ensembl = out.genes.ensembl[!is.na(names(out.genes.ensembl))]
for(i in names(out.genes.ensembl)){
  p = barplot(assay(se)[out.genes.ensembl[i],], main = paste0(i,"-",out.genes.ensembl[i]))
  print(p)
}
d2.out = out.genes.ensembl
# Day 5:
out.genes.symbol = as.data.frame(res_d5)[as.data.frame(res_d5)$baseMean !=0 & is.na(as.data.frame(res_d5)$padj),]$symbol
out.genes.ensembl = rownames(as.data.frame(res_d5))[as.data.frame(res_d5)$baseMean !=0 & is.na(as.data.frame(res_d5)$padj)]
names(out.genes.ensembl) = out.genes.symbol
out.genes.ensembl = out.genes.ensembl[!is.na(names(out.genes.ensembl))]
for(i in names(out.genes.ensembl)){
  p = barplot(assay(se)[out.genes.ensembl[i],], main = paste0(i,"-",out.genes.ensembl[i]))
  print(p)
}
d5.out = out.genes.ensembl
setdiff(d5.out, d2.out)

genes.to.remove = out.genes.ensembl[names(out.genes.ensembl) != "Irf7"]
genes.to.remove

# Remove genes with no baseMean expression
res_d2 = res_d2[res_d2$baseMean != 0,]
res_d5 = res_d5[res_d5$baseMean != 0,]

# Remove genes with no symbol annotation
res_d2 = res_d2[!is.na(res_d2$symbol),]
res_d5 = res_d5[!is.na(res_d5$symbol),]

# Remove Cook outliers (except Irf7)
res_d2 = res_d2[!(res_d2$symbol %in% names(genes.to.remove)),]
res_d5 = res_d5[!(res_d5$symbol %in% names(genes.to.remove)),]

# Write DE lists
# write.table(x = as.data.frame(res_d2)[order(res_d2$log2FoldChange, decreasing = T),],file = "~/MGN/Results/DE_Day2_vs_HI_Day2_mergeGenomes_FINAL.txt", quote = F, row.names = T, col.names = T)
# write.table(x = as.data.frame(res_d5)[order(res_d5$log2FoldChange, decreasing = T),],file = "~/MGN/Results/DE_Day5_vs_HI_Day5_mergeGenomes_FINAL.txt", quote = F, row.names = T, col.names = T)


# Plot volcano ----------------------------------------------------------
# For Day2
dat = as.data.frame(res_d2)
dat = dat[complete.cases(dat),]
dat = dat[order(dat$padj, decreasing = F),]
l = rep(NA, nrow(dat))
top10 = dat$symbol[1:10]
dn10 = dat[dat$log2FoldChange < -1 & dat$padj < 0.15,]
dn10 = dn10[order(dn10$log2FoldChange, decreasing = F),]$symbol[1:10]
cond = dat$symbol %in% c(top10,"Fgf21","Fap","Gli1","Stac3","Hmcn2","Til2","Adam33","Ccl24","Sbk2","Sbk3")

#cond = -log10(dat$padj) > -log10(0.15) & abs(dat$log2FoldChange) > 1 | dat$log2FoldChange < -1 & dat$padj < 0.15
l[cond] = dat$symbol[cond]

p1 = ggplot(dat, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(color = "lightgrey") +
  geom_point(data = dat[dat$log2FoldChange < -1 & dat$padj < 0.15,], color = "darkblue")+
  geom_point(data = dat[dat$log2FoldChange > 1 & dat$padj < 0.15,], color = "firebrick2")+
  theme_classic() +
  labs(x= "log2(FC)", y = "-log10(FDR)", title = "Day 2 post-infection")+
  geom_hline(yintercept = -log10(0.15), lty=2) +
  geom_vline(xintercept = c(-1,1), lty=2) +
  scale_color_manual(values = c("lightgrey","firebrick2")) +
  geom_text_repel(label = l, nudge_y = 10, size = 3, max.overlaps = 50) +
  xlim(-max(abs(dat$log2FoldChange)), max(abs(dat$log2FoldChange)))
print(p1)

# pdf(file = "~/MGN/Results/Volcano_Day2_final.pdf", width = 3.13, height = 3.66)
# print(p1)
# dev.off()

# For Day5
dat = as.data.frame(res_d5)
dat = dat[complete.cases(dat),]
dat = dat[order(dat$padj, decreasing = F),]

l = rep(NA, nrow(dat))
cond = -log10(dat$padj) > 10 & abs(dat$log2FoldChange) > 2 | dat$log2FoldChange < -5 & dat$padj < 0.05

top10 = dat$symbol[1:10]
dn10 = dat[dat$log2FoldChange < -1 & dat$padj < 0.15,]
dn10 = dn10[order(dn10$log2FoldChange, decreasing = F),]$symbol[1:10]
cond = dat$symbol %in% c(top10,"Fgf21","Fap","Gli1","Stac3","Hmcn2","Til2","Adam33","Ccl24","Sbk2","Sbk3")

l[cond] = dat$symbol[cond]
dat$padj[-log10(dat$padj) > 50] = 1e-50
dat$log2FoldChange[dat$log2FoldChange < -10] = -10

p2 = ggplot(dat, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(color = "lightgrey") +
  geom_point(data = dat[dat$log2FoldChange < -1 & dat$padj < 0.15,], color = "darkblue")+
  geom_point(data = dat[dat$log2FoldChange > 1 & dat$padj < 0.15,], color = "firebrick2")+
  theme_classic() +
  geom_hline(yintercept = -log10(0.15), lty=2) +
  geom_vline(xintercept = c(-1,1), lty=2) +
  geom_text_repel(label = l, nudge_y = 10, size = 3, max.overlaps = 50) +
  xlim(-max(abs(dat$log2FoldChange)), max(abs(dat$log2FoldChange)))
print(p2)
#ggsave(p2, file = "~/MGN/Results/Volcano_Day5_final.pdf", height = 4.7,width = 4.1)

# Create expression matrix for heatmap
exp = counts(dds,normalized=TRUE)[order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE),]
exp = exp[rownames(res_d2),]
rownames(exp) = res_d2$symbol
colnames(exp) = c(paste0("IOL-2dpi_",1:4),paste0("HI_2dpi_",1:2),paste0("IOL-5dpi_",1:4), paste0("HI_5dpi_",1:2))

mat <- exp[,c("HI_2dpi_1","HI_2dpi_2","HI_5dpi_1","HI_5dpi_2","IOL-2dpi_1","IOL-2dpi_2","IOL-2dpi_3","IOL-2dpi_4","IOL-5dpi_1","IOL-5dpi_2","IOL-5dpi_3","IOL-5dpi_4")]

genes = unique(c(res_d2$symbol[res_d2$padj < 0.01 & abs(res_d2$log2FoldChange) > 1],res_d5$symbol[res_d5$padj < 0.01 & abs(res_d5$log2FoldChange) > 1] ))
genes= genes[complete.cases(genes)]
mat = mat[genes,]

idx = t(apply(mat,1,scale))
idx = !(idx[,3] > 1.5)
mat = mat[idx,]

paletteLength <- 20
myColor <- colorRampPalette(c("darkblue", "white", "darkred"))(paletteLength)
myBreaks = c(seq(-2,2,length.out = 20))
annotation.col = data.frame(row.names = colnames(mat), Treatment = c(rep("Ctrl",4), rep("CHIKV 2dpi",4), rep("CHIKV 5dpi",4)))

#pdf("~/MGN/Results/All_samples_mergedGenome_Heatmap.pdf", width = 4, height = 8)
pheatmap::pheatmap(mat, cluster_cols = F, show_colnames = T, scale = "row", show_rownames = F, color = myColor, breaks = myBreaks, annotation_col = annotation.col)
#dev.off()

# Pathway analysis --------------------------------------------------------
# HALLMARK -  Day2  -------------------------------------------------------
msigdbr_show_species()
m_df = msigdbr(species = "Mus musculus", category = "H") # Pick the C2 hallmark gene pathways from mDBsig
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

geneList = as.data.frame(res_d2[,c("log2FoldChange","padj","symbol")])
geneList = geneList[!duplicated(geneList$symbol),]
geneList = geneList[complete.cases(geneList),]

geneList$rank_value <- NA
geneList$rank_value[geneList$log2FoldChange > 0] = -log10(geneList$padj)[geneList$log2FoldChange > 0]
geneList$rank_value[geneList$log2FoldChange < 0] = -log10(geneList$padj)[geneList$log2FoldChange < 0]*-1
geneList = geneList[order(geneList$rank_value, decreasing = T),]

geneRank = geneList$rank_value
names(geneRank) = geneList$symbol
geneRank = geneRank[complete.cases(geneRank)]

register(SerialParam())

res <- fgsea(pathways = m_list,
             stats = geneRank,
             minSize=10,
             maxSize=1000,
             nperm=100000,
             nproc = 18)

for(j in c("HALLMARK_INTERFERON_ALPHA_RESPONSE")){
  if(sum(j %in% res$pathway) == 1){
    p = plotEnrichment(pathway = m_list[[j]], stats = geneRank, ticksSize = 0.5) +
      labs(title = paste0("Day2 - ",j))
    print(p)
  }
}

genes.pathway = unlist(res$leadingEdge[res$pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE"])
paletteLength <- 20
myColor <- colorRampPalette(c("darkblue", "white", "darkred"))(paletteLength)
myBreaks = c(seq(-2,2,length.out = 20))

mat <- exp[,c("HI_2dpi_1","HI_2dpi_2","HI_5dpi_1","HI_5dpi_2","IOL-2dpi_1","IOL-2dpi_2","IOL-2dpi_3","IOL-2dpi_4","IOL-5dpi_1","IOL-5dpi_2","IOL-5dpi_3","IOL-5dpi_4")]
mat = mat[res_d5[order(res_d5$log2FoldChange),]$symbol,]
mat = mat[rownames(mat) %in% genes.pathway,]
pheatmap::pheatmap(mat, cluster_cols = F,cluster_rows = T, show_colnames = T, scale = "row",color = myColor, breaks = myBreaks,)

temp = (res[order(res$padj, decreasing = F),])
temp$pathway = factor(temp$pathway, levels = rev(temp$pathway))

# Save pathway table results for HALLMARK Day2.
out = temp[,1:7]
write.csv(x = out, file = paste0("~/MGN/ProcessedData/Day2_HALLMARK_statistics.csv"))

day2.levels = levels(temp$pathway)
temp = temp[1:50,]
temp$cat = temp$NES >= 0
ggplot(temp, aes(y = -log10(padj), x = pathway, fill = NES)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.15), lty=2)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0)+
  coord_flip()+
  facet_wrap(~cat, nrow = 1) +
  theme_classic()

# HALLMARK -  Day5  -------------------------------------------------------

msigdbr_show_species()m_df = msigdbr(species = "Mus musculus", category = "H") # Pick the C2 hallmark gene pathways from mDBsig
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

geneList = as.data.frame(res_d5[,c("log2FoldChange","padj","symbol")])
geneList = geneList[!duplicated(geneList$symbol),]
geneList = geneList[complete.cases(geneList),]

geneList$rank_value <- NA
geneList$rank_value[geneList$log2FoldChange > 0] = -log10(geneList$padj)[geneList$log2FoldChange > 0]
geneList$rank_value[geneList$log2FoldChange < 0] = -log10(geneList$padj)[geneList$log2FoldChange < 0]*-1
geneList = geneList[order(geneList$rank_value, decreasing = T),]

geneRank = geneList$rank_value
names(geneRank) = geneList$symbol
geneRank = geneRank[complete.cases(geneRank)]

register(SerialParam())

res <- fgsea(pathways = m_list,
             stats = geneRank,
             minSize=10,
             maxSize=1000,
             nperm=100000,
             nproc = 18)

for(j in c("HALLMARK_INTERFERON_ALPHA_RESPONSE")){
  if(sum(j %in% res$pathway) == 1){
    p = plotEnrichment(pathway = m_list[[j]], stats = geneRank, ticksSize = 0.5) +
      labs(title = paste0("Day5 - ",j))
    print(p)
  }
}

temp = (res[order(res$padj, decreasing = F),])
temp$cat = temp$NES >= 0
temp$pathway = factor(temp$pathway, levels = rev(temp$pathway))

# Save table for pathways HALLMARK Day5
out = temp[,1:7]
write.csv(x = out, file = paste0("~/MGN/ProcessedData/Day5_HALLMARK_statistics.csv"))

ggplot(temp, aes(y = -log10(padj), x = pathway, fill = NES)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.15), lty=2)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0)+
  coord_flip()+
  facet_wrap(~cat, nrow = 1) +
  theme_classic()

temp$pathway = factor(temp$pathway, levels= day2.levels) # To preserve same order as in Day2
ggplot(temp, aes(y = -log10(padj), x = pathway, fill = NES)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.15), lty=2)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0)+
  coord_flip()+
  facet_wrap(~cat, nrow = 1) +
  theme_classic()

# For REACTOME and KEGG ---------------------------------------------------
for(i in c("CP:REACTOME","CP:KEGG")){
  # HALLMARK -  Day2  -------------------------------------------------------
  msigdbr_show_species()
  m_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = i) # Pick the C2 hallmark gene pathways from mDBsig
  m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

  geneList = as.data.frame(res_d2[,c("log2FoldChange","padj","symbol")])
  geneList = geneList[!duplicated(geneList$symbol),]
  geneList = geneList[complete.cases(geneList),]

  geneList$rank_value <- NA
  geneList$rank_value[geneList$log2FoldChange > 0] = -log10(geneList$padj)[geneList$log2FoldChange > 0]
  geneList$rank_value[geneList$log2FoldChange < 0] = -log10(geneList$padj)[geneList$log2FoldChange < 0]*-1
  geneList = geneList[order(geneList$rank_value, decreasing = T),]

  geneRank = geneList$rank_value
  names(geneRank) = geneList$symbol
  geneRank = geneRank[complete.cases(geneRank)]

  register(SerialParam())

  res <- fgsea(pathways = m_list,
               stats = geneRank,
               minSize=10,
               maxSize=1000,
               nperm=100000,
               nproc = 18)


  for(j in c("KEGG_VIRAL_MYOCARDITIS","KEGG_CARDIAC_MUSCLE_CONTRACTION","KEGG_ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC","REACTOME_INNATE_IMMUNE_SYSTEM","REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION")){
    if(sum(j %in% res$pathway) == 1){
      p = plotEnrichment(pathway = m_list[[j]], stats = geneRank, ticksSize = 0.5) +
        labs(title = paste0("Day2 - ",j))
      print(p)
    }
  }
  #save(res,file =  "~/MGN/ProcessedData/Pathways_HIvsIOL_Day2Only_mergedGenome_REACTOME_v4_FINAL.Rdata")
  #load("~/MGN/ProcessedData/Pathways_HIvsIOL_Day2Only_mergedGenome_REACTOME_v4.Rdata")

  #write.table(as.data.frame(res[,1:5]), file = "~/MGN/ProcessedData/Pathways_Day2Only_REACTOME_v3.txt", quote = F, row.names = T, col.names = T)

  temp = (res[order(res$padj, decreasing = F),])
  temp$cat = temp$NES >= 0

  out = temp[,1:7]
  write.csv(x = out, file = paste0("~/MGN/ProcessedData/Day2_",gsub(i,pattern = "[:]",replacement = "_"),"_FULL_TABLE.csv"))

  temp = temp[temp$padj<0.15,]
  temp$pathway = factor(temp$pathway, levels = rev(temp$pathway))
  day2.levels = levels(temp$pathway)

  p1 =ggplot(temp, aes(y = -log10(padj), x = pathway, fill = NES)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.15), lty=2)+
    scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0)+
    coord_flip()+
    labs(title = "Day2")+
    facet_wrap(~cat, nrow = 1) +
    theme_classic() +
    ylim(0,3)
  print(p1)

  # HALLMARK -  Day5  -------------------------------------------------------

  geneList = as.data.frame(res_d5[,c("log2FoldChange","padj","symbol")])
  geneList = geneList[!duplicated(geneList$symbol),]
  geneList = geneList[complete.cases(geneList),]

  geneList$rank_value <- NA
  geneList$rank_value[geneList$log2FoldChange > 0] = -log10(geneList$padj)[geneList$log2FoldChange > 0]
  geneList$rank_value[geneList$log2FoldChange < 0] = -log10(geneList$padj)[geneList$log2FoldChange < 0]*-1
  geneList = geneList[order(geneList$rank_value, decreasing = T),]

  geneRank = geneList$rank_value
  names(geneRank) = geneList$symbol
  geneRank = geneRank[complete.cases(geneRank)]

  register(SerialParam())

  res <- fgsea(pathways = m_list,
               stats = geneRank,
               minSize=10,
               maxSize=1000,
               nperm=100000,
               nproc = 18)

  for(j in c("KEGG_VIRAL_MYOCARDITIS","KEGG_CARDIAC_MUSCLE_CONTRACTION","KEGG_ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC","REACTOME_INNATE_IMMUNE_SYSTEM","REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION")){
    if(sum(j %in% res$pathway) == 1){
      p = plotEnrichment(pathway = m_list[[j]], stats = geneRank, ticksSize = 0.5) +
        labs(title = paste0("Day5 - ",j))
      print(p)
    }
  }

  temp = (res[order(res$padj, decreasing = F),])
  temp$cat = temp$NES >= 0

  out = temp[,1:7]
  write.csv(x = out, file = paste0("~/MGN/ProcessedData/Day5_",gsub(i,pattern = "[:]",replacement = "_"),"_FULL_TABLE.csv"))

  temp = as.data.frame(temp[temp$padj<0.15,])
  temp$pathway = factor(temp$pathway, levels= day2.levels) # To preserve same order as in Day2

  p2 = ggplot(temp, aes(y = -log10(padj), x = pathway, fill = NES)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.15), lty=2)+
    scale_fill_gradient2(low = "blue", mid = "white", high = "darkred", midpoint = 0)+
    coord_flip()+
    labs(title = "Day5")+
    facet_wrap(~cat, nrow = 1) +
    ylim(0,3) +
    theme_classic()
  print(p2)
}

# Plot heatmaps for specific pathways -------------------------------------
msigdbr_show_species()
m_df = msigdbr(species = "Mus musculus", category = "C2") # Pick the C2 hallmark gene pathways from mDBsig
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#path.genes = unlist(temp[temp$pathway == "KEGG_VIRAL_MYOCARDITIS",]$leadingEdge)
#path.genes = unlist(res[res$pathway == "KEGG_CARDIAC_MUSCLE_CONTRACTION",]$leadingEdge)
path.genes = unlist(temp[temp$pathway == "KEGG_ARRHYTHMOGENIC_RIGHT_VENTRICULAR_CARDIOMYOPATHY_ARVC",]$leadingEdge)

paletteLength <- 20
myColor <- colorRampPalette(c("darkblue", "white", "darkred"))(paletteLength)
myBreaks = c(seq(-2,2,length.out = 20))

mat <- exp[,c("HI_2dpi_1","HI_2dpi_2","HI_5dpi_1","HI_5dpi_2","IOL-2dpi_1","IOL-2dpi_2","IOL-2dpi_3","IOL-2dpi_4","IOL-5dpi_1","IOL-5dpi_2","IOL-5dpi_3","IOL-5dpi_4")]
mat = mat[res_d5[order(res_d5$log2FoldChange),]$symbol,]
mat = mat[rownames(mat) %in% path.genes,]
pheatmap::pheatmap(mat, cluster_cols = F,cluster_rows = T, show_colnames = T, scale = "row",color = myColor, breaks = myBreaks,)

m_df = msigdbr(species = "Mus musculus", category = "H") # Pick the C2 hallmark gene pathways from mDBsig
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)
names(m_list)

path.genes = unlist(res[res$pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE",]$leadingEdge)

paletteLength <- 20
myColor <- colorRampPalette(c("darkblue", "white", "darkred"))(paletteLength)
myBreaks = c(seq(-2,2,length.out = 20))

mat <- exp[,c("HI_2dpi_1","HI_2dpi_2","HI_5dpi_1","HI_5dpi_2","IOL-2dpi_1","IOL-2dpi_2","IOL-2dpi_3","IOL-2dpi_4","IOL-5dpi_1","IOL-5dpi_2","IOL-5dpi_3","IOL-5dpi_4")]
mat = mat[res_d5[order(res_d5$log2FoldChange),]$symbol,]
mat = mat[rownames(mat) %in% path.genes,]
pheatmap::pheatmap(mat, cluster_cols = F,cluster_rows = T, show_colnames = T, scale = "row",color = myColor, breaks = myBreaks,)

path.genes = unlist(res[res$pathway == "HALLMARK_FATTY_ACID_METABOLISM",]$leadingEdge)

paletteLength <- 20
myColor <- colorRampPalette(c("darkblue", "white", "darkred"))(paletteLength)
myBreaks = c(seq(-2,2,length.out = 20))

mat <- exp[,c("HI_2dpi_1","HI_2dpi_2","HI_5dpi_1","HI_5dpi_2","IOL-2dpi_1","IOL-2dpi_2","IOL-2dpi_3","IOL-2dpi_4","IOL-5dpi_1","IOL-5dpi_2","IOL-5dpi_3","IOL-5dpi_4")]
mat = mat[res_d5[order(res_d5$log2FoldChange),]$symbol,]
mat = mat[rownames(mat) %in% path.genes,]
pheatmap::pheatmap(mat, cluster_cols = F,cluster_rows = T, show_colnames = T, scale = "row",color = myColor, breaks = myBreaks,)

path.genes = unlist(res[res$pathway == "HALLMARK_APOPTOSIS",]$leadingEdge)
paletteLength <- 20
myColor <- colorRampPalette(c("darkblue", "white", "darkred"))(paletteLength)
myBreaks = c(seq(-2,2,length.out = 20))

mat <- exp[,c("HI_2dpi_1","HI_2dpi_2","HI_5dpi_1","HI_5dpi_2","IOL-2dpi_1","IOL-2dpi_2","IOL-2dpi_3","IOL-2dpi_4","IOL-5dpi_1","IOL-5dpi_2","IOL-5dpi_3","IOL-5dpi_4")]
mat = mat[res_d5[order(res_d5$log2FoldChange),]$symbol,]
mat = mat[rownames(mat) %in% path.genes,]
pheatmap::pheatmap(mat, cluster_cols = F,cluster_rows = T, show_colnames = T, scale = "row",color = myColor, breaks = myBreaks,)


# To check expression of inflammasome specific genes -----------------------------------
path.genes = c("Nlrp3","Casp1","Il18","Il1b")
paletteLength <- 20
myColor <- colorRampPalette(c("darkblue", "white", "darkred"))(paletteLength)
myBreaks = c(seq(-2,2,length.out = 20))

mat <- exp[,c("HI_2dpi_1","HI_2dpi_2","HI_5dpi_1","HI_5dpi_2","IOL-2dpi_1","IOL-2dpi_2","IOL-2dpi_3","IOL-2dpi_4","IOL-5dpi_1","IOL-5dpi_2","IOL-5dpi_3","IOL-5dpi_4")]
mat = mat[res_d5[order(res_d5$log2FoldChange),]$symbol,]
mat = mat[rownames(mat) %in% path.genes,]
pheatmap::pheatmap(mat, cluster_cols = F,cluster_rows = T, show_colnames = T, scale = "row",color = myColor, breaks = myBreaks,)


# Heatmap of genes in the leading edges day 2 and 5 -----------------------

msigdbr_show_species()
m_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") # Pick the C2 hallmark gene pathways from mDBsig
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

geneList = as.data.frame(res_d2[,c("log2FoldChange","padj","symbol")])
geneList = geneList[!duplicated(geneList$symbol),]
geneList = geneList[complete.cases(geneList),]

geneList$rank_value <- NA
geneList$rank_value[geneList$log2FoldChange > 0] = -log10(geneList$padj)[geneList$log2FoldChange > 0]
geneList$rank_value[geneList$log2FoldChange < 0] = -log10(geneList$padj)[geneList$log2FoldChange < 0]*-1
geneList = geneList[order(geneList$rank_value, decreasing = T),]

geneRank = geneList$rank_value
names(geneRank) = geneList$symbol
geneRank = geneRank[complete.cases(geneRank)]

register(SerialParam())

res_day2 <- fgsea(pathways = m_list,
             stats = geneRank,
             minSize=10,
             maxSize=1000,
             nperm=100000,
             nproc = 18)

# For day 5
geneList = as.data.frame(res_d5[,c("log2FoldChange","padj","symbol")])
geneList = geneList[!duplicated(geneList$symbol),]
geneList = geneList[complete.cases(geneList),]

geneList$rank_value <- NA
geneList$rank_value[geneList$log2FoldChange > 0] = -log10(geneList$padj)[geneList$log2FoldChange > 0]
geneList$rank_value[geneList$log2FoldChange < 0] = -log10(geneList$padj)[geneList$log2FoldChange < 0]*-1
geneList = geneList[order(geneList$rank_value, decreasing = T),]

geneRank = geneList$rank_value
names(geneRank) = geneList$symbol
geneRank = geneRank[complete.cases(geneRank)]

register(SerialParam())

res_day5 <- fgsea(pathways = m_list,
                  stats = geneRank,
                  minSize=10,
                  maxSize=1000,
                  nperm=100000,
                  nproc = 18)

path.genes = unique(c(unlist(res_day2[res_day2$pathway == "REACTOME_INTERFERON_SIGNALING",]$leadingEdge),
               unlist(res_day5[res_day5$pathway == "REACTOME_INTERFERON_SIGNALING",]$leadingEdge)))

paletteLength <- 20
myColor <- colorRampPalette(c("darkblue", "white", "darkred"))(paletteLength)
myBreaks = c(seq(-2,2,length.out = 20))

mat <- exp[,c("HI_2dpi_1","HI_2dpi_2","HI_5dpi_1","HI_5dpi_2","IOL-2dpi_1","IOL-2dpi_2","IOL-2dpi_3","IOL-2dpi_4","IOL-5dpi_1","IOL-5dpi_2","IOL-5dpi_3","IOL-5dpi_4")]
mat = mat[res_d5[order(res_d5$log2FoldChange),]$symbol,]
mat = mat[rownames(mat) %in% path.genes,]
pheatmap::pheatmap(mat, cluster_cols = F,cluster_rows = T, show_colnames = T, scale = "row",color = myColor, breaks = myBreaks,)

# Count matrix for GEO
load(file = "~/MGN/ProcessedData/se.transcripts.merge.genomes.Rdata")
mtx = se@assays@data$counts
rownames(mtx) = mapIds(Mus.musculus,
                       keys=unlist(lapply(strsplit(rownames(se),"[.]"), function(x) x[1])),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
mtx = mtx[rowSums(mtx) > 0,] # remove genes with no counts
mtx = mtx[!is.na(rownames(mtx)),] # remove genes that are not mapped
dim(mtx)
View(mtx)
mtx = rbind(mtx,assay(se)[grep(rownames(se), pattern = "CHIK", value = T),]) # Add CHIKV genome counts
colnames(mtx) = unlist(lapply(strsplit(colnames(mtx),"_"), function(x) x[1])) # Clean up column names
mtx = as.data.frame(mtx)
mtx$GeneSymbol = rownames(mtx)
mtx = mtx[,c("GeneSymbol","M1-IOL-2dpi","M2-IOL-2dpi","M3-IOL-2-dpi","M4-IOL-2dpi","M5-HI-control-2dpi","M6-HI-control-2dpi","M7-IOL-5dpi","M8-IOL-5dpi","M9-IOL-5dpi","M10-IOL-5dpi","M11-HI-control-5dpi","M12-HI-control-5dpi")]

write.csv(x = mtx, file = "~/MGN/ProcessedData/GeneCounts.csv",quote = F, row.names = F, col.names = T)

