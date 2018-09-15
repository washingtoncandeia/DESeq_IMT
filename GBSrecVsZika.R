##----------------------------------------
# IMT 
# GBS-rec Vs ZIKV com Controls
# Data: 12/09/2018
# Viviane Nogueira, Diego, Washington
##----------------------------------------
library(DESeq2)
library(edgeR)
library(dplyr)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
source("cooltable.R") #function to retrieve Up and Down genes at clusterprofiler data

###### preparing data set ####

#reading counts table
#countdata <- read.table("~/Desktop/GBS_NEW/featureCounts_out2.txt", header=TRUE, row.names=1, check.names = FALSE)
countdata <- read.table("featureCounts_out.txt", header=TRUE, row.names=1, check.names = FALSE, stringsAsFactors=F)


#remove sexual chromosomes
countdata_sex <-countdata[!(countdata$Chr %in% c("X","Y")), ]

#remove the 5 first colunms 
countdata <- countdata_sex[ ,6:ncol(countdata_sex)]
rm(countdata_sex)
#tira os .bam do nome das amostras
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

#filtro pros dados a serem analisados
countdata <- countdata %>% dplyr::select(grep("07.", names(countdata)), #GBS_rec_pair
                                         grep("16.", names(countdata)), #GBS_rec_pair
                                         
                                         grep("12.", names(countdata)), #12_ZIKA
                                         grep("24.", names(countdata)), #24_ZIKA
                                         grep("36.", names(countdata)), #36_ZIKA
                                         grep("48.", names(countdata)), #48_ZIKA
                                         
                                         grep("10.", names(countdata)), #control
                                         grep("11.", names(countdata)), #control
                                         grep("19.", names(countdata)), #control
                                         grep("20.", names(countdata)), #control
                                         grep("33.", names(countdata)), #control
                                         grep("34.", names(countdata))  #control
)


countdata <- as.matrix(countdata)


#condition - especificar replicatas e grupos
#(condition <- factor(c(rep("disease", 48), rep("postrecovery", 48))))
(condition <- factor(c(rep("REC",16), #1
                       rep("ZIKA",32), #9
                       rep("CTL",48) #10
)
)
)


#gerar o coldata
(coldata <- data.frame(row.names=colnames(countdata), condition))

#gerar o input pro DESeq
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

#juntar as replicatas
dds$sample <-  factor(c(rep("07.lane",8), 
                        rep("16.lane",8), 
                        
                        rep("12.lane",8), 
                        rep("24.lane",8), 
                        rep("36.lane",8), 
                        rep("48.lane",8),
                        
                        rep("10.lane",8),
                        rep("11.lane",8),
                        rep("19.lane",8),
                        rep("20.lane",8),
                        rep("33.lane",8),
                        rep("34.lane",8)
)
)


dds$run <- paste0("run",1:96) 

ddsColl <- collapseReplicates(dds, dds$sample, dds$run, renameCols = TRUE)


#filtrar counts insignificantes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

rm(keep)


#Rodar Real&Oficial
dds <- DESeq(ddsColl)



# creating the object rld which transform the data into a log2 
# scale to minimizes differences between the samples
rld <- rlogTransformation(dds)

#performing the PCA
pca <- plotPCA(rld, ntop = nrow(counts(dds)),returnData=F)
pca

#extracting the results from Differential expression analysis 
contr_RECzika <- as.data.frame(results(dds, contrast=c('condition','REC','ZIKA')))

#creating a new column with the gene SYMBOL name
contr_RECzika$genes <- rownames(contr_RECzika)
#removing the NAs from padj column
contr_RECzika$padj[is.na(contr_RECzika$padj)] <- 1

DEGs_RECzika <- subset(contr_RECzika, padj <= 0.05 & abs(log2FoldChange) > 1)


#Volcanoplot
with(as.data.frame(contr_RECzika[!(-log10(contr_RECzika$padj) == 0), ]), plot(log2FoldChange,-log10(padj), pch=16, axes=T,
                                                                              xlim = c(-6,6), ylim = c(0,4),                                    
                                                                              xlab = NA, ylab = "-Log10(Pvalue-Adjusted)", main = "GBS_REC vs ZIKA"
                                                                              
)
)

with(subset(subset(as.data.frame(contr_RECzika), padj <= 0.05), log2FoldChange <= -1), points(log2FoldChange,-log10(padj), pch=21, col="black",bg = "#69B1B7"))
with(subset(subset(as.data.frame(contr_RECzika), padj <= 0.05), log2FoldChange >= 1), points(log2FoldChange,-log10(padj),pch=21, col="black",bg = "tomato3"))
abline(h=1.3,col="green", lty = 2, cex= 3)
abline(v=1,col="green", lty = 2, cex= 3)
abline(v=-1,col="green", lty = 2, cex= 3)



#enrichment
ENSEMBL_dict <- toTable(org.Hs.egENSEMBL)
SYMBOL_dict <- toTable(org.Hs.egSYMBOL)

dict <- merge(SYMBOL_dict, ENSEMBL_dict, by="gene_id", all=T)

colnames(DEGs_RECzika)[7] <- "ensembl_id"

#DEGs_RECzika <- merge(DEGs_RECzika, dict, by="ensembl_id", all.x=T, no.dups = T)

DEGs_RECzika <- join(DEGs_RECzika, dict, by = "ensembl_id", type = "left", match= "first")
# Criar tabela que possa ser usada com fgsea:
write.csv(DEGs_RECzika, 'DEGs_RECzika_14-09-2018.csv', row.names = FALSE)


GO_RECzika <- enrichGO(DEGs_RECzika$symbol,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = 'BP',
                       qvalueCutoff = 0.05
                       )

View(summary(GO_RECzika))
# Criar uma tabela para análises posteriores.
write.csv(GO_RECzika, 'GO_RECzika_14-09-2018.csv', row.names = FALSE)

#recovering the expression status for each GO 
GO_GBSrec_ZIKV_tab <- cool_table_GO(GO_RECzika, 
                                    rownames(subset(DEGs_RECzika, log2FoldChange > 0)),
                                    rownames(subset(DEGs_RECzika, log2FoldChange < 0)) 
                               ) 

View(GO_GBSrec_ZIKV_tab)
# Criar uma tabela para análises posteriores.
write.csv(GO_GBSrec_ZIKV_tab, 'GO_GBSrec_ZIKV_tab_14-09-2018.csv', row.names = FALSE)

#keeping only those therms with qvalue < 0.05 
GO_GBSctl_tab <- subset(GO_GBSrec_ZIKV_tab, GO_GBSrec_ZIKV_tab$qvalue < 0.05) #creating a orer for GO therms GO_GBSctl_tab$ord <- c(1:nrow(GO_GBSctl_tab))


GO_GBSctl_tab$ord <- c(1:nrow(GO_GBSrec_ZIKV_tab))
write.csv(GO_GBSctl_tab, 'GO_GBSctl_tab_14-09-2018.csv', row.names = FALSE)

# Formas de salvar as tabelas de GO.
write.table(as.data.frame(contr_RECzika), sep = ",", "contr_RECzika.txt", row.names = FALSE)
