##----------------------------------------
# IMT 
# GBS-recuperados Vs ZIKV com Controls
# Data: 31/10/2019
# Washington C. Araujo
# GSEA - fgsea Bioconductor
##----------------------------------------
library(DESeq2)
library(edgeR)
library(dplyr)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
source("cooltable.R") #function to retrieve Up and Down genes at clusterprofiler data

# 1. Carregar o arquivo de contagem de genes (featurecounts):
countdata <- read.table("featureCounts_out.txt", header=TRUE, row.names=1, check.names = FALSE, stringsAsFactors=F)


# 2. Remover os cromossomos sexuais, na coluna de cromossomos.
countdata_sex <-countdata[!(countdata$Chr %in% c("X","Y")), ]

# 3. Remover as 5 primeiras colunas, que não serão usadas. 
countdata <- countdata_sex[ ,6:ncol(countdata_sex)]
rm(countdata_sex)

# 4. Manter os nomes das amostras, sem a extensão .BAM.
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

# 5. Filtrar amostras.
countdata <- countdata %>% dplyr::select(grep("09.", names(countdata)), #GBS_rec_pair
                                         grep("07.", names(countdata)), #GBS_rec_pair
                                         grep("16.", names(countdata)), #GBS_rec_pair
                                         grep("30.", names(countdata)), #GBS_rec_pair
                                         grep("27.", names(countdata)), #GBS_rec
                                         grep("28.", names(countdata)), #GBS_rec
                                         grep("08.", names(countdata)), #GBS_rec
                                         grep("38.", names(countdata)), #GBS_rec
                                         grep("39.", names(countdata)), #GBS_rec
                                         
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


# 6. Condition - especificar replicatas e grupos.
(condition <- factor(c(rep("GBS_REC",72), #1
                       rep("ZIKA",32), #9
                       rep("CTL",48) #10
)
)
)


# 7. Gerar o coldata.
(coldata <- data.frame(row.names=colnames(countdata), condition))

# 8. Gerar entradas para DESeq.
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

# 9. Juntar as 8 replicatas de cada amostra.
dds$sample <-  factor(c(rep("09.lane",8), 
                        rep("07.lane",8), 
                        rep("16.lane",8), 
                        rep("30.lane",8), 
                        rep("27.lane",8),
                        rep("28.lane",8),
                        rep("08.lane",8),
                        rep("38.lane",8),
                        rep("39.lane",8),
                        
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

### Observação: n Total x 8
dds$run <- paste0("run",1:152) 

ddsColl <- collapseReplicates(dds, dds$sample, dds$run, renameCols = TRUE)


# 10. Filtrar por counts insignificantes.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

rm(keep)


# 11. Rodar Real&Oficial
dds <- DESeq(ddsColl)



# 12. Criar objeto rld, que transforma os dados para log2FC. 
# Scale/fator para minimizar diferenças entre amostras.
rld <- rlogTransformation(dds)

# 13. Gerar PCA.
pca <- plotPCA(rld, ntop = nrow(counts(dds)),returnData=F)

# 14. Visualizar a PCA.
pca

# 15. Extrair os resultados da análise de DE.
contr_RECzika <- as.data.frame(results(dds, contrast=c('condition','GBS_REC','ZIKA')))

###----------------------- GSEA - Arquivo 1 ---------------------
# Criar csv para fgsea para contr_GBS_zika:
write.csv(contr_RECzika, 'contr_1_GBS-recuperados_xZika_GSEA_2019.csv')
###---------------------------------------------------------------

# 16. Criar uma nova coluna com os nomes (SYMBOLS) dos genes.
contr_RECzika$genes <- rownames(contr_RECzika)

# 17. Remoção de NAs na coluna de padj.
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


#### ------ Gene Ontology Enrichment ------ ####    
ENSEMBL_dict <- toTable(org.Hs.egENSEMBL)
SYMBOL_dict <- toTable(org.Hs.egSYMBOL)

dict <- merge(SYMBOL_dict, ENSEMBL_dict, by="gene_id", all=T)

colnames(DEGs_RECzika)[7] <- "ensembl_id"

#DEGs_RECzika <- merge(DEGs_RECzika, dict, by="ensembl_id", all.x=T, no.dups = T)

DEGs_RECzika <- join(DEGs_RECzika, dict, by = "ensembl_id", type = "left", match= "first")


GO_RECzika <- enrichGO(DEGs_RECzika$symbol,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      qvalueCutoff = 0.05
)

GO_RECzika_tab <- cool_table_GO(GO_RECzika,
                                subset(DEGs_RECzika$symbol, DEGs_RECzika$log2FoldChange > 0),
                                subset(DEGs_RECzika$symbol, DEGs_RECzika$log2FoldChange < 0)
                                )


GO_RECzika_tab <- subset(GO_RECzika_tab, qvalue < 0.05)
