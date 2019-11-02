##----------------------------------------
# IMT 
# CHIKV Vs CHIKV recuperados 
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
countdata <- countdata %>% dplyr::select(grep("22.", names(countdata)), #CHIKV
                                         grep("31.", names(countdata)), #CHIKV
                                         grep("32.", names(countdata)), #CHIKV
                                         grep("40.", names(countdata)), #CHIKV
                                         grep("18.", names(countdata)), #CHIKV_REC
                                         grep("35.", names(countdata)), #CHIKV_REC
                                         grep("41.", names(countdata)), #CHIKV_REC
                                         grep("42.", names(countdata)), #CHIKV_REC
                                         grep("46.", names(countdata)), #CHIKV_REC
                                         grep("47.", names(countdata)), #CHIKV_REC
                                         grep("10.", names(countdata)), #control
                                         grep("11.", names(countdata)), #control
                                         grep("19.", names(countdata)), #control
                                         grep("20.", names(countdata)), #control
                                         grep("33.", names(countdata)), #control
                                         grep("34.", names(countdata))  #control
)


countdata <- as.matrix(countdata)

# 6. Condition - especificar replicatas e grupos.
(condition <- factor(c(rep("CHIKV", 32),     # n=4
                       rep("CHIKV_REC", 48),  # n=6
                       rep("CTL", 48)        # n=6
)
)
)


# 7. Gerar o coldata.
(coldata <- data.frame(row.names=colnames(countdata), condition))

# 8. Gerar entradas para DESeq.
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

# 9. Juntar as 8 replicatas de cada amostra.
dds$sample <-  factor(c(rep("22.lane", 8), 
                        rep("31.lane", 8), 
                        rep("32.lane", 8),  
                        rep("40.lane", 8), 
                        rep("18.lane", 8),
                        rep("35.lane", 8),
                        rep("41.lane", 8), 
                        rep("42.lane", 8), 
                        rep("46.lane", 8), 
                        rep("47.lane", 8),
                        rep("10.lane", 8),
                        rep("11.lane", 8),
                        rep("19.lane", 8),
                        rep("20.lane", 8),
                        rep("33.lane", 8),
                        rep("34.lane", 8)
)
)

### Observação: n Total x 8
dds$run <- paste0("run",1:128) 

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
contr_CHIK_CHIKVrec <- as.data.frame(results(dds, contrast=c('condition','CHIKV','CHIKV_REC')))

###----------------------- GSEA - Arquivo 1 ---------------------
# Criar csv para fgsea para contr_CHIK_CHIKVrec:
write.csv(contr_CHIK_CHIKVrec, 'contr_4_CHIK_vs_CHIKVrec_GSEA_2019.csv')
###---------------------------------------------------------------

# 16. Criar uma nova coluna com os nomes (SYMBOLS) dos genes.
contr_CHIK_CHIKVrec$genes <- rownames(contr_CHIK_CHIKVrec)

# 17. Remoção de NAs na coluna de padj.
contr_CHIK_CHIKVrec$padj[is.na(contr_CHIK_CHIKVrec$padj)] <- 1

DEGs_GBS_zika <- subset(contr_CHIK_CHIKVrec, padj <= 0.05 & abs(log2FoldChange) > 1)

# 18.Volcanoplot
with(as.data.frame(contr_CHIK_CHIKVrec[!(-log10(contr_CHIK_CHIKVrec$padj) == 0), ]), plot(log2FoldChange,-log10(padj), pch=16, axes=T,
                                        xlim = c(-6,6), ylim = c(0,4),                                    
                                        xlab = NA, ylab = "-Log10(Pvalue-Adjusted)", main = "CHIKV vs CHIKV Recuperados"
                                                                            
)
)

with(subset(subset(as.data.frame(contr_CHIK_CHIKVrec), padj <= 0.05), log2FoldChange <= -1), points(log2FoldChange,-log10(padj), pch=21, col="black",bg = "#69B1B7"))
with(subset(subset(as.data.frame(contr_CHIK_CHIKVrec), padj <= 0.05), log2FoldChange >= 1), points(log2FoldChange,-log10(padj),pch=21, col="black",bg = "tomato3"))
abline(h=1.3,col="green", lty = 2, cex= 3)
abline(v=1,col="green", lty = 2, cex= 3)
abline(v=-1,col="green", lty = 2, cex= 3)
