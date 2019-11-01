##----------------------------------------
# IMT 
# GBS Vs ZIKV com Controls
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
countdata <- countdata %>% dplyr::select(grep("01.", names(countdata)), #GBS
                                         grep("02.", names(countdata)), #GBS
                                         grep("04.", names(countdata)), #GBS
                                         grep("13.", names(countdata)), #GBS
                                         grep("14.", names(countdata)), #GBS
                                         grep("25.", names(countdata)), #GBS
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
(condition <- factor(c(rep("GBS", 48),
                       rep("ZIKA", 32),
                       rep("CTL", 48)
)
)
)


# 7. Gerar o coldata.
(coldata <- data.frame(row.names=colnames(countdata), condition))

# 8. Gerar entradas para DESeq.
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

# 9. Juntar as 8 replicatas de cada amostra.
dds$sample <-  factor(c(rep("01.lane", 8), 
                        rep("02.lane", 8), 
                        rep("04.lane", 8),  
                        rep("13.lane", 8), 
                        rep("14.lane", 8),
                        rep("25.lane", 8),
                        rep("12.lane", 8), 
                        rep("24.lane", 8), 
                        rep("36.lane", 8), 
                        rep("48.lane", 8),
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
contr_GBS_zika <- as.data.frame(results(dds, contrast=c('condition','GBS','ZIKA')))

###----------------------- GSEA - Arquivo 1 ---------------------
# Criar csv para fgsea para contr_GBS_zika:
write.csv(contr_GBS_zika, 'contr_GBS_vs_zika_GSEA_2019.csv')
###---------------------------------------------------------------

# 16. Criar uma nova coluna com os nomes (SYMBOLS) dos genes.
contr_GBS_zika$genes <- rownames(contr_GBS_zika)

# 17. Remoção de NAs na coluna de padj.
contr_GBS_zika$padj[is.na(contr_GBS_zika$padj)] <- 1

DEGs_GBS_zika <- subset(contr_GBS_zika, padj <= 0.05 & abs(log2FoldChange) > 1)


# 18.Volcanoplot
with(as.data.frame(contr_GBS_zika[!(-log10(contr_GBS_zika$padj) == 0), ]), plot(log2FoldChange,-log10(padj), pch=16, axes=T,
                                        xlim = c(-6,6), ylim = c(0,4),                                    
                                        xlab = NA, ylab = "-Log10(Pvalue-Adjusted)", main = "GBS vs ZIKA"
                                                                            
)
)

with(subset(subset(as.data.frame(contr_GBS_zika), padj <= 0.05), log2FoldChange <= -1), points(log2FoldChange,-log10(padj), pch=21, col="black",bg = "#69B1B7"))
with(subset(subset(as.data.frame(contr_GBS_zika), padj <= 0.05), log2FoldChange >= 1), points(log2FoldChange,-log10(padj),pch=21, col="black",bg = "tomato3"))
abline(h=1.3,col="green", lty = 2, cex= 3)
abline(v=1,col="green", lty = 2, cex= 3)
abline(v=-1,col="green", lty = 2, cex= 3)


#### ------ Gene Ontology Enrichment ------ #### 
ENSEMBL_dict <- toTable(org.Hs.egENSEMBL)
SYMBOL_dict <- toTable(org.Hs.egSYMBOL)
dict <- merge(SYMBOL_dict, ENSEMBL_dict, by="gene_id", all=T)
colnames(DEGs_GBS_zika)[7] <- "ensembl_id"

# Adicionar novas colunas com genes_ids.
DEGs_GBS_zika <- join(DEGs_GBS_zika, dict, by = "ensembl_id", type = "left", match= "first")


GO_RECzika <- enrichGO(DEGs_GBS_zika$symbol,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      qvalueCutoff = 0.05
)

GO_RECzika_tab <- cool_table_GO(GO_RECzika,
                                subset(DEGs_GBS_zika$symbol, DEGs_GBS_zika$log2FoldChange > 0),
                                subset(DEGs_GBS_zika$symbol, DEGs_GBS_zika$log2FoldChange < 0)
                                )


GO_RECzika_tab <- subset(GO_RECzika_tab, qvalue < 0.05) # Mantendo apenas termos com qvalue < 0.05

#creating a orer for GO therms GO_GBSctl_tab$ord <- c(1:nrow(GO_GBSctl_tab))
GO_GBSctl_tab <- subset(GO_GBSrec_ZIKV_tab, GO_GBSrec_ZIKV_tab$qvalue < 0.05) #creating a orer for GO therms GO_GBSctl_tab$ord <- c(1:nrow(GO_GBSctl_tab))


GO_GBSctl_tab$ord <- c(1:nrow(GO_GBSrec_ZIKV_tab))
write.csv(GO_GBSctl_tab, 'GO_GBSctl_tab_14-09-2018.csv', row.names = FALSE)


# Retirando as colunas gene_id e 'symbol' para fgsea:
DEGs_GBS_zika$gene_id <- NULL
DEGs_GBS_zika$symbol <- NULL

###----------------------- GSEA - Arquivo 3 -----------------------------------
# Criar csv para fgsea para DEGs_GBS_zika:
# write.csv(DEGs_GBS_zika, 'DEGs_GBS_zika_18-09-2018.csv', row.names = FALSE)
###----------------------------------------------------------------------------
