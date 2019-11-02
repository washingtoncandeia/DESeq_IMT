##---------------------------
# Usando fgsea em GBS
# Data: 22/09/2018
# Washington Candeia
# GSEA: GO Biological Process
# GBS (doentes) x ZIKV
##---------------------------
library(tidyverse)
library(fgsea)
library(org.Hs.eg.db)
library(ggplot2)
library(DT)

# Arquivo com ENSEMBL IDs.
## 1. Criar coluna SYMBOL contendo símbolos dos genes 
#  associados aos IDs Ensembl.
res <- read_csv('tab/contr_GBS_vs_zika_GSEA.csv')

# Anotações dos símbolos a partir do Ensembl.
# Cria um data frame (ens2symbol) contendo colunas SYMBOL e ENSEMBL.
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$ensembl_id, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")

# Criar a coluna de símbolos a partir dos IDs Ensembl da primeira coluna.
# as_tibble: Coerce lists and matrices to data frames.
# https://tibble.tidyverse.org/
# Tibble é uma forma modernizada de data frames.
ens2symbol <- as_tibble(ens2symbol)

# Confirmar:
head(ens2symbol, 10)

# Unir a coluna SYMBOL ao data frame a ser analisado (res).
# inner_join: Funçao de dplyr usada para juntar tiblle 'ens2symbol' com 'res'.
res <- inner_join(res, ens2symbol, by = c('ensembl_id'='ENSEMBL'))

# Confirmar:
head(res, 10)

# Pegar SYMBOL e stat e remover NAs
res2 <- res %>%  
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>%  
  group_by(SYMBOL) 

# Confirmando:
head(res2, 3)
summary(res2)

## 2. Usando do fgsea.
# A função fgsea requer uma lista de conjunto de genes para checar, e
# um vetor nomeado de estatísticas em nível gênico.
ranks <- tibble::deframe(res2)

# Carregar as vias em uma lista de nomes a partir do MSigDB.
# Seria do data('examplePathways')
pathways.GObp <- gmtPathways('symbols/c5.bp.v6.2.symbols.gmt')

# Se quiser ver todos de uma vez, descomente abaixo:
head(pathways.GObp, 3)

# Mostrar as vias e, dentro delas, os primeiros genes. 
pathways.GObp %>% head() %>% lapply(head)

# Usando a função fgsea
fgseaRes <- fgsea(pathways=pathways.GObp, stats=ranks, nperm=10000) #maxSize=500)

# Funcionou
head(fgseaRes, 3)


# Estatisticas:
## Fazendo contraste entre Up e Down
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n= 15), pathway]

topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n= 15), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

# Plot
plotGseaTable(pathways = pathways.GObp[topPathways], 
              fgseaRes, stats=ranks, gseaParam = 0.2)


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n= 15), pathway]
topPathwaysUp
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n= 15), pathway]
topPathwaysDown

# Plots de vias específicas.
# Up:
plotEnrichment(pathways.GObp[["GO_REGULATION_OF_TRANSPORT"]],
               stats=ranks) + labs(title='8-GO_REGULATION_OF_TRANSPORT')

plotEnrichment(pathways.GObp[["GO_NEUROGENESIS"]],
               stats=ranks) + labs(title="6-GO_NEUROGENESIS" )

# Down:
plotEnrichment(pathways.GObp[["GO_REGULATION_OF_TELOMERASE_RNA_LOCALIZATION_TO_CAJAL_BODY"]],
               stats=ranks) + labs(title="GO_REGULATION_OF_TELOMERASE_RNA_LOCALIZATION_TO_CAJAL_BODY"  )


plotEnrichment(pathways.GObp[["GO_POSITIVE_REGULATION_OF_VIRAL_TRANSCRIPTION"]],
               stats=ranks) + labs(title="GO_POSITIVE_REGULATION_OF_VIRAL_TRANSCRIPTION")




# collapsedPathways <- collapsePathways(fgseaRes2[order(pval)][padj < 0.01], athways=pathways.GObp, ranks)

# mainPathways <- fgseaRes2[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

# plotGseaTable(examplePathways[mainPathways], exampleRanks, fgseaRes, gseaParam = 0.5)



##------------------ Se quiser organizar os resultados ------------------
# Tidy de resultados (apenas para observar em uma tabela DT).
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Uma tabela mais interessante.
# Resultará em uma tabela que se utiliza de javascript.
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% DT::datatable()
  
## Testes sem relaçao com análise:
# Exibir colunas pathway com NES e padj:
fgseaRes[1:10,c("pathway", "NES","padj")]
fgseaRes[1:10][['pathway']]
fgseaResTidy[1:10,c("pathway", "NES","padj")]
fgseaResTidy[fgseaResTidy$pval < 0.01, ]
fgseaResTidy[fgseaResTidy$pval < 0.01, ]$pathway

# Selecionar por pval < 0.01:
fgseaResTidy[fgseaResTidy$pval < 0.01, ]$pathway
head(fgseaRes[order(pval), ], 50)

# Selecionar por padj < 0.05:
fgseaResTidy[fgseaResTidy$padj < 0.05, ]$pathway
# Ordenar até 61 por pval:
fgsea.table <- head(fgseaRes[order(pval), ], 61)

df <- data.frame(matrix(unlist(fgsea.table), byrow=T))
df <- data.frame(matrix(unlist(fgsea.table), byrow=T),stringsAsFactors=FALSE)
head(df, 60)
write.table(df, 'dataframe.txt')

# Tentar usar somente os 60 primeiros:
fgseaResTidy[fgseaResTidy$pval <= 0.001338688, ]$pathway
fgseaResTidy[fgseaResTidy$pval, ]$pathway
##----------------------------------------------------------------
