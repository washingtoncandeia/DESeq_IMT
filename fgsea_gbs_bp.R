##---------------------
# Usando fgsea em GBS
# Data: 15/09/2018
# Washington Candeia
# GSEA: GO Biologica Process
##---------------------
library(tidyverse)
library(fgsea)
library(org.Hs.eg.db)
library(ggplot2)
library(DT)

# Arquivo com ENSEMBL IDs.
## 1. Criar coluna SYMBOL contendo símbolos dos genes 
#  associados aos IDs Ensembl.
res <- read_csv('contr_Zika_GBS_2pacientes_15-09-2018.csv')

# Anotações dos símbolos a partir do Ensembl.
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$ensembl_id, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")

# Criar a coluna de símbolos a partir dos IDs Ensembl, da primeira coluna.
ens2symbol <- as_tibble(ens2symbol)

# Confirmar:
ens2symbol
head(ens2symbol, 10)

# Unir a coluna SYMBOL ao data frame:
res <- inner_join(res, ens2symbol, by = c('ensembl_id'='ENSEMBL'))

# Confirmar:
head(res, 10)

# Pegar SYMBOL e stat e remover NAs
res2 <- res %>%  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% distinct() %>%  group_by(SYMBOL)

res2

## 2. Usando do fgsea.
# A função fgsea requer uma lista de conjunto de genes para checar, e
# um vetor nomeado de estatísticas em nível gênico.
# 2.1. Criar um vetor nomeado de estatísticas do teste
# tibble::deframe() para converter dados das duas colunas
# de um data frame em um vetor nomeado ou lista:
# a. Primeira coluna: nome;
# b. Segunda coluna: valor;
ranks <- tibble::deframe(res2)

# Carregar as vias em uma lista de nomes a partir do MSigDB.
# Molecular Signatures Database (MSigDB) é uma coleção 
# de conjunto de genes anotados para uso com GSEA (software).

# Ver: http://software.broadinstitute.org/gsea/downloads.jsp
# Ver: http://software.broadinstitute.org/gsea/msigdb/collection_details.jsp#H
pathways.hallmark <- gmtPathways('gsea_symbols/c5.bp.v6.2.symbols.gmt')

# Se quiser ver todos de uma vez, descomente abaixo:
head(pathways.hallmark)

# Mostrar as vias e, dentro delas, os primeiros genes. 
pathways.hallmark %>% head() %>% lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

# 10 vias iniciais:
fgseaRes[['pathway']][1:10]

# Tidy de resultados:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


# Uma tabela mais interessante:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% DT::datatable()

# Exibir colunas pathway com NES e padj:
fgseaRes[1:10,c("pathway", "NES","padj")]
fgseaRes[1:10][['pathway']]
fgseaResTidy[1:10,c("pathway", "NES","padj")]
fgseaResTidy[fgseaResTidy$pval < 0.01, ]
fgseaResTidy[fgseaResTidy$pval < 0.01, ]$pathway


# Estatisticas:
# Fazer tabela para várias vias selecionadas.
topPathwaysUP <- fgseaRes[ES > 0][head(order(pval), n = 30), pathway]
topPathwaysDOWN <- fgseaRes[ES < 0][head(order(pval), n = 30), pathway]

# Plot
topPathways <- fgseaRes[head(order(pval), n=30)][order(NES), pathway]
plotGseaTable(pathways = pathways.hallmark[topPathways], 
              fgseaRes, stats=ranks, gseaParam = 0.5)

topPathwaysUP
topPathwaysDOWN


###------------------
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

