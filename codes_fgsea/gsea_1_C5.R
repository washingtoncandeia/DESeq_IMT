##---------------------------
# Usando fgsea em GBS
# Data: 01/11/2019
# Washington Candeia
# GSEA: C5 GO Biological Process
# GBS (recuperados) x ZIKV
##---------------------------
library(tidyverse)
library(fgsea)
library(org.Hs.eg.db)
library(ggplot2)
library(DT)

# Arquivo com ENSEMBL IDs.
## 1. Criar coluna SYMBOL contendo símbolos dos genes 
#  associados aos IDs Ensembl.
res <- read_csv('tab/contr_1_GBS-recuperados_x_Zika_GSEA_2019.csv')

# A. Anotações dos símbolos a partir do Ensembl.
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$ensembl_id, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")

# B. Criar a coluna de símbolos a partir dos IDs Ensembl, da primeira coluna.
ens2symbol <- as_tibble(ens2symbol)
ens2symbol

# C. Confirmar:
head(ens2symbol, 10)

# D. Unir a coluna SYMBOL ao data frame:
res <- inner_join(res, ens2symbol, by = c('ensembl_id'='ENSEMBL'))

# E. Confirmar:
head(res, 10)

# F. Pegar SYMBOL e stat e remover NAs
res2 <- res %>%  
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>%  
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
res2

# G. Confirmando:
head(res2, 10)
summary(res2)

## 2. Usando do fgsea.
# A. A função fgsea requer uma lista de conjunto de genes para checar, e
# um vetor nomeado de estatísticas em nível gênico.
ranks <- tibble::deframe(res2)
head(ranks, 20)

# B. Carregar as vias em uma lista de nomes a partir do MSigDB.
# Análise 1: H
pathways.GObp <- gmtPathways('symbols/c5.bp.v7.0.symbols.gmt')

# C. Se quiser ver todos de uma vez, descomente abaixo (alternativa)
head(pathways.GObp, 3)

# D. Mostrar as vias e, dentro delas, os primeiros genes. 
pathways.GObp %>% 
  head() %>% 
  lapply(head)

# E. Usando a função fgsea
fgseaRes <- fgsea(pathways=pathways.GObp, stats=ranks, nperm=1000)

# F. Verificar
head(fgseaRes, 3)

# G. Tidy dos resultados
# Notar que na tabela ggplot2 pode-se adicionar a variável
# fgseaRes (seção 2.E) ou fgseaResTidy (seção 2.G, código abaixo)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Tabela javascript:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

# 4. Estatisticas fgsea:
# A. Fazer tabela para várias vias selecionadas.
# Fazendo contraste entre Up e Down
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = 30), pathway]

topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = 30), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

# B. Plotar
plotGseaTable(pathways = pathways.GObp[topPathways], 
              fgseaRes, stats=ranks, gseaParam = 0.2)


# C. Mostrar vias UP e DOWN
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = 30), pathway]
topPathwaysUp

topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = 30), pathway]
topPathwaysDown

# Plots de vias específicas.
# Up:
plotEnrichment(pathways.GObp[["GO_BRAIN_MORPHOGENESIS" ]],
               stats=ranks) + labs(title="GO_BRAIN_MORPHOGENESIS" )

plotEnrichment(pathways.GObp[["GO_DENDRITE_EXTENSION" ]],
               stats=ranks) + labs(title="GO_DENDRITE_EXTENSION" )

# Down:
plotEnrichment(pathways.GObp[[ "GO_APOPTOTIC_PROCESS" ]],
               stats=ranks) + labs(title= "GO_APOPTOTIC_PROCESS"  )

plotEnrichment(pathways.GObp[[ "GO_REGULATION_OF_IMMUNE_RESPONSE"  ]],
               stats=ranks) + labs(title= "GO_REGULATION_OF_IMMUNE_RESPONSE"  )



