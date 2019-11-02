##---------------------------
# Usando fgsea em GBS
# Data: 02/11/2019
# Washington Candeia
# GSEA: Hallmark
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
res <- read_csv('tab/contr_2_GBS_vs_zika_GSEA_2019.csv')
res

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
# Análise 1: H (Estados biológicos ou processos bem definidos)
pathways.hallmark <- gmtPathways('symbols/h.all.v7.0.symbols.gmt')

# C. Se quiser ver todos de uma vez, descomente abaixo (alternativa)
head(pathways.hallmark, 3)

# D. Mostrar as vias e, dentro delas, os primeiros genes. 
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

# E. Usando a função fgsea
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

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
plotGseaTable(pathways = pathways.hallmark[topPathways], 
              fgseaRes, stats=ranks, gseaParam = 0.2)


# C. Mostrar vias UP e DOWN
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = 30), pathway]
topPathwaysUp

topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = 30), pathway]
topPathwaysDown


## 5. Plotar os resultados com ggplot2:
# A. Variavel: fgseaRes
ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x="Vias", y="Escore Enriquecido Normalizado (NES)",
       title="Principais Vias de Afetados por GBS (por ZIKV) Vs ZIKV (recuperados)") + 
  theme_minimal()

# B. Variavel: fgseaResTidy
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x="Vias", y="Escore Enriquecido Normalizado (NES)",
       title="Principais Vias de Afetados por GBS (por ZIKV) Vs ZIKV (recuperados)") + 
  theme_minimal()


#### https://stephenturner.github.io/deseq-to-fgsea/ ####
# Que genes estão nestas vias
# a. Uma tabela (tibble) com todas as vias e seus genes
# b. Opcionalmente, filtrar a lista para incluir os significantes etc
pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest() %>% 
  inner_join(res, by="SYMBOL")

####################################################################################
# Plots de vias específicas (ver seção 2.H)
# Up:
plotEnrichment(pathways.hallmark[["HALLMARK_PI3K_AKT_MTOR_SIGNALING"]],
               stats=ranks) + labs(title="HALLMARK_PI3K_AKT_MTOR_SIGNALING")

# Down:
plotEnrichment(pathways.hallmark[["HALLMARK_APOPTOSIS"]],
               stats=ranks) + labs(title="HALLMARK_APOPTOSIS")



# collapsedPathways <- collapsePathways(fgseaRes2[order(pval)][padj < 0.01], athways=pathways.GObp, ranks)

# mainPathways <- fgseaRes2[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

# plotGseaTable(examplePathways[mainPathways], exampleRanks, fgseaRes, gseaParam = 0.5)


