##---------------------------
# Usando fgsea em GBS
# Data: 18/09/2018
# Washington Candeia
# GSEA: Hallmark
# GBS (2 pacientes) x ZIKV
##---------------------------
library(tidyverse)
library(fgsea)
library(org.Hs.eg.db)
library(ggplot2)
library(DT)

# Arquivo com ENSEMBL IDs.
## 1. Criar coluna SYMBOL contendo símbolos dos genes 
#  associados aos IDs Ensembl.
res <- read_csv('tab/contr_GBSrecZ_GSEA.csv')

# Anotações dos símbolos a partir do Ensembl.
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$ensembl_id, 
                                    columns="SYMBOL",
                                    keytype="ENSEMBL")

# Criar a coluna de símbolos a partir dos IDs Ensembl, da primeira coluna.
ens2symbol <- as_tibble(ens2symbol)

# Confirmar:
head(ens2symbol, 10)

# Unir a coluna SYMBOL ao data frame:
res <- inner_join(res, ens2symbol, by = c('ensembl_id'='ENSEMBL'))

# Confirmar:
head(res, 10)

# Pegar SYMBOL e stat e remover NAs
res2 <- res %>%  
  dplyr::select(SYMBOL, stat, padj) %>% 
  na.omit() %>% 
  distinct() %>%  
  group_by(SYMBOL) 

# Confirmando:
head(res2, 10)

## 2. Usando do fgsea.
# A função fgsea requer uma lista de conjunto de genes para checar, e
# um vetor nomeado de estatísticas em nível gênico.
ranks <- tibble::deframe(res2)

# Carregar as vias em uma lista de nomes a partir do MSigDB.
# Seria do data('examplePathways')
pathways.hallmark <- gmtPathways('symbols/h.all.v6.2.symbols.gmt')

# Se quiser ver todos de uma vez, descomente abaixo:
head(pathways.hallmark, 3)

# Mostrar as vias e, dentro delas, os primeiros genes. 
pathways.hallmark %>% head() %>% lapply(head)

# Usando a função fgsea
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=10000) #maxSize=500)

# Funcionou
head(fgseaRes, 3)

# Estatisticas:
# Fazer tabela para várias vias selecionadas.
## Fazendo contraste entre Up e Down
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = 15), pathway]

topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = 15), pathway]

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

# Plot
plotGseaTable(pathways = pathways.GObp[topPathways], 
              fgseaRes, stats=ranks, gseaParam = 0.2)


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = 30), pathway]
topPathwaysUp

topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = 30), pathway]
topPathwaysDown
# Plots de vias específicas.
# Up:

plotEnrichment(pathways.GObp[["HALLMARK_KRAS_SIGNALING_DN"]],
               stats=ranks) + labs(title="HALLMARK_KRAS_SIGNALING_DN")

plotEnrichment(pathways.GObp[["HALLMARK_NOTCH_SIGNALING"]],
               stats=ranks) + labs(title="HALLMARK_NOTCH_SIGNALING")

plotEnrichment(pathways.GObp[["HALLMARK_PANCREAS_BETA_CELLS"]],
               stats=ranks) + labs(title= "HALLMARK_PANCREAS_BETA_CELLS")

# Down:
plotEnrichment(pathways.GObp[["HALLMARK_COMPLEMENT" ]],
               stats=ranks) + labs(title="HALLMARK_COMPLEMENT" )


plotEnrichment(pathways.GObp[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
               stats=ranks) + labs(title="HALLMARK_INTERFERON_GAMMA_RESPONSE")


## 3. Plotar os resultados com ggplot2:
ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x="Vias", y="Escore Enriquecido Normalizado (NES)",
       title="Principais Vias (NES) GBS (ZIKV) Vs ZIKV (recuperados)") + 
  theme_minimal()


# collapsedPathways <- collapsePathways(fgseaRes2[order(pval)][padj < 0.01], athways=pathways.GObp, ranks)

# mainPathways <- fgseaRes2[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]

# plotGseaTable(examplePathways[mainPathways], exampleRanks, fgseaRes, gseaParam = 0.5)


