##---------------------
# Usando fgsea em GBS
# Data: 14/09/2018
##---------------------
library(tidyverse)
library(fgsea)
library(DESeq2)
library(edgeR)
library(dplyr)
library(plyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(fgsea)
library(DT)

# Arquivo com ENSEMBL IDs.
## 1. Criar coluna SYMBOL contendo símbolos dos genes 
#  associados aos IDs Ensembl.
res <- read_csv('DEGs_RecZika.csv')

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

# 
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
pathways.hallmark <- gmtPathways("h.all.v6.2.symbols.gmt")

# Se quiser ver todos de uma vez, descomente abaixo:
pathways.hallmark

# Mostrar as vias e, dentro delas, os primeiros genes. 
pathways.hallmark %>% head() %>% lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)

# Tidy de resultados:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


# Uma tabela mais interessante:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% DT::datatable()

## 3. Plotar os resultados com ggplot2:
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj <0.05)) +
  coord_flip() +
  labs(x="Vias", y="Escore Enriquecido Normalisado",
       title="Principais Vias NES de GSEA") + 
  theme_minimal()



