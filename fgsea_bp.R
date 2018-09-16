##---------------------------
# Usando fgsea em GBS
# Data: 16/09/2018
# Washington Candeia
# GSEA: GO Biological Process
##---------------------------
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
res2 <- res %>%  
        dplyr::select(SYMBOL, stat) %>% 
        na.omit() %>% 
        distinct() %>%  
        group_by(SYMBOL) 

# Confirmando:
summary(res2)

## 2. Usando do fgsea.
# A função fgsea requer uma lista de conjunto de genes para checar, e
# um vetor nomeado de estatísticas em nível gênico.
ranks <- tibble::deframe(res2)

# Carregar as vias em uma lista de nomes a partir do MSigDB.
# Seria do data('examplePathways')
pathways.GObp <- gmtPathways('gsea_symbols/c5.bp.v6.2.symbols.gmt')

# data('examplePathways')
system.file("gen_gene_ranks.R", package="fgsea")

# Se quiser ver todos de uma vez, descomente abaixo:
head(pathways.GObp, 3)

# Mostrar as vias e, dentro delas, os primeiros genes. 
pathways.GObp %>% head() %>% lapply(head)

# Usando a função fgsea
fgseaRes <- fgsea(pathways=pathways.GObp, stats=ranks, nperm=1000)

is.list(fgseaRes)

# Tidy de resultados:
# Usar dplyr:: para:
# a. Deixar fgseaRes como tabela organizada como data frame (as_tibble);
# b. Organizar por valores descrescentes de NES [arrange(desc())].
fgseaResTidy <- fgseaRes %>%
                as_tibble() %>%
                arrange(desc(NES))   # dplyr::arrange(desc()) -> ordem descrescente.

head(fgseaResTidy, 3)

# Organizar por NES (forma crescente).
fgseaResTidy2 <- fgseaRes %>%
                as_tibble()

# Organizar por vias e pval.
fgseaResTidy_arr <- arrange(fgseaResTidy, pathway, pval)
head(fgseaResTidy_arr, 2)

# Organizar por valores descrescentes de pval.
fgseaResTidy_arr2 <- arrange(fgseaResTidy, desc(pval))
head(fgseaResTidy_arr2, 2)


# Agrupar por vias.
by_pathway <- fgseaResTidy %>% 
              group_by(pathway)

head(by_pathway, 3)

# Agrupar por NES. 
by_pathway2 <- by_pathway %>% arrange(desc(ES))
head(by_pathway2, 10)

# Agrupar por pval.
by_pathway3 <- by_pathway %>% arrange(pval)
head(by_pathway3, 3)

# Observar as 10 primeiras vias.
fgseaRes[['pathway']][1:10]

# Observar as 10 primeiras vias (Tidy).
fgseaResTidy[['pathway']][1:10]