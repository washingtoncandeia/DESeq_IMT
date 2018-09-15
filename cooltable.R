cool_table_GO <- function(clusterprofiler,up,down){
  
  # if (id == "ensembl") {
  #   col = dict$ensembl_id
  # } else if (id == "geneID") {
  #   col = dict$gene_id
  # } else {
  #   return("ERROR: id values are ensembl or geneID")
  # }
  # 
  # ENSEMBL <- toTable(org.Cf.egENSEMBL)
  # SYMBOL <- toTable(org.Cf.egSYMBOL)
  # dict <- merge(SYMBOL, ENSEMBL, by="gene_id", all=T)
  # 
  
  #splits all the rows form the colunm geneID based on the /
  splited_col <- strsplit(clusterprofiler@result$geneID, '/')
  
  #cratind tab to store the values 
  tab <- data.frame(matrix(nrow = length(splited_col), ncol = 4))
  colnames(tab) <- c("up","down", "up_genes", "down_genes")
  
  for (i in 1:length(splited_col)){
    # dict_gg <- subset(dict, col %in% splited_col[[i]])
    # dict_gg <- subset(dict_gg, !duplicated(dict_gg[,3]))
    # col_symbol <- unique(subset(dict$symbol, dict$ensembl_id %in% splited_col[[i]]))
    # up_symbol <- unique(subset(dict$symbol, col %in% up))
    # down_symbol <- unique(subset(dict$symbol, col %in% down))
    # 
    up2 <- subset(up, up %in% splited_col[[i]])
    down2 <- subset(down, down %in% splited_col[[i]])
    tab[i,] <- c(length(up2), length(down2),paste(up2, collapse = '/') , paste(down2, collapse = '/'))
    
    
  }
  
  #subsetting dict values present into the 
  #geneid <- unique(subset(dict$symbol, c(blood_up$genes,blood_down$genes) %in% dict$ensembl_id))
  new_tab <- cbind(clusterprofiler@result,tab)
  return(new_tab)
  
}
 



cool_table_KEGG <- function(clusterprofiler,up,down, id){
  
  ENSEMBL <- toTable(org.Hs.egENSEMBL)
  SYMBOL <- toTable(org.Hs.egSYMBOL)
  dict <- merge(SYMBOL, ENSEMBL, by="gene_id", all=T)
  
  if (id == "ensembl") {
    col = dict$ensembl_id
  } else if (id == "geneID") {
    col = dict$gene_id
  } else if (id == "symbol") {
    col = dict$symbol
  } else {
    return("ERROR: id values are ensembl or geneID")
  }
  
  
  #splits all the rows form the colunm geneID based on the /
  splited_col <- strsplit(clusterprofiler@result$geneID, '/')
  
  #cratind tab to store the values 
  tab <- data.frame(matrix(nrow = length(splited_col), ncol = 4))
  colnames(tab) <- c("up","down", "up_genes", "down_genes")
  
  for (i in 1:length(splited_col)){
    # dict_gg <- subset(dict, col %in% splited_col[[i]])
    # dict_gg <- subset(dict_gg, !duplicated(dict_gg[,3]))
    col_symbol <- unique(subset(dict$symbol, dict$gene_id %in% splited_col[[i]]))
    up_symbol <- unique(subset(dict$symbol, col %in% up))
    down_symbol <- unique(subset(dict$symbol, col %in% down))
    up2 <- subset(up_symbol, up_symbol %in% col_symbol)
    down2 <- subset(down_symbol, down_symbol %in% col_symbol)
    tab[i,] <- c(length(up2), length(down2),paste(up2, collapse = '/') , paste(down2, collapse = '/'))
    
  }
  
  #subsetting dict values present into the 
  #geneid <- unique(subset(dict$symbol, c(blood_up$genes,blood_down$genes) %in% dict$ensembl_id))
  new_tab <- cbind(clusterprofiler@result,tab)
  return(new_tab)
  
}



# up2 <- length(subset(up, up %in% splited_col[[i]]))
# down2 <- length(subset(down, down %in% splited_col[[i]]))
# tab[i,] <- c(up2, down2, up_symbol, down_symbol)



