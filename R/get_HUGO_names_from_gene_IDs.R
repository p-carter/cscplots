#' Get HUGO gene names for a list of Ensembl gene IDs.
#'
#' This function uses biomaRt to look up HUGO gene names for a list of Ensembl gene IDs.
#'
#' @param gene_IDs List of gene IDs; an example of an Ensembl gene ID is 'ENSG00000000003'.
#' @return A dataframe where column 1 = ensembl gene IDs and column 2 = hgnc_symbols.
#' @export 
get_HUGO_names_from_gene_IDs <- function(gene_IDs){

  # Similar to code by Artem Solokov

  #library("biomaRt")
  
  mart        <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  ID_mappings <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = gene_IDs, mart = mart)

  no_gene_names_found <- which(ID_mappings[,2] == "")
  if(length(no_gene_names_found) > 0){
      ID_mappings <- ID_mappings[-no_gene_names_found,]
  }
  
  return(ID_mappings)
}
