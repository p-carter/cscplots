#' For all HUGO gene groups find all their members positions in a list of genes.
#'
#' This function takes a list of gene names, and a named list of vectors, with each 
#'     vector describing a HUGO gene group by containing the HUGO gene symbols for each
#'     groups' members. The function converts the named list of vectors, so that the gene names are
#'     replaced with their corresponding positions in the gene list provided as the first param.
#'     This mapping of HUGO gene groups to their positions in the list of genes is
#'     used in the function printExpressionBarGraphs_for_families to produce plots of preferential
#'     expression of HUGO gene groups.
#'
#' @param gene_names List of HUGO gene names, e.g. list of column names in an scRNA-Seq
#'     expression matrix. 
#' @param HUGO_groups_list Named list of vectors, with each vector describing a HUGO gene group 
#'     by containing the HUGO gene symbols for each groups' members.
#' @return HUGO_groups_positions Named list of vectors, with each vector containing the 
#'      positions of the HUGO gene group members in the gene_names list.
#' @export
map_gene_positions_to_HUGO_gene_groups <- function(gene_names, HUGO_groups_list){

    HUGO_groups_positions <- list()
    for(i in 1:length(HUGO_groups_list)){
        these_positions <- vector()
    
        this_enumerated_group_name <- names(HUGO_groups_list)[i]
    
        for(j in 1:length(HUGO_groups_list[[i]])){
            this_enumerated_name <- HUGO_groups_list[[i]][j]
            for(k in 1:length(gene_names)){
                this_transcript_name <- gene_names[k]
    
                if(this_enumerated_name == this_transcript_name){
    		these_positions[length(these_positions)+1] <- k
                }	    
            }
        }
        if(length(these_positions) > 0){
    	HUGO_groups_positions[[(length(HUGO_groups_positions) + 1)]] <- these_positions
            names(HUGO_groups_positions)[length(HUGO_groups_positions)] <- this_enumerated_group_name
        } else {
            #cat(paste0(this_enumerated_group_name, " not found\n"))
        }
    }

    return(HUGO_groups_positions)
}
