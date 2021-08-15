#' Convert a vector of gene names into their corresponding full HUGO gene names.
#'
#' This function uses hgnc_complete_set.txt from the EBI to convert a vector of gene names
#'     into their corresponding full HUGO gene names.
#'
#' @param gene_names Vector of gene names to be converted into full gene names. 
#' @return Returned value genes_full_names is a vector of full genes names for the gene_names
#'     parameter provided.
#' @export
getGenesFullNames <- function(genes_names, all_HUGO_abbreviations_and_fullnames = NULL, filesDir = NULL){
    
    if(is.null(all_HUGO_abbreviations_and_fullnames)){
        if(!file.exists(paste0(filesDir, "hgnc_complete_set.txt"))){
            HUGO_ftp_address <- 'ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt'
            download.file(HUGO_ftp_address, paste0(filesDir, "hgnc_complete_set.txt"))
        }
        all_HUGO_info <- read.csv(paste0(filesDir, "hgnc_complete_set.txt"), header=TRUE, sep="\t", stringsAsFactors = FALSE)
        all_HUGO_abbreviations_and_fullnames <- all_HUGO_info[,c(2,3)]
    }
    if(file.exists(paste0(filesDir, "custom_gene_fullnames.txt"))){ # manual additions not detected automatically
        custom_gene_fullnames <- read.csv(file = paste0(filesDir, "custom_gene_fullnames.txt"), header = T, stringsAsFactors = FALSE, sep = "\t")
        all_HUGO_abbreviations_and_fullnames <- all_HUGO_abbreviations_and_fullnames[-which(all_HUGO_abbreviations_and_fullnames[,1] %in% custom_gene_fullnames[,1]),] # remove duplicates
        #custom_gene_fullnames <- custom_gene_fullnames[-which(custom_gene_fullnames[,1] %in% all_HUGO_abbreviations_and_fullnames[,1]),] # use this instead to give HUGO preference over custom
        all_HUGO_abbreviations_and_fullnames_to_use <- rbind(all_HUGO_abbreviations_and_fullnames, custom_gene_fullnames)
        all_HUGO_abbreviations_and_fullnames_to_use <- all_HUGO_abbreviations_and_fullnames_to_use[order(all_HUGO_abbreviations_and_fullnames_to_use[,1]),]
    } else {
        all_HUGO_abbreviations_and_fullnames_to_use <- all_HUGO_abbreviations_and_fullnames
    }
    
    # Find full names for pCSCs subset genes
    match_sequence_tmp <- match(genes_names, all_HUGO_abbreviations_and_fullnames_to_use[,1]) # matches of sorted gene names from pCSCs subset found in full names
    match_sequence_not_na_pos <- which(!is.na(match_sequence_tmp))
    match_pos <- match_sequence_tmp[match_sequence_not_na_pos]
    
    genes_abbrevs_found <- all_HUGO_abbreviations_and_fullnames_to_use[match_sequence_tmp,1]
    genes_full_names <- genes_names
    genes_full_names[match_sequence_not_na_pos] <-
        paste0(all_HUGO_abbreviations_and_fullnames_to_use[match_pos, 2], " (", genes_abbrevs_found[match_sequence_not_na_pos], ")")
    # If wanting to manually add missing full names to custom_gene_fullnames.txt, can print the missing out
    missing_names_genes_full_names_print_str <- paste0(sort(genes_full_names[which(is.na(match_sequence_tmp))]), collapse = '\n')
    #cat(paste0("No matching full names in genes_full_names for:\n", missing_names_genes_full_names_print_str, "\n", collapse = ", "))

    return(genes_full_names)
}
    
#' Convert a vector of gene names into a data frame containing gene names and matching full HUGO gene names.
#'
#' This function uses hgnc_complete_set.txt from the EBI to convert a vector of gene names
#'     into a DF containing the gene names and their corresponding full HUGO gene names (or '-' if not found).
#'
#' @param gene_names Vector of gene names to be converted into full gene names. 
#' @return Returned value genes_full_names is a data frame containing the gene names in the first column
#'     and the full genes names in the second column.
#' @export
getGenesAndFullNamesDF <- function(genes_names, all_HUGO_abbreviations_and_fullnames = NULL, filesDir = NULL){

    # Get the full gene names/functions for printing
    if(is.null(all_HUGO_abbreviations_and_fullnames)){
        if(!file.exists(paste0(filesDir, "hgnc_complete_set.txt"))){
            HUGO_ftp_address <- 'ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt'
            download.file(HUGO_ftp_address, paste0(filesDir, "hgnc_complete_set.txt"))
        }
        all_HUGO_info <- read.csv(paste0(filesDir, "hgnc_complete_set.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        all_HUGO_abbreviations_and_fullnames <- all_HUGO_info[,c(2,3)]
    }
    if(file.exists(paste0(filesDir, "custom_gene_fullnames.txt"))){ # manual additions not detected automatically
        custom_gene_fullnames <- read.csv(file = paste0(filesDir, "custom_gene_fullnames.txt"), header = T, stringsAsFactors = FALSE, sep = "\t")
        all_HUGO_abbreviations_and_fullnames <- all_HUGO_abbreviations_and_fullnames[-which(all_HUGO_abbreviations_and_fullnames[,1] %in% custom_gene_fullnames[,1]),] # remove duplicates
        #custom_gene_fullnames <- custom_gene_fullnames[-which(custom_gene_fullnames[,1] %in% all_HUGO_abbreviations_and_fullnames[,1]),] # use this instead to give HUGO preference over custom
        all_HUGO_abbreviations_and_fullnames_to_use <- rbind(all_HUGO_abbreviations_and_fullnames, custom_gene_fullnames)
        all_HUGO_abbreviations_and_fullnames_to_use <- all_HUGO_abbreviations_and_fullnames_to_use[order(all_HUGO_abbreviations_and_fullnames_to_use[,1]),]
    } else {
        all_HUGO_abbreviations_and_fullnames_to_use <- all_HUGO_abbreviations_and_fullnames
    }

    match_pos <- match(gene_names_to_print, all_HUGO_abbreviations_and_fullnames_to_use[,1])
    gene_names_to_print_and_functions <- data.frame(
        hgnc_symbol = gene_names_to_print,
        description = all_HUGO_abbreviations_and_fullnames_to_use[match_pos,2],
        stringsAsFactors = FALSE
    )
    gene_names_to_print_and_functions[which(is.na(gene_names_to_print_and_functions[,2])),2] <- '-'

    return(gene_names_to_print_and_functions)
}
