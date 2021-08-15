#' Get gene lengths using HUGO gene names.
#'
#' This function uses biomaRt to obtain chromosomomal location which are used to obtain the
#'     full gene lengths for each of a vector of gene names i.e. HGNC symbols.
#'
#' @param gene_names Vector of gene names to obtain the gene lengths for. 
#' @return Returned value genes_lengths is a named vector of integers containing
#'     the full gene lengths for the continuous region between the chromsomal start and
#'     stop locations provided by biomaRt/Ensembl.
#' @export
get_gene_lengths_using_gene_names <- function(gene_names, type = NULL){

    library("biomaRt")
    mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

    # Chromosomal locations. From biomaRt documentation: Gene Start (bp);  Gene End (bp)
    gene_coords_df = getBM(attributes = c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters = "hgnc_symbol", values = gene_names, mart = mart)
    gene_coords_df$size = (gene_coords_df$end_position - gene_coords_df$start_position) + 1

    gene_lengths <- vector()
    gene_lengths <- rep(NA, length(gene_names))
 
    gene_lengths <- gene_coords_df[(match(gene_names, gene_coords_df[,1])),5]
    names(gene_lengths) <- gene_names

    return(gene_lengths)
}

#' Get transcript lengths using HUGO gene names.
#'
#' This function uses biomaRt to obtain transcript lengths for each of a vector of gene names
#'     (HGNC symbols).
#'
#' @param gene_names Vector of gene names to obtain the transcript lengths for. 
#' @param type Optional: 4 possible choices are:
#'     1) "oldest" = get the lengths of the oldest transcript lengths - this uses the external
#'     transcript names to indicate which transcript was first publicly deposited e.g. out of A1BG-201,
#'     A1BG-202, A1BG-203, A1BG-204 and A1BG-205, the oldest transcript is taken as A1BG-201.
#'     2) "longest" = get always the longest transcript out of all possibilites.
#'     3) "appris or oldest" = if any of the transcripts for a gene has a BioMart APPRIS annotation
#'     of 'principal1' use that transcript length, or if more than one has that annotation use the
#'     average of them, otherwise use the oldest transcript. This is the default.
#'     4) "appris or longest" = use transcripts annotated by APPRIS as 'principal1' otherwise
#'     use the longest transcript.
#' @return Returned value transcript_lengths is a named vector of integers.
#' @export
get_transcript_lengths_using_gene_names <- function(gene_names, type = NULL){

    library("biomaRt")
    mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

    # From biomaRt documentation: Transcript length (including UTRs and CDS); Transcript Start (bp); Transcript End (bp); APPRIS annotation; Associated Transcript Name
    trans_info <- getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_length", "transcript_start", "transcript_end", "transcript_appris",
                                       "external_transcript_name"),
                                       filters = "hgnc_symbol", values = gene_names, mart = mart)

    trans_info$transcript_coords_size = ((trans_info$transcript_end - trans_info$transcript_start) + 1)
    
    unique_gene_names <- unique(trans_info[,1])
    unique_gene_names_df <- data.frame(unique_gene_names, stringsAsFactors = FALSE)
    
    for(i in 1:nrow(unique_gene_names_df)){
    
        this_gene <- unique_gene_names_df[i,1]
    
        this_set_of_transcripts = '';
        this_set_of_transcripts = which(trans_info$hgnc_symbol == this_gene)    
    
        oldest_transcript_pos = ''
        oldest_transcript_order_pos = ''
        oldest_transcript_order_pos <- order(trans_info[this_set_of_transcripts,7])[1]
        oldest_transcript_pos = this_set_of_transcripts[oldest_transcript_order_pos]
        oldest_transcript_len <- trans_info[oldest_transcript_pos,3]
        unique_gene_names_df[i,2] <- oldest_transcript_len
    
        longest_transcript_pos = ''
        longest_transcript_order_pos = ''
        longest_transcript_order_pos <- rev(order(trans_info[this_set_of_transcripts,3]))[1]
        longest_transcript_pos = this_set_of_transcripts[longest_transcript_order_pos]
        longest_transcript_len <- trans_info[longest_transcript_pos,3]
        unique_gene_names_df[i,3] <- longest_transcript_len
    
        trans_info[this_set_of_transcripts,]
        appris_principal1_transcript_pos = ''
        appris_principal1_transcript_pos <- which(trans_info[this_set_of_transcripts,6] == "principal1")
    
        appris_principal1_or_oldest_transcript_len = ''
        appris_principal1_or_longest_transcript_len = ''
        if(length(appris_principal1_transcript_pos) == 0){
            appris_principal1_or_oldest_transcript_len = oldest_transcript_len
            appris_principal1_or_longest_transcript_len = longest_transcript_len
        } else if(length(appris_principal1_transcript_pos) == 1){
            appris_principal1_or_oldest_transcript_len <- trans_info[this_set_of_transcripts[appris_principal1_transcript_pos[1]],3]
            appris_principal1_or_longest_transcript_len <- trans_info[this_set_of_transcripts[appris_principal1_transcript_pos[1]],3]
        } else {
            #appris_principal1_or_oldest_transcript_len  <- round(mean(trans_info[this_set_of_transcripts[appris_principal1_transcript_pos],3]), 2)
            #appris_principal1_or_longest_transcript_len <- round(mean(trans_info[this_set_of_transcripts[appris_principal1_transcript_pos],3]), 2)
            appris_principal1_or_oldest_transcript_len  <- round(mean(trans_info[this_set_of_transcripts[appris_principal1_transcript_pos],3]), 0)
            appris_principal1_or_longest_transcript_len <- round(mean(trans_info[this_set_of_transcripts[appris_principal1_transcript_pos],3]), 0)
        }
        unique_gene_names_df[i,4] <- appris_principal1_or_oldest_transcript_len
        unique_gene_names_df[i,5] <- appris_principal1_or_longest_transcript_len    
    }
    
    colnames(unique_gene_names_df) <- c("hgnc_symbol", "oldest_transcript_len", "longest_transcript_len", "appris_transcript_len_or_oldest", "appris_transcript_len_or_longest")
    unique_gene_names_transcript_lengths_df <- unique_gene_names_df

    #save(unique_gene_names_transcript_lengths_df, file = paste0("../data/unique_gene_names_transcript_lengths_df.RData"))    
    #write.table(unique_gene_names_transcript_lengths_df, file = paste0("gene_lengths-oldest_longest_appris.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", append = FALSE, quote = FALSE)

    if(is.null(type)){
        unique_gene_names_df_col_index <- 4
    } else if(type == "oldest"){
        unique_gene_names_df_col_index <- 2
    } else if(type == "longest"){
        unique_gene_names_df_col_index <- 3
    } else if(type == "appris or oldest"){
        unique_gene_names_df_col_index <- 4
    } else if(type == "appris or longest"){
        unique_gene_names_df_col_index <- 5
    } else {
        unique_gene_names_df_col_index <- 4
    }

    transcript_lengths <- vector()
    transcript_lengths <- rep(NA, length(gene_names))
    
    transcript_lengths <- unique_gene_names_df[(match(gene_names, unique_gene_names_df[,1])),unique_gene_names_df_col_index]
    names(transcript_lengths) <- gene_names
    
    return(transcript_lengths)
}
