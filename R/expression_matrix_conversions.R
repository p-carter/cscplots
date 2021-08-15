#' Convert a matrix of expression values to a matrix of count fraction values.
#'
#' This function converts an expression matrix (e.g. of TPM values), to
#' a matrix of count fractions (e.g. TPM Count Fractions, if the input was TPMs). The 
#' count fraction matrix can be helpful in identifying preferential expression e.g. in
#' a group of cancer stem cells compared to a tumour population with/without cells of
#' others types.
#'
#' @param counts Expression matrix describing a cancer single-cell population; column
#'     names should be gene names, row names should be cancer cell names.
#' @return Expression matrix containing count fraction values.
#' @export
count_fractions_transcriptwise <- function(counts){

    count_fractions <- counts
    count_fractions[count_fractions > 0] <- 0

    for(i in 1:ncol(counts)){
        #if(i %% 1000 == 0){
        #    if(i == 1000){
        #        cat("Progress: ")
        #    }
        #    cat(paste0(i, "."))
        #}
        this_col_vec_of_counts <- counts[,i]
        this_col_vec_of_count_fractions <- vector()        
        if(sum(this_col_vec_of_counts) == 0){
            this_col_vec_of_count_fractions <- this_col_vec_of_counts
        } else {
            N <- sum(this_col_vec_of_counts)            
            for(j in 1:length(this_col_vec_of_counts)){
                this_count_fraction <- ((this_col_vec_of_counts[j] / N) * 1e6)
                this_col_vec_of_count_fractions[j] <- this_count_fraction
            }
        }
        count_fractions[,i] <- this_col_vec_of_count_fractions
    }
    #cat("..Completed\n")
    
    return(count_fractions)
}

#' Convert a matrix of read counts to a matrix of Cells Per Million values.
#'
#' This function converts an expression matrix (e.g. raw read counts), to
#' a matrix of Cells Per Million. The CellsPM matrix can be helpful in identifying
#' preferential expression e.g. in a group of cancer stem cells compared to a
#' tumour population with/without cells of others types.
#'
#' @param counts Expression matrix describing a cancer single-cell population; column
#'     names should be gene names, row names should be cancer cell names.
#' @param lengths Vector containing the lengths of all the genes or their transcripts.
#' @return Expression matrix containing Cells Per Million values.
#' @export
count_fractions_per_base_transcriptwise <- function(counts, lengths){
    
    if(ncol(counts) != length(lengths)){
        cat("ERROR: LENGTH OF COUNTS AND LENGTH OF LENGTHS NOT THE SAME\n")
    }
    count_fractions_per_base <- counts
    count_fractions_per_base[count_fractions_per_base > 0] <- 0

    for(i in 1:ncol(counts)){
        #if(i %% 1000 == 0){
        #    if(i == 1000){
        #        cat("Progress: ")
        #    }
        #    cat(paste0(i, "."))
        #}
        this_col_vec_of_counts <- counts[,i]
        this_col_vec_of_count_fractions_per_base <- vector()        
        if(sum(this_col_vec_of_counts) == 0){
            this_col_vec_of_count_fractions_per_base <- this_col_vec_of_counts
        } else {
            N <- sum(this_col_vec_of_counts)        
            for(j in 1:length(this_col_vec_of_counts)){
                this_count_fraction          <- (this_col_vec_of_counts[j] / N)
                this_count_fraction_per_base <- (this_count_fraction / (lengths[i]/1000) * 1e6)
                this_col_vec_of_count_fractions_per_base[j] <- this_count_fraction_per_base
            }        
        }
        count_fractions_per_base[,i] <- this_col_vec_of_count_fractions_per_base
    }
    #cat("..Completed\n")

    return(count_fractions_per_base)
}

#' Convert a matrix of read counts to a matrix of Transcripts Per Million values.
#'
#' This function converts an expression matrix of raw read counts, to
#' a matrix of TPM values.
#'
#' @param counts Expression matrix e.g. describing a cancer single-cell population where
#'     column names should be gene names, row names should be cancer cell names.
#' @param lengths Vector containing the lengths of all the transcripts.
#' @return Expression matrix containing TPM values.
#' @export
tpm_samplewise <- function(counts, lengths){

    if(ncol(counts) != length(lengths)){
        cat("ERROR: LENGTH OF COUNTS AND LENGTH OF LENGTHS NOT THE SAME\n")
    }

    tpms <- counts
    tpms[tpms > 0] <- 0

    for(i in 1:nrow(counts)){
        #if(i %% 1000 == 0){
        #    if(i == 1000){
        #        cat("Progress: ")
        #    }
        #    cat(paste0(i, "."))
        #}
        this_row_vec_of_counts <- counts[i,]
        this_row_vec_of_tpms <- vector()        

        if(sum(this_row_vec_of_counts) == 0){
            this_row_vec_of_tpms <- this_row_vec_of_counts
        } else {

            rates <- vector()
            
            for(j in 1:length(this_row_vec_of_counts)){
                this_rate <- -1
                this_rate <- (this_row_vec_of_counts[j] / lengths[j]) # by length
                rates[j] <- this_rate
            }
            
            #sum_rates            <- sum(rates)
            #sum_rates_tens      <- sum(rates) / 10 # same as with 1000
            sum_rates_thousands <- sum(rates) / 1000
            
            for(j in 1:length(rates)){
                #this_row_vec_of_tpms[j] <- (rates[j] / sum_rates)
                #this_row_vec_of_tpms[j] <- (rates[j] / sum_rates_tens)
                this_row_vec_of_tpms[j] <- 1e3 * (rates[j] / sum_rates_thousands)
                #this_row_vec_of_tpms[j] <- ((rates[j] / sum_rates) * 1e6)
            }

        }
        tpms[i,] <- this_row_vec_of_tpms
    }
    #cat("..Completed\n")

    return(tpms)
}
