#' Function to add zeros to the start of perplexity values to make them nicer for printing.
#'
#' This function takes a vector of numbers and adds zeros to the start of each to make
#'     them nicer for printing.
#'
#' @param Vector of numbers (perplexity values) to have zeros added to if they are <100. 
#' @return Vector of numbers with zeros now added.
#' @export
add_zeros_to_perplexity_values <- function(perplexityVals){

    perplexityVals_modified <- vector()

    for(i in 1:length(perplexityVals)){
        if(perplexityVals[i] < 10){
            perplexityVals_modified[i] <- paste0('00', perplexityVals[i])
        } else if(perplexityVals[i] < 100){
            perplexityVals_modified[i] <- paste0('0', perplexityVals[i])
        } else {
            perplexityVals_modified[i] <- perplexityVals[i]
        }
    }

    return(perplexityVals_modified)
}
