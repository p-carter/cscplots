#' Add scores to list of t-SNE results.
#'
#' This function adds scores to a list of tSNE results for later printing e.g. add stemness predictions,  
#' total transcription, total transcribed genes, etc. 
#'
#' @param tSNE_results_with_annotations List of data frames. Each data frame column 1 and 2 should be
#'     Y1 and Y2 t-SNE info respectively and additional columns containing t-SNE scores/annotations e.g.
#'     stemness score or total transcription; row names should be single-cell IDs.
#' @param perplexityValues List of numbers corresponding to the perplexity value used for each t-SNE
#'     in tSNE_results_with_annotations.
#' @param scores List of scores to add to tSNE_results_with_annotations.
#' @param scoresName Name to use for this score in tSNE_results_with_annotations.
#' @return List of dataframes, i.e. tSNE_results_with_annotations but with score now added as an additional
#'     field.
#' @export
add_scores_to_tSNE_results <- function(tSNE_results_with_annotations, perplexityValues, scores, scoresName){

    options(scipen=999)

    perplexityNames <- list()
    i <- 1
    for(p in perplexityValues){
        perplexityNames[[i]] <- paste0("Perplexity=", p)
        i <- i + 1
    }

    ## Add Scores
    for(i in 1:length(perplexityValues)){
        numberOfColumns <- length(tSNE_results_with_annotations[[i]])
        #tSNE_results_with_annotations_copy <- tSNE_results_with_annotations[[i]]
        #tSNE_results_with_annotations[[i]] <- data.frame(tSNE_results_with_annotations_copy, scores)
        tSNE_results_with_annotations[[i]] <- data.frame(tSNE_results_with_annotations[[i]], scores)

        colnames(tSNE_results_with_annotations[[i]])[(numberOfColumns + 1)] <- scoresName
    }

    return (tSNE_results_with_annotations)
}

