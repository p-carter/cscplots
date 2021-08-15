#' Get correlations between a single-cell expression matrix and a vector of scores for each cell.
#'
#' This function uses default correlation (Pearson) to give correlation values between a single-cell
#'     RNA-Seq expression matrix and a vector of scores which contains one score for each cell.
#'     For example, this could be a vector of predicted stemness scores, or a vector of
#'     transcription summary values e.g. total transcribed genes, or total transcription for each cell.
#'
#' @param tExprMatrix Data frame describing t-SNE or PCA and score/s including the score to print.
#' @param scores Vector of scores to correlate with each row of tExprMatrix i.e. the expression
#'     levels for all genes for each single-cell.
#' @return List containing two named vectors that respectively contain the positive & negative
#'     correlation results.
#' @export
get_expression_scores_correlations <- function(tExprMatrix, scores){
    
    sample_size <- nrow(tExprMatrix)
    score_order <- order(scores)
    dfX <- tExprMatrix[score_order,]
    score_corr_results <- apply(dfX, 2, function(col) suppressWarnings(cor(col, 1:sample_size))) # correlate each column (gene) ordered by score with 1:nrows
    
    positive_threshold <- 0.01
    negative_threshold <- -0.01
    
    score_corr_results_positive       <- score_corr_results[which(score_corr_results > positive_threshold)]
    score_corr_results_positive_order <- rev(order(score_corr_results[which(score_corr_results > positive_threshold)]))
    sc_positive_ordered               <- score_corr_results_positive[score_corr_results_positive_order]
    
    score_corr_results_negative       <- score_corr_results[which(score_corr_results < negative_threshold)]
    score_corr_results_negative_order <- order(score_corr_results[which(score_corr_results < negative_threshold)])
    sc_negative_ordered               <- score_corr_results_negative[score_corr_results_negative_order]
    
    return(list(sc_positive_ordered, sc_negative_ordered))
}
