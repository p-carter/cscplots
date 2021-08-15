#' Run t-SNEs for cancer single-cells expression matrix for all perplexity values provided.
#'
#' This function runs one or more t-SNEs for a cancer single-cell expression matrix i.e. each 
#' using a different perplexity value, provided perplexityValues; recommended = 70. This
#' function calls runTSNE which generates the t-SNE.
#'
#' @param perplexityValues List of perplexity values to use e.g. c(5, 10, 20, 30, 40, 50, 60, 70);
#'     recommended = 70.
#' @param tExpressionMatrix Expression matrix for cancer single-cell population or sample; where
#'     column names should be gene names, and row names should be sample names.
#' @param runName Label for run, e.g. 'rawCounts_GSM4104152_SC4.LB17009'; can be anything.
#' @param outDir Directory name for printing to.
#' @param montage If TRUE use montage system command to combine previously generated plots in tif
#'     files.
#' @return A list of data frames, each containing a t-SNE result.
#' @export
run_all_tSNEs <- function(perplexityValues, tExpressionMatrix, runName, outDir, montage = FALSE){

    if(is.null(montage)){ montage <- FALSE }
    
    raw_tSNE_results <- list()
    filesToMontage <- ''

    i <- 1
    for(j in perplexityValues){
        cat(paste0("Perplexity = ", j, "\n"))
        if(j < 10){
            jString <- paste0("00", j)
        } else if(j < 100){
            jString <- paste0("0", j)
        } else {
            jString <- j
        }
        raw_tSNE_results[[i]] <- run_tSNE(j, tExpressionMatrix, runName, outDir)
        filesToMontage <- paste0(filesToMontage, " ", runName, "_tSNE-", jString, "-no_colours.tif ")
        i <- i + 1
    }

    perplexityNames <- list()
    i <- 1
    for(p in perplexityValues){
        perplexityNames[[i]] <- paste0("Perplexity=", p)
        i <- i + 1
    }
    names(raw_tSNE_results) <- perplexityNames

    if((montage == TRUE) && (length(perplexityValues > 1))){ ## Combine plots from multiple perplexity runs
        montageCommand <- paste0("cd ", outDir, "; montage -geometry +0+0 ", filesToMontage, " ", runName, "_tSNEs-no_colours-montage.tif")
        cat(paste0("montageCommand = ", montageCommand, "\n"))
        system(montageCommand)
        rmFilesCommand <- paste0("cd ", outDir, "; rm ", filesToMontage)
        system(rmFilesCommand)
    }
    
    return(raw_tSNE_results)
}

#' Run a t-SNE for cancer single-cells expression matrix for a perplexity value provided.
#'
#' This function runs a t-SNE for a cancer single-cell expression matrix using a provided perplexity 
#'     value e.g. 70.
#'
#' @param perplexityValue Perplexity value to use e.g. 5-130; recommended = 70.
#' @param tExprMatr Expression matrix for cancer single-cell population; where
#'     column names should be gene names, and row names should be sample names.
#' @param runName Label for run, e.g. 'rawCounts_GSM4104152_SC4.LB17009'; can be anything.
#' @param outDir Directory name for printing to.
#' @return A data frame describing the t-SNE result i.e. Y1 & Y2.
#' @export
run_tSNE <- function(perplexityValue, tExprMatr, runName, outDir){

    #library("ggplot2")
    #library("Rtsne")

    if(perplexityValue < 10){
        perplexityValueString <- paste0("00", perplexityValue)
    } else if(perplexityValue < 100){
        perplexityValueString <- paste0("0", perplexityValue)
    } else {
        perplexityValueString <- perplexityValue
    }

    if(nrow(tExprMatr) < 1000){
        point_size <- 0.3
    } else if((nrow(tExprMatr) >= 1000) && (nrow(tExprMatr) <= 10000)){
        point_size <- 0.2
    } else {
        point_size <- 0.1
    }
    
    tSNE_out <- Rtsne::Rtsne(as.matrix(tExprMatr), perplexity = perplexityValue, max_iter = 1000)
    tSNE_DF <- data.frame(tSNE_out$Y[,1], tSNE_out$Y[,2])
    colnames(tSNE_DF) <- c("Y1", "Y2")
    rownames(tSNE_DF) <- rownames(tExprMatr)

    head(tSNE_DF)
    outFile <- paste0(outDir, runName, "_tSNE-", perplexityValueString, "-no_colours.tif")
    p <- ggplot2::ggplot(tSNE_DF, aes(x = Y1, y = Y2)) +
    geom_point(size = point_size, alpha = 0.7) + 
    theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.spacing = unit(0.1, "cm"),
          plot.title = element_text(size = 9)
    )
    plot_title <- paste0("t-SNE: ", runName, "; p", perplexityValue)

    p <- p + ggtitle(plot_title)
    p <- p + coord_fixed(ratio = 1/1)

    print(p)      
    ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 200, width = 5, height = 4.5, units = "in")

    return(tSNE_DF)
}
