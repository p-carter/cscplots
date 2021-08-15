#' Plot a score overlaid on a supplied t-SNE or PCA generated using a cancer
#'     single-cell population expression matrix.
#'
#' This function is a generic function for plotting a t-SNE or PCA previously generated using 
#' cancer single-cell population expression matrix. It can be used to label or colour 
#' individual or small groups of cells to assist identification of subpopulations of interest, 
#' e.g. groups of cancer stem cells.
#'
#' @param myData Data frame describing t-SNE or PCA and score/s including the score to print.
#' @param xAxisName Name of column in myData that describes x-axis info for plot e.g. "Y1" in a t-SNE.
#' @param yAxisName Name of column in myData that describes y-axis info for plot e.g. "Y2" in a t-SNE.
#' @param runName Run name, e.g. 'tSNE_070'; will be used in the title of the plot and the output tif name.
#' @param scoresName Name of column in myData that contains the score info for plotting e.g. 'ranger_hipsci_raw_IPS'.
#' @param datasetName Name of dataset, e.g. 'mela_4558_TPM'.
#' @param outDir Output directory to use.
#' @param labelPositions Optional: Row number/s of myData i.e. cell positions, to print cell labels for.
#' @param colourPositions Optional: Row number/s of myData i.e. cell positions, to highlight cells (in black).
#' @param highestOnTop Optional: Set to TRUE to print highest scoring cells on top of lower scoring cells.
#' @return No return value.
#' @export
plot_scores <- function(myData, xAxisName, yAxisName, runName, scoresName, datasetName, outDir, labelPositions = NULL, colourPositions = NULL,
                        highestOnTop = NULL){

    #library("ggplot2")

    if(is.null(highestOnTop)){ highestOnTop <- 0 }
    #cat(paste0("Highest On Top = ", highestOnTop, "\n"))
    if(highestOnTop == TRUE){
        scoreColNumber <- which(colnames(myData) %in% scoresName)
        new_row_order_for_highest_on_top <- order(myData[,scoreColNumber])
        myData <- myData[new_row_order_for_highest_on_top,] # reorder
        if(length(labelPositions) > 0){
            labelPositions_for_highest_on_top_order <- match(labelPositions, new_row_order_for_highest_on_top)
            labelPositions <- labelPositions_for_highest_on_top_order
        }
        if(length(colourPositions) > 0){
            colourPositions_for_highest_on_top_order <- match(colourPositions, new_row_order_for_highest_on_top)
            colourPositions <- colourPositions_for_highest_on_top_order
        }
        outFile <- paste0(outDir, runName, "-", scoresName, "_", xAxisName, "vs", yAxisName, "-highestOnTop.tif")
    } else {
        outFile <- paste0(outDir, runName, "-", scoresName, "_", xAxisName, "vs", yAxisName, ".tif")
    }

    printScoresName <- scoresName
    if(length(grep("-", scoresName)) >= 1){
       scoreColNumber <- which(colnames(myData) %in% scoresName)
       scoresName <- gsub("-", "_", scoresName)
       colnames(myData)[scoreColNumber] <- scoresName
    }

    #point_size <- 0.7
    if(nrow(myData) < 100){
        point_size <- 0.7
    } else if((nrow(myData) >= 100) && (nrow(myData) < 1000)){
        point_size <- 0.7
    } else if((nrow(myData) >= 1000) && (nrow(myData) < 10000)){
        point_size <- 0.7
    } else {
        point_size <- 0.5
    }

    p <- ggplot2::ggplot(myData, aes_string(x=xAxisName, y=yAxisName)) +
    geom_point(size=point_size, alpha=0.7, aes_string(colour=scoresName)) + 
    theme(
          axis.line = element_line(colour = "black", size = 1), 
          axis.title.x = element_text(size = 12), 
          axis.title.y = element_text(size = 12), 
          axis.text.x = element_text(size = 10), 
          axis.text.y = element_text(size = 10), 
          axis.ticks.length = unit(.2, "cm"), 
          legend.text = element_text(size = 8),
          text = element_text(size = 12),
          legend.box = "horizontal",
#          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.spacing = unit(0.1, "cm"),
    ) + scale_colour_gradientn(colours=rev(c("red", "orange", "yellow", "green", "cyan", "blue"))) 
    plot_title <- paste0(datasetName, "; ", runName)

    if(nchar(scoresName) > 10){
        scoresName_label <- gsub("(.{10,}?)_", "\\1\n", printScoresName) # insert some new lines if long string
    } else {
        scoresName_label <- printScoresName
    }

    p <- p + labs(colour=scoresName_label, size="", fill="", shape="") + ggtitle(plot_title)
    p <- p + coord_fixed(ratio = 1/1)

    min_x <- min(round(myData[[xAxisName]]))
    max_x <- max(round(myData[[xAxisName]]))
    min_y <- min(round(myData[[yAxisName]]))
    max_y <- max(round(myData[[yAxisName]]))

    min_x_break <- round_down(min_x, by <- 10)
    max_x_break <- round_up(max_x, 10)
    min_y_break <- round_down(min_y, 10)
    max_y_break <- round_up(max_y, 10)
        
    breaks_x <- c(seq(min_x_break, max_x_break, by = 10))
    breaks_y <- c(seq(min_y_break, max_y_break, by = 10))

    p <- p + scale_x_continuous(breaks = breaks_x)
    p <- p + scale_y_continuous(breaks = breaks_y)
    
    # add some cell labels if specified
    if(length(labelPositions) > 0){
        #library("ggrepel")
        dataSubset <- myData[labelPositions,]
        subsetNames <- rownames(dataSubset)
        options(ggrepel.max.overlaps = Inf)
        p <- p + ggrepel::geom_text_repel(data = dataSubset, aes(label = subsetNames), size = 1.6,
                                          segment.color = 'black', segment.size = 0.1, 
                                          force = 1 #2
                                          )
    }
    # add some black points if specified
    if(length(colourPositions) > 0){
        dataSubset <- myData[colourPositions,]
        p <- p + geom_point(data = dataSubset, colour = "black", size = point_size, alpha = 1)
    }

    
    
    print(p)
    ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 300, width = 15, height = 13.5, units = "in", limitsize = FALSE)

    return()
}


#' Plot a score overlaid on a supplied t-SNE or PCA generated using a cancer
#'     single-cell population expression matrix, with grid in background.
#'
#' This function is a generic function for plotting a t-SNE or PCA previously generated using 
#' cancer single-cell population expression matrix. The grid in the background can help when
#' looking at exact positions of points (cells). It can also be used to label or colour 
#' individual or small groups of cells to assist identification of subpopulations of interest 
#' e.g. groups of cancer stem cells.
#'
#' @param myData Data frame describing t-SNE or PCA and additional info including the score to print.
#' @param xAxisName Name of column in myData that describes x-axis info for plot e.g. "Y1" in a t-SNE.
#' @param yAxisName Name of column in myData that describes y-axis info for plot e.g. "Y2" in a t-SNE.
#' @param runName Run name, e.g. 'tSNE_070'; will be used in the title of the plot and the output tif name.
#' @param scoresName Name of column in myData that contains the score info for plotting e.g. 'ranger_hipsci_raw_IPS'.
#' @param datasetName Name of dataset e.g. 'mela_4558_TPM'.
#' @param outDir Output directory to use.
#' @param labelPositions Optional: Row number/s of myData i.e. cell positions, to print cell labels for.
#' @param colourPositions Optional: Row number/s of myData i.e. cell positions, to highlight cells (in black).
#' @param highestOnTop Optional: Set to TRUE to print highest scoring cells on top of lower scoring cells.
#' @return No return value.
#' @export
plot_scores_with_grid <- function(myData, xAxisName, yAxisName, runName, scoresName, datasetName, outDir, labelPositions = NULL, colourPositions = NULL,
                                  highestOnTop = NULL){

    #library("ggplot2")

    if(is.null(highestOnTop)){
        highestOnTop <- FALSE
    }
    
    if(highestOnTop == TRUE){
        scoreColNumber <- which(colnames(myData) %in% scoresName)
        new_row_order_for_highest_on_top <- order(myData[,scoreColNumber])
        myData <- myData[new_row_order_for_highest_on_top,] # reorder
        if(length(labelPositions) > 0){
            labelPositions_for_highest_on_top_order <- match(labelPositions, new_row_order_for_highest_on_top)
            labelPositions <- labelPositions_for_highest_on_top_order
        }
        if(length(colourPositions) > 0){
            colourPositions_for_highest_on_top_order <- match(colourPositions, new_row_order_for_highest_on_top)
            colourPositions <- colourPositions_for_highest_on_top_order
        }
        outFile <- paste0(outDir, runName, "-", scoresName, "_", xAxisName, "vs", yAxisName, "-highestOnTop-with_grid.tif")
    } else {
        outFile <- paste0(outDir, runName, "-", scoresName, "_", xAxisName, "vs", yAxisName, "-with_grid.tif")
    }

    if((length(labelPositions) > 0) && (length(colourPositions) > 0)){
        outFile <- gsub(".tif", "-selected_cells_and_annotated.tif", outFile)
    } else if(length(colourPositions) > 0){
        outFile <- gsub(".tif", "-selected_cells.tif", outFile)
    } else if(length(labelPositions) > 0){
        outFile <- gsub(".tif", "-annotated.tif", outFile)
    }
    
    printScoresName <- scoresName
    if(length(grep("-", scoresName)) >= 1 ){
       scoreColNumber <- which(colnames(myData) %in% scoresName)
       scoresName <- gsub("-", "_", scoresName)
       colnames(myData)[scoreColNumber] <- scoresNameNew
    }

    #point_size <- 0.1
    if(nrow(myData) < 100){
        point_size <- 0.5
    } else if((nrow(myData) >= 100) && (nrow(myData) < 1000)){
        point_size <- 0.4
    } else {
        point_size <- 0.3
    }

    p <- ggplot2::ggplot(myData, aes_string(x=xAxisName, y=yAxisName)) +

    geom_point(size = point_size, alpha = 0.7, aes_string(colour = scoresName)) + 

#    theme( axis.line = element_line(colour = "black") ) +
    theme(
          axis.line = element_line(colour = "black", size = 1), 
          axis.title.x = element_text(size = 12), 
          axis.title.y = element_text(size = 12), 
          axis.text.x = element_text(size = 7), 
          axis.text.y = element_text(size = 7), 
          axis.ticks.length = unit(.2, "cm"),
          legend.text = element_text(size = 8),
          text = element_text(size = 12),
          legend.box = "horizontal"
         ) + 

    scale_colour_gradientn(colours=rev(c("red", "orange", "yellow", "green", "cyan", "blue")))
    #scale_x_continuous(breaks = round(seq(min(myData[[xAxisName]]-2), max(myData[[xAxisName]]+2), by = 1), 0) ) +
    #scale_y_continuous(breaks = round(seq(min(myData[[yAxisName]]-2), max(myData[[yAxisName]]+2), by = 1), 0) ) 

    min_x <- min(round(myData[[xAxisName]]))
    max_x <- max(round(myData[[xAxisName]]))
    min_y <- min(round(myData[[yAxisName]]))
    max_y <- max(round(myData[[yAxisName]]))

    min_x_break <- round_down(min_x, by <- 10)
    max_x_break <- round_up(max_x, 10)
    min_y_break <- round_down(min_y, 10)
    max_y_break <- round_up(max_y, 10)
        
    breaks_x <- c(seq((min_x_break - 5), (max_x_break + 5), by = 1))
    breaks_y <- c(seq((min_y_break - 5), (max_y_break + 5), by = 1))

    p <- p + scale_x_continuous(breaks = breaks_x)
    p <- p + scale_y_continuous(breaks = breaks_y)
    
    plot_title <- paste0(datasetName, "; ", runName)
    if(nchar(scoresName) > 10){
        scoresName_label <- gsub("(.{10,}?)_", "\\1\n", printScoresName) # insert some new lines if long string
    } else {
        scoresName_label <- printScoresName
    }

    p <- p + labs(colour=scoresName_label, size="", fill="", shape="") + ggtitle(plot_title)
    p <- p + coord_fixed(ratio = 1/1)

    # add some cell labels if specified
    if(length(labelPositions) > 0){
        #library("ggrepel")
        
        dataSubset <- myData[labelPositions,]
        subsetNames <- rownames(dataSubset)
        options(ggrepel.max.overlaps = Inf)
        p <- p + ggrepel::geom_text_repel(data = dataSubset, aes(label = subsetNames), size = 1.6,
                                          segment.color = 'black', segment.size = 0.1, 
                                          force = 2)
    }
    # add some small black points to cells if specified
    if(length(colourPositions) > 0){
        dataSubset <- myData[colourPositions,]
        p <- p + geom_point(data = dataSubset, colour = "black", size = 0.1, alpha = 0.7)
    }

    print(p)
    ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 300, width = 15, height = 13.5, units = "in", limitsize = FALSE)

    return()
}


#' Plot a cell category overlaid on a supplied t-SNE or PCA generated using a cancer
#'     single-cell population expression matrix.
#'
#' This function is a generic function for plotting a t-SNE or PCA previously generated using 
#' cancer single-cell population expression matrix. It can be used to label or colour 
#' individual or groups of cells by categories of information (fields containing categorical variables/factors),
#' e.g. to assist identification of subpopulations of interest such as groups of cancer stem cells.
#'
#' @param myData Data frame describing t-SNE or PCA and additional info including the category to print.
#' @param xAxisName Name of column in myData that describes x-axis info for plot e.g. "Y1" in a t-SNE.
#' @param yAxisName Name of column in myData that describes y-axis info for plot e.g. "Y2" in a t-SNE.
#' @param runName Run name, e.g. 'tSNE_070'; will be used in the title of the plot and the output tif name.
#' @param categoryName Name of column in myData that contains the category info for plotting e.g. 'cell_type'.
#' @param datasetName Name of dataset e.g. 'mela_4558_TPM'.
#' @param outDir Output directory to use.
#' @param labelPositions Optional: Row number/s of myData i.e. cell positions, to print cell labels for.
#' @param colourPositions Optional: Row number/s of myData i.e. cell positions, to highlight cells (in black).
#' @param labelsVec Optional: Vector containing custom labels to display in legend for the printed category.
#' @param labelField Optional: If labelPositions parameter has been provided, labelField can provide the name of
#'     another categorical variable within a column of myData to label the defined labelPositions.
#' @return No return value.
#' @export
plot_categories <- function(myData, xAxisLabel, yAxisLabel, runName, categoryName, datasetName, outDir, labelPositions = NULL,
                            colourPositions = NULL, labelsVec = NULL, labelField = NULL){

    #library("ggplot2")

    if(!is.factor(myData[[categoryName]])){
        if(length(unique(myData[[categoryName]])) > 50){ return() } # probably a continuous field
        #cat(paste0("plot_categories is for printing factor fields; converting but not appropriate for continuous fields"))
        myData[[categoryName]] <- as.factor(myData[[categoryName]])
    }

    if(is.numeric(myData[[categoryName]])){
        myData[[categoryName]] <- factor(myData[[categoryName]], levels = sort(as.numeric(levels(myData[[categoryName]]))))
    }    
    if(is.null(labelsVec)){
        labelsVec <- levels(myData[[categoryName]])
    }
    
    numberOfLevels <- length(levels(myData[[categoryName]]))
    numberOfFieldsPresentInThisData <- length(unique(myData[[categoryName]]))

    if(numberOfLevels != numberOfFieldsPresentInThisData){
        cat(paste0("This subset doesn't contain all of the fields, so the labelsVec will be adjusted.\n"))

        fields_present_numbers <- as.integer(sort(unique(myData[[categoryName]])))
        labelsVec <- labelsVec[fields_present_numbers]
    }

    rainbowTmp <- rainbow(length(unique(myData[[categoryName]])))

    palette(rainbowTmp)
    rainbow_palette <- palette(rainbowTmp)

    if(nrow(myData) < 100){
        point_size <- 0.5
    } else if((nrow(myData) >= 100) && (nrow(myData) < 1000)){
        point_size <- 0.4
    } else {
        point_size <- 0.3
    }
    
    #    outFile <- paste0(outDir, plotTitle, "-Categories_", categoryName, "_", xAxisLabel, "vs", yAxisLabel, ".tif")
    outFile <- paste0(outDir, runName, "-", categoryName, "_", xAxisLabel, "vs", yAxisLabel, ".tif")
    
    p <- ggplot2::ggplot(myData, aes_string(x = xAxisLabel, y = yAxisLabel), label = rownames(myData))
    p <- p + geom_point(size = point_size, alpha = 0.7, aes_string(colour = categoryName))
    p <- p + theme(
                axis.line = element_line(colour = "black", size = 1), 
                axis.title.x = element_text(size = 12), 
                axis.title.y = element_text(size = 12), 
                axis.text.x = element_text(size = 10), 
                axis.text.y = element_text(size = 10), 
#                axis.line = element_line(colour = "black", size = 1.5), 
#                axis.title.x = element_text(size = 15), 
#                axis.title.y = element_text(size = 15), 
#                axis.text.x = element_text(size = 12), 
#                axis.text.y = element_text(size = 12), 
                axis.ticks.length = unit(.2, "cm"), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.spacing = unit(0.1, "cm"),
                legend.text = element_text(size = 8),
                text = element_text(size = 12),
                legend.box = "horizontal"
    ) 
    #print(rainbow_palette)
    #rainbow_palette[1] <- "lightgray"
    #print(rainbow_palette)

    p <- p + scale_color_manual(name = categoryName, 
                                labels = labelsVec,
                                values = rainbow_palette
                                )
    #    this_plot_title <- paste0(plotTitle, " with ", categoryName)
    plot_title <- paste0(datasetName, '; ', runName, ": ", categoryName) 

    if(nchar(categoryName) > 10){
        categoryName_label <- gsub("(.{10,}?)_", "\\1\n", categoryName) # insert some new lines if long string
    } else {
        categoryName_label <- categoryName
    }
    
    p <- p + labs(colour = categoryName_label, size = "", fill = "", shape = "") + ggtitle(plot_title)
    p <- p + coord_fixed(ratio = 1/1)

    min_x <- min(round(myData[[xAxisLabel]]))
    max_x <- max(round(myData[[xAxisLabel]]))
    min_y <- min(round(myData[[yAxisLabel]]))
    max_y <- max(round(myData[[yAxisLabel]]))

    min_x_break <- round_down(min_x, by <- 10)
    max_x_break <- round_up(max_x, 10)
    min_y_break <- round_down(min_y, 10)
    max_y_break <- round_up(max_y, 10)
        
    breaks_x <- c(seq(min_x_break, max_x_break, by = 10))
    breaks_y <- c(seq(min_y_break, max_y_break, by = 10))

    p <- p + scale_x_continuous(breaks = breaks_x)
    p <- p + scale_y_continuous(breaks = breaks_y)

    #if(!is.null(labelPositions)){
        if(length(labelPositions) > 0){
            #library("ggrepel")
            
            dataSubset <- myData[labelPositions,]
        
            if(!is.null(labelField)){
                subsetNames <- dataSubset[[labelField]]
            } else {
                subsetNames <- rownames(dataSubset)
            }
            options(ggrepel.max.overlaps = Inf)
            p <- p + #geom_text(data = dataSubset, size = 0.8, aes(label = subsetNames))
                     ggrepel::geom_text_repel(data = dataSubset, aes(label = subsetNames), size = 1.6,
                                              segment.color = 'black', segment.size = 0.1, 
                                              force = 1 #2
                                             )
        }
    #}
    #if(!is.null(colourPositions)){    
        if(length(colourPositions) > 0){
            dataSubset <- myData[colourPositions,]
            #p <- p + geom_point(data = dataSubset, colour = "green", alpha = 0.5)
            p <- p + geom_point(data = dataSubset, colour = "black", size = point_size, alpha = 1)
        }
    #}

    print(p)
    ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 300, width = 15, height = 13.5, units = "in")

    return()
}


#' Plot a cell category overlaid on a supplied t-SNE or PCA generated using a cancer
#'     single-cell population expression matrix, with grid in background.
#'
#' This function is a generic function for plotting a t-SNE or PCA previously generated using 
#' cancer single-cell population expression matrix. It can be used to label or colour 
#' individual or groups of cells by categories of information (fields containing categorical variables/factors),
#' e.g. to assist identification of subpopulations of interest such as groups of cancer stem cells.
#'
#' @param myData Data frame describing t-SNE or PCA and additional info including the category to print.
#' @param xAxisName Name of column in myData that describes x-axis info for plot e.g. "Y1" in a t-SNE.
#' @param yAxisName Name of column in myData that describes y-axis info for plot e.g. "Y2" in a t-SNE.
#' @param runName Run name, e.g. 'tSNE_070'; will be used in the title of the plot and the output tif name.
#' @param categoryName Name of column in myData that contains the category info for plotting e.g. 'cell_type'.
#' @param datasetName Name of dataset e.g. 'mela_4558_TPM'.
#' @param outDir Output directory to use.
#' @param labelPositions Optional: Row number/s of myData i.e. cell positions, to print cell labels for.
#' @param colourPositions Optional: Row number/s of myData i.e. cell positions, to highlight cells (in black).
#' @param labelsVec Optional: Vector containing custom labels to display in legend for the printed category.
#' @param labelField Optional: If labelPositions parameter has been provided, labelField can provide the name of
#'     another categorical variable within a column of myData to label the defined labelPositions.
#' @return No return value.
#' @export
plot_categories_with_grid <- function(myData, xAxisLabel, yAxisLabel, runName, categoryName, datasetName, outDir, labelPositions = NULL,
                                      colourPositions = NULL, labelsVec = NULL, labelField = NULL){

    #library("ggplot2")

    if(!is.factor(myData[[categoryName]])){
        if(length(unique(myData[[categoryName]])) > 50){ return() } # probably a continuous field
        #cat(paste0("plot_categories is for printing factor fields; converting but not appropriate for continuous fields"))
        myData[[categoryName]] <- as.factor(myData[[categoryName]])
    }

    if(is.numeric(myData[[categoryName]])){
        myData[[categoryName]] <- factor(myData[[categoryName]], levels = sort(as.numeric(levels(myData[[categoryName]]))))
    }    
    if(is.null(labelsVec)){
        labelsVec <- levels(myData[[categoryName]])
    }
    
    numberOfLevels <- length(levels(myData[[categoryName]]))
    numberOfFieldsPresentInThisData <- length(unique(myData[[categoryName]]))

    if(numberOfLevels != numberOfFieldsPresentInThisData){
        cat(paste0("This subset doesn't contain all of the fields, so the labelsVec will be adjusted.\n"))

        fields_present_numbers <- as.integer(sort(unique(myData[[categoryName]])))
        labelsVec <- labelsVec[fields_present_numbers]
    }

    rainbowTmp <- rainbow(length(unique(myData[[categoryName]])))

    palette(rainbowTmp)
    rainbow_palette <- palette(rainbowTmp)

    #point_size <- 0.1
    if(nrow(myData) < 100){
        point_size <- 0.5
    } else if((nrow(myData) >= 100) && (nrow(myData) < 1000)){
        point_size <- 0.4
    } else {
        point_size <- 0.3
    }
    
    #outFile <- paste0(outDir, runName, "-", categoryName, "_", xAxisLabel, "vs", yAxisLabel, ".tif")
    outFile <- paste0(outDir, runName, "-", categoryName, "_", xAxisLabel, "vs", yAxisLabel, "-with_grid.tif")
    if((length(labelPositions) > 0) && (length(colourPositions) > 0)){
        outFile <- gsub(".tif", "-selected_cells_and_annotated.tif", outFile)
    } else if(length(colourPositions) > 0){
        outFile <- gsub(".tif", "-selected_cells.tif", outFile)
    } else if(length(labelPositions) > 0){
        outFile <- gsub(".tif", "-annotated.tif", outFile)
    }
    
    p <- ggplot2::ggplot(myData, aes_string(x = xAxisLabel, y = yAxisLabel), label = rownames(myData))
    p <- p + geom_point(size = point_size, alpha = 0.7, aes_string(colour = categoryName))
    p <- p + theme(axis.line = element_line(colour = "black", size = 1), 
                   axis.title.x = element_text(size = 12), 
                   axis.title.y = element_text(size = 12), 
                   axis.text.x = element_text(size = 7), 
                   axis.text.y = element_text(size = 7), 
                   axis.ticks.length = unit(.2, "cm"), 
                   panel.spacing = unit(0.1, "cm"),
                   legend.text = element_text(size = 8),
                   text = element_text(size = 12),
                   legend.box = "horizontal")

    #print(rainbow_palette)
    #rainbow_palette[1] <- "lightgray"
    #print(rainbow_palette)
    p <- p + scale_color_manual(name = categoryName, 
                                labels = labelsVec,
                                values = rainbow_palette) #+ 

    #scale_x_continuous(breaks = round(seq(min(myData[[xAxisLabel]]-2), max(myData[[xAxisLabel]]+2), by = 1), 0) ) +
    #scale_y_continuous(breaks = round(seq(min(myData[[yAxisLabel]]-2), max(myData[[yAxisLabel]]+2), by = 1), 0) ) 

    min_x <- min(round(myData[[xAxisLabel]]))
    max_x <- max(round(myData[[xAxisLabel]]))
    min_y <- min(round(myData[[yAxisLabel]]))
    max_y <- max(round(myData[[yAxisLabel]]))

    min_x_break <- round_down(min_x, by <- 10)
    max_x_break <- round_up(max_x, 10)
    min_y_break <- round_down(min_y, 10)
    max_y_break <- round_up(max_y, 10)
        
    breaks_x <- c(seq((min_x_break - 5), (max_x_break + 5), by = 1))
    breaks_y <- c(seq((min_y_break - 5), (max_y_break + 5), by = 1))

    p <- p + scale_x_continuous(breaks = breaks_x)
    p <- p + scale_y_continuous(breaks = breaks_y)
    
    #this_plot_title <- paste0(plotTitle, " with ", categoryName)    
    plot_title <- paste0(datasetName, '; ', runName, ": ", categoryName) 
    p <- p + labs(colour = categoryName, size = "", fill = "", shape = "") + ggtitle(plot_title)
    p <- p + coord_fixed(ratio = 1/1)

    if(length(labelPositions) > 0){
        #library("ggrepel")
        
        dataSubset <- myData[labelPositions,]
    
        if(!is.null(labelField)){
            subsetNames <- dataSubset[[labelField]]
        } else {
            subsetNames <- rownames(dataSubset)
        }
        options(ggrepel.max.overlaps = Inf)
        p <- p + #geom_text(data = dataSubset, size = 0.8, aes(label = subsetNames))
                 ggrepel::geom_text_repel(data = dataSubset, aes(label = subsetNames), size = 1.6,
                                          segment.color = 'black', segment.size = 0.1, 
                                          force = 1) #2
    }
    if(length(colourPositions) > 0){
        dataSubset <- myData[colourPositions,]
        #p <- p + geom_point(data = dataSubset, colour = "green", alpha = 0.5)
        p <- p + geom_point(data = dataSubset, colour = "black", size = point_size, alpha = 1)
    }

    print(p)
    ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 300, width = 15, height = 13.5, units = "in")

    return()
}


#' Plot a score overlaid on a supplied t-SNE or PCA generated using a cancer single-cell population
#' expression matrix; a range in the score defined by 2 thresholds is coloured grey and light grey below
#' the lower threshold.
#'
#' For printing tSNEs or PCAs with additional field (e.g. individual gene) specified by
#' fieldName overlaid with an additional colour range in grey specified by thresholdOne &
#' thresholdTwo parameters (e.g. useful for lower gene expression range such as 0 to 1).
#' In grey >thresholdOne && < thresholdTwo; lgrey <= thresholdOne; coloured >= thresholdTwo.
#' It can be also be used to label or colour individual or small groups of cells to assist
#' identification of subpopulations of interest e.g. groups of cancer stem cells.
#'
#' @param myData Data frame describing t-SNE or PCA and score/s including the score to print.
#' @param xAxisName Name of column in myData that describes x-axis info for plot e.g. "Y1" in a t-SNE.
#' @param yAxisName Name of column in myData that describes y-axis info for plot e.g. "Y2" in a t-SNE.
#' @param fieldName Name of column in myData that contains the field for colouring the plot with.
#' @param fullFieldName Full name of fieldName to use in the plot title e.g. full HUGO gene name.
#' @param thresholdOne Lower value to print cells in grey. Cells <= this are in light grey.
#' @param thresholdTwo Upper value to print cells in grey. Cells > thresholdOne && < thresholdTwo are grey.
#' @param runName Run name, e.g. 'tSNE_070'; will be used in the title of the plot and the output tif name.
#' @param outDir Output directory to use.
#' @param labelPositions Optional: Row number/s of myData i.e. cell positions, to print cell labels for.
#' @param colourPositions Optional: Row number/s of myData i.e. cell positions, to highlight cells (in black).
#' @param labelField Optional: If labelPositions parameter has been provided, labelField can provide the name of
#'     another categorical variable within a column of myData to label the defined labelPositions. E.g.
#'     use labelPositions to define malignant cell clusters to be labelled, and use labelField to
#'     label those cells with additional clinical info if available such as patient age, sex, or treament/s given. 
#' @param highestOnTop Optional: Set to TRUE to print highest scoring cells on top of lower scoring cells.
#' @param point_size Optional: Size of points in plot. 
#' @return No return value.
#' @export
plot_field_with_thresholds <- function(myData, xAxisName, yAxisName, fieldName, fullFieldName, thresholdOne, thresholdTwo, runName, outDir, 
                                       labelPositions = NULL, colourPositions = NULL, labelField = NULL, highestOnTop = NULL, point_size = NULL){

    #library("ggplot2")

    if((is.null(thresholdOne)) || (is.null(thresholdTwo))){
        cat(paste0("The code should be better but one or both thresholds were not found\n"))
        return()
    }
    if(is.null(point_size)){
        point_size <- 0.2
    }
    if(is.null(highestOnTop)){
        highestOnTop <- FALSE
    }
    if(thresholdOne > thresholdTwo){ # thresholdOne should be the lower threshold
        thresholdOne_tmp <- thresholdOne
        thresholdTwo_tmp <- thresholdTwo
        thresholdOne <- thresholdTwo_tmp
        thresholdTwo <- thresholdOne_tmp         
    }

    scoresName <- fieldName
    printFieldName <- fieldName
    if(length(grep("-", scoresName)) >= 1 ){
       scoreColNumber <- which(colnames(myData) %in% scoresName)
       scoresName  <- gsub("-", "_", scoresName)
       colnames(myData)[scoreColNumber] <- scoresName
       fieldName <- scoresName
    }

    if(highestOnTop == TRUE){
        scoreColNumber <- which(colnames(myData) %in% scoresName)
        new_row_order_for_highest_on_top <- order(myData[,scoreColNumber])
        myData <- myData[new_row_order_for_highest_on_top,] # reorder
        if(length(labelPositions) > 0){
            labelPositions_for_highest_on_top_order <- match(labelPositions, new_row_order_for_highest_on_top)
            labelPositions <- labelPositions_for_highest_on_top_order
        }
        if(length(colourPositions) > 0){
            colourPositions_for_highest_on_top_order <- match(colourPositions, new_row_order_for_highest_on_top)
            colourPositions <- colourPositions_for_highest_on_top_order
        }
        outFile <- paste0(outDir, runName, "-", printFieldName, "-with_thresholds_", xAxisName, "vs", yAxisName, "-highestOnTop.tif")
    } else {
	outFile <- paste0(outDir, runName, "-", printFieldName, "-with_thresholds_", xAxisName, "vs", yAxisName, ".tif")
    }

    over_thresholds_subsetPositions <- which(myData[[fieldName]] >= thresholdTwo)
    over_thresholds_dataSubset <- myData[over_thresholds_subsetPositions,]

    if(length(unique(over_thresholds_dataSubset[[fieldName]])) == 1){
        only_one_unique_val_above_threshold <- 1
    } else {
        only_one_unique_val_above_threshold <- 0
    }
    under_thresholds_subsetPositions <- which(myData[[fieldName]] <= thresholdOne)
    under_thresholds_dataSubset <- myData[under_thresholds_subsetPositions,]

    aboveThresholdOne_Positions <- which(myData[[fieldName]] > thresholdOne)
    belowThresholdTwo_Positions <- which(myData[[fieldName]] < thresholdTwo)

    within_thresholds_subsetPositions <- aboveThresholdOne_Positions[which(aboveThresholdOne_Positions %in% belowThresholdTwo_Positions)]

    within_thresholds_dataSubset <- myData[within_thresholds_subsetPositions,]
    within_thresholds_subsetNames <- rownames(within_thresholds_dataSubset)

    # Print the plot
    p <- ggplot2::ggplot(myData, aes_string(x=xAxisName, y=yAxisName)) +

    geom_point(data = under_thresholds_dataSubset, size = point_size, alpha = 0.7, colour = "grey92") + 
    geom_point(data = within_thresholds_dataSubset, size = point_size, alpha = 0.7, colour = "grey80") + 
    geom_point(data = over_thresholds_dataSubset, size = point_size, alpha = 0.7, aes_string(colour = fieldName)) + 

    theme(
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.spacing = unit(0.1, "cm") #,
          #legend.text = element_text(size = 1)
    ) 

    if(only_one_unique_val_above_threshold == 1){ # if only one value above the threshold range, just print it as blue as easier to see
        p <- p + scale_colour_gradientn(colours="blue") +
             geom_point(data = over_thresholds_dataSubset, size = point_size, alpha = 0.9, aes_string(colour = fieldName))
    } else {
        p <- p + scale_colour_gradientn(colours=rev(c("red", "orange", "yellow", "green", "cyan", "blue"))) +
             geom_point(data = over_thresholds_dataSubset, size = point_size, alpha = 0.9, aes_string(colour = fieldName))
    }

    plot_title <- paste0(printFieldName, " ", runName, " (", printFieldName, " = ", fullFieldName, ")\n",
                         "Colour thresholds are light grey: <=", thresholdOne, "; grey: >", thresholdOne, " and <", thresholdTwo)
    p <- p + labs(colour=printFieldName, size="", fill="", shape="") + ggtitle(plot_title)
    p <- p + coord_fixed(ratio = 1/1)
    #p <- p + geom_text(size = 1, aes(label = myData[[fieldName]])) # display value as text inside point

    # add cell labels
    if(length(labelPositions) > 0){
        #library("ggrepel")
        dataSubset <- myData[labelPositions,]
        if(!is.null(labelField)){
            subsetNames <- dataSubset[[labelField]] # include the cell names also: subsetNames <- c(rownames(dataSubset), dataSubset[[labelField]])
        } else {
            subsetNames <- rownames(dataSubset)
        }
        options(ggrepel.max.overlaps = Inf)
        p <- p + ggrepel::geom_text_repel(data = dataSubset, aes(label = subsetNames), size = 1.6,
                                          segment.color = 'black', segment.size = 0.1, 
                                          force = 2)
    }
    if(length(colourPositions) > 0){
        dataSubset <- myData[colourPositions,]
        p <- p + geom_point(data = dataSubset, colour = "green", alpha = 0.5)
    }

    # axis breaks
    minValue_1 <- min(round(myData[[xAxisName]]))
    maxValue_1 <- max(round(myData[[xAxisName]]))
    minValue_2 <- min(round(myData[[yAxisName]]))
    maxValue_2 <- max(round(myData[[yAxisName]]))
    minV1_mod <- minValue_1 %% 10
    minValue_1 <- minValue_1 - minV1_mod
    maxV1_mod <- maxValue_1 %% 10
    maxValue_1 <- maxValue_1 - maxV1_mod
    minV2_mod <- minValue_2 %% 10
    minValue_2 <- minValue_2 - minV2_mod
    maxV2_mod <- maxValue_2 %% 10
    maxValue_2 <- maxValue_2 - maxV2_mod
    breaks_x <- c(seq(minValue_1, maxValue_1, by = 10))
    breaks_y <- c(seq(minValue_2, maxValue_2, by = 10))

    p <- p + scale_y_continuous(breaks = breaks_y)
    p <- p + scale_x_continuous(breaks = breaks_x)

    print(p)
    ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 300, width = 15, height = 13.5, units = "in", limitsize = FALSE)

    return()
}

#' Print a list of t-SNEs each overlaid with expression of a specic genes; multiple t-SNE prints can generated
#' using a list of genes.
#'
#' This function is a wrapper for calling plotFieldWithThresholds to print a list of t-SNEs described 
#' in tSNE_list_with_annotations, each overlaid with expression for one gene taken from CDX_tRawCounts_list.
#'
#' @param gene_names_to_print List of one or more gene names to print expression levels for on one or more t-SNEs.
#' @param tSNE_with_annotations Each data frame column 1 and 2 should be Y1 and Y2 t-SNE info respectively and
#'     additional columns containing t-SNE scores/annotations e.g. stemness score or total transcription; row names
#'     should be single-cell IDs.
#' @param tExprMatr List of expression matrices describing cancer single-cell population; for each
#'     matrix, column names should be gene names, and row names should be cell IDs.
#' @param run_name Run name, e.g. 'tSNE_070'; will be used in the title of the plot and the output tif name.
#' @param outDir Output directory to use.
#' @param point_size Size of points in plot.
#'     value used.
#' @return No return value.
#' @export
print_tSNE_coloured_by_gene <- function(gene_names_to_print, tSNE_with_annotations, tExprMatr, run_name, outDir, point_size = NULL,
                                        HUGO_abbreviations_and_fullnames = NULL, HUGOFullNameFilesDir = NULL, use_biomart = FALSE){

    options(scipen=999) # To print non-scientific numbers

    if(is.null(point_size)){
        point_size <- 1.2
    }
    
    # Get gene functions
    #use_biomart <- FALSE
    if(use_biomart == TRUE){
        if(!exists(mart)){
            #library("biomaRt")
            mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
        }
        gene_names_to_print_and_functions <- biomaRt::getBM(attributes = c("hgnc_symbol", "description"), filters = 'hgnc_symbol', gene_names_to_print, mart = mart)
    } else {
        gene_names_to_print_and_functions <- getGenesAndFullNamesDF(gene_names_to_print, HUGO_abbreviations_and_fullnames, HUGOFullNameFilesDir)
    }
    
    # Print the t-SNE coloured by genes of interest
    dir.create(file.path(outDir), showWarnings = FALSE)

    for(i in 1:length(gene_names_to_print)){
        #i <- 1
        this_gene_name <- gene_names_to_print[[i]]
        this_full_gene_name <- gene_names_to_print_and_functions[i,2]
        cat(paste0("  ", this_gene_name, " [", this_full_gene_name, "]", "\n"))
        this_gene_col <- which(colnames(tExprMatr) %in% this_gene_name)
        if(length(this_gene_col) > 0){
            # Add the gene's expression values to an extra field in a temp DF created just for printing
            #this_gene_expr_name <- paste0(this_gene_name, "_expr")
            this_gene_expr <- as.vector(tExprMatr[,this_gene_col]) # expression values for gene to print
            numberOfColumns <- length(tSNE_with_annotations)
            tSNE_with_annotations <- data.frame(tSNE_with_annotations, this_gene_expr)            
            colnames(tSNE_with_annotations)[(numberOfColumns + 1)] <- this_gene_name
    
            # Print the plot
            plot_field_with_thresholds(tSNE_with_annotations,
                                       "Y1", "Y2",
                                       this_gene_name, this_full_gene_name,
                                       thresholdOne <- 0, # lower threshold for cells in grey
                                       thresholdTwo <- 1, # upper threshold
                                       run_name, outDir, 
                                       NULL, NULL, NULL, highestOnTop <- TRUE,
                                       point_size)
            tSNE_with_annotations[(numberOfColumns + 1)] <- NULL
        } else {
            cat(paste0(this_gene_name, " not found in this dataset so skipping.\n"))
        }
    }

    return()
}


#' Function to round up numbers to nearest numbers divisible by a specified value.
#'
#' This function rounds up to a nearest number, e.g. value = 7, by = 10, will round up to 10.
#'
#' @param value Value to round up
#' @param by Round up to the nearest value divisible by this number
#' @return Rounded up value
#' @export
round_up <- function(value, by){

    rounded_up_val <- by + (value - (value %% by))
    #rounded_up_val <- value + (by - value %% by)
    
    return(rounded_up_val)
}


#' Function to round down numbers to nearest numbers divisible by a specified value.
#'
#' This function rounds down to a nearest number, e.g. value = -7, by = 10, will round down to -10.
#'
#' @param value Value to round down
#' @param by Round down to the nearest value divisible by this number
#' @return Rounded down value
#' @export
round_down <- function(value, by){

    rounded_down_val <- value - (value %% by)

    return(rounded_down_val)
}
