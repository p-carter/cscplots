#' Run machine learning stemness predictors.
#'
#' This function runs ML classifiers on single-cell expression matrices to predict 
#' biological characteristics using annotationed expression data as the training data
#' for the classifier. In particular, use this function to predict stemness using one 
#' of three random forest methods ('ranger', 'rf', or 'randomForest'), training 
#' them with either 'HipSci' or 'PCBC' stemness expression data to predict stemness
#' probability for each cell in the cancer single-cell populations.
#'
#' @param model_type Type of classifier to use: 'ranger', 'rf', or 'randomForest'.
#' @param sample_tCounts_to_use_list List of single-cell cancer samples in expression matrix where column headers are gene names, and
#'     row names are cancer cell IDs.
#' @param stemData Stemness expression data for training where column one must be 'Cell_type' and the rest are gene names, and row names are sample names.
#' @param stemness_type Name of stemness data; choices are 'hipsci' or 'PCBC'.
#' @param number_of_genes_to_compare Number of genes to compare between the samples and the stemness data; recommended = 10000 to 16000.
#' @param sample_data_type Type of single-cell expression data being used e.g. 'raw' or 'tpm'; just a name though, can be anything.
#' @param stemness_out_dir Directory to save results from this function.
#' @param show_warnings Optional: Show warnings when building predictor model, FALSE is default.
#' @return A list of prediction results i.e. one for each scRNA-Seq sample.
#' @export
run_ML_stemness_predictors <- function(model_type, sample_tCounts_to_use_list, stemData, stemness_type, number_of_genes_to_compare, sample_data_type, stemness_out_dir, show_warnings = FALSE){

    #library("caret")
    #library("ranger")
    #library("e1071")
    #library("randomForest")
    ##library("gbm")

    results_list <- list()

    #cat(paste0("Model:", model_type, "\n"))
    for(i in 1:length(sample_tCounts_to_use_list)){
        #i <- 1
        this_sample_name <- names(sample_tCounts_to_use_list)[i]
        cat(paste0("[", i, "]:", this_sample_name, "\n"))

        sample_tCounts <- sample_tCounts_to_use_list[[i]]        
        if(model_type == "randomForest"){
            colnames(sample_tCounts) <- gsub("-", "_", colnames(sample_tCounts))
            colnames(stemData) <- gsub("-", "_", colnames(stemData))
        }
        sample_gene_IDs <- colnames(sample_tCounts)
        stemData_gene_IDs <- colnames(stemData)

        # Find overlapping genes between sample and stemness data
        sample_in_stemData_pos <- which(sample_gene_IDs %in% stemData_gene_IDs)
        stemData_in_sample_pos <- which(stemData_gene_IDs %in% sample_gene_IDs)        
        #print(identical(sort(sample_gene_IDs[sample_in_stemData_pos]), sort(stemData_gene_IDs[stemData_in_sample_pos])))
        cat(paste0("Overlapping genes between sample and stemness data = ", length(sample_in_stemData_pos), "\n"))

        # Subset on overlapping genes
        sample_tCounts_subset <- sample_tCounts[,sample_in_stemData_pos]
        stemData_subset <- stemData[,c(1, stemData_in_sample_pos)]        
        sample_tCounts_subset_order <- order(colnames(sample_tCounts_subset))
        stemData_subset_order <- order(colnames(stemData_subset[2:ncol(stemData_subset)]))
        sample_tCounts_subset_ordered <- sample_tCounts_subset[,sample_tCounts_subset_order]
        stemData_subset_ordered <- stemData_subset[,c(1, 1+stemData_subset_order)]
        #print(identical(colnames(sample_tCounts_subset_ordered), colnames(stemData_subset_ordered[2:ncol(stemData_subset_ordered)])))

        # Subset sample and stemness data to max number of genes
        if(length(sample_in_stemData_pos) > number_of_genes_to_compare){
            #sample_subset_pos <- 1:number_of_genes_to_compare # take the first n genes
            #stem_subset_pos <- 1:(number_of_genes_to_compare + 1)
            sample_subset_pos <- sort(sample(1:ncol(sample_tCounts_subset_ordered), number_of_genes_to_compare)) # use random positions
            #sample_subset_pos_file <- paste0(stemness_out_dir, "sample_subset_pos-", this_sample_name, "-", train_type, "_", model_type, "-", number_of_genes_to_compare, "_", sample_data_type, ".txt")
            #write.table(sample_subset_pos, file = sample_subset_pos_file, col.names = F, row.names = F, sep = "\t", append = FALSE, quote = FALSE)
            stem_subset_pos <- c(1,(sample_subset_pos + 1))
            sample_tCounts_subset_to_use <- sample_tCounts_subset_ordered[,sample_subset_pos]
            stemData_subset_to_use <- stemData_subset_ordered[,stem_subset_pos]
        } else {
            sample_tCounts_subset_to_use <- sample_tCounts_subset_ordered
            stemData_subset_to_use <- stemData_subset_ordered
        }
        #cat("Genes match between subsets to use:"); print(identical(colnames(sample_tCounts_subset_to_use), colnames(stemData_subset_to_use[2:ncol(stemData_subset_to_use)])))

        #inTrain <-  caret::createDataPartition(stemData_subset_to_use[,1], p=0.8, list=FALSE)[,1]
        inTrain  <-  caret::createDataPartition(stemData_subset_to_use[,1], p=0.6, list=FALSE)[,1]
        training <- stemData_subset_to_use[inTrain,]
        testing  <-  stemData_subset_to_use[-inTrain,]

        cat("Start of predictions\n")        
        #controlProbs <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 5, classProbs = TRUE)
        #controlProbs <- caret::trainControl(method = "repeatedcv", number = 2, repeats = 1, classProbs = TRUE)
        #controlProbs <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs = TRUE)
        controlProbs <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 2, classProbs = TRUE)

        #cat("Fitting\n")
        #if(model_type == "gbm"){
        #    thisFit <- caret::train(Cell_type ~ ., data = training, method = model_type, trControl = controlProbs, train.fraction = 0.8)        
        #}
        if(model_type == "randomForest"){
            if(show_warnings == TRUE){
                thisFit <- randomForest::randomForest(Cell_type ~ ., data = training)
            } else {
                thisFit <- suppressWarnings(randomForest::randomForest(Cell_type ~ ., data = training))
            }
        } else {
            if(show_warnings == TRUE){
                thisFit <- caret::train(Cell_type ~ ., data = na.omit(training), method = model_type, trControl = controlProbs)
            } else {
                thisFit <- suppressWarnings(caret::train(Cell_type ~ ., data = na.omit(training), method = model_type, trControl = controlProbs))
            }
        }
        #save(thisFit, file = paste0(stemness_out_dir, this_sample_name, "-", stemness_type, "_", model_type, "_fit-", number_of_genes_to_compare, "_", sample_data_type, ".RData"))

        cat("Predictions: testing \n")       
        thisTestingPredict <- stats::predict(thisFit, newdata = testing)
        #print(thisTestingPredict)
        caret::confusionMatrix(thisTestingPredict, testing$Cell_type)
        
        cat("Probabilities: testing\n")
        thisTestingPredict_prob <- stats::predict(thisFit, newdata = testing, type = "prob")
        #print(thisTestingPredict_prob)

        cat(paste0("Predictions: ", this_sample_name, "\n"))
        thisSamplePredict <- stats::predict(thisFit, newdata = sample_tCounts_subset_to_use)
        #print(thisSamplePredict)
        
        cat(paste0("Probabilities: ", this_sample_name, "\n"))
        thisSamplePredict_prob <- stats::predict(thisFit, newdata = sample_tCounts_subset_to_use, type = "prob")
        rownames(thisSamplePredict_prob) <- rownames(sample_tCounts_subset_to_use)
        #print(thisSamplePredict_prob)

        print_df <- data.frame("Cell_names"=rownames(thisSamplePredict_prob),thisSamplePredict_prob)
        if(stemness_type == 'hipsci'){
            write.table(print_df, file = paste0(stemness_out_dir, "/", this_sample_name, "-", stemness_type, "_", model_type, "_predictions-", number_of_genes_to_compare, "_", sample_data_type, ".txt"),
                        col.names = TRUE, row.names = FALSE, sep = "\t", append = FALSE, quote = FALSE)
        } else if(stemness_type == 'PCBC'){
            write.table(print_df, file = paste0(stemness_out_dir, "/", this_sample_name, "-", stemness_type, "_", model_type, "_predictions-", "H9_IPS_SC", "_", sample_data_type, ".txt"),
                        col.names = TRUE, row.names = FALSE, sep = "\t", append = FALSE, quote = FALSE)
        } else {  cat("Stemness name not hipsci or PCBC\n")  }

        if(model_type == "randomForest"){
            thisSamplePredict_prob <- as.data.frame(thisSamplePredict_prob)
        }

        results_list[[i]] <- thisSamplePredict_prob
        names(results_list)[i] <- this_sample_name
    }

    cat("ML stemness predictions completed\n")

    return(results_list)
}


#' Run machine learning stemness predictors.
#'
#' This function runs ML classifiers on single-cell expression matrices to predict 
#' biological characteristics using annotationed expression data as the training data
#' for the classifier. In particular, use this function to predict stemness using one 
#' of three random forest methods ('ranger', 'rf', or 'randomForest') or gbm, training 
#' them with expressionData with cell type annotations to predict cell type probability
#' for each cell in the cancer single-cell populations.
#'
#' @param model_type Type of classifier to use: 'ranger', 'rf', 'randomForest' or 'gbm'.
#' @param sample_tCounts_to_use_list List of single-cell cancer samples in expression matrix where column headers are gene names, and
#'     row names are cancer cell IDs.
#' @param cellTypeData Expression data for training where column one must be 'Cell_type' and the rest are gene names, and row names are sample names.
#' @param expressionData_type Name of expression data e.g. 'NSCLC', 'melanoma', 'hipsci', etc.; can be anything.
#' @param number_of_genes_to_compare Number of genes to compare between the samples and the expression data; recommended = 16000 to 10000.
#' @param sample_data_type Type of single-cell expression data being used e.g. 'raw' or 'tpm'; just a name though, can be anything.
#' @param cellType_out_dir Directory to save results from this function.
#' @param show_warnings Optional: Show warnings when building predictor model, FALSE is default.
#' @return A list of prediction results i.e. one for each scRNA-Seq sample.
#' @export
run_ML_cellType_predictors <- function(model_type, sample_tCounts_to_use_list, cellTypeData, expressionData_type, number_of_genes_to_compare, sample_data_type, cellType_out_dir, show_warnings = FALSE){

    #library("caret")
    #library("ranger")
    #library("e1071")
    #library("randomForest")
    ##library("gbm")

    results_list <- list()

    #cat(paste0("Model:", model_type, "\n"))
    for(i in 1:length(sample_tCounts_to_use_list)){
        #i <- 1
        this_sample_name <- names(sample_tCounts_to_use_list)[i]
        cat(paste0("[", i, "]:", this_sample_name, "\n"))

        sample_tCounts <- sample_tCounts_to_use_list[[i]]        
        if(model_type == "randomForest"){
            colnames(sample_tCounts) <- gsub("-", "_", colnames(sample_tCounts))
            colnames(cellTypeData) <- gsub("-", "_", colnames(cellTypeData))
        }
        sample_gene_IDs <- colnames(sample_tCounts)
        cellTypeData_gene_IDs <- colnames(cellTypeData)

        # Find overlapping genes between sample and stemness data
        sample_in_cellTypeData_pos <- which(sample_gene_IDs %in% cellTypeData_gene_IDs)
        cellTypeData_in_sample_pos <- which(cellTypeData_gene_IDs %in% sample_gene_IDs)        
        #print(identical(sort(sample_gene_IDs[sample_in_cellTypeData_pos]), sort(cellTypeData_gene_IDs[cellTypeData_in_sample_pos])))
        cat(paste0("Overlapping genes between sample and stemness data = ", length(sample_in_cellTypeData_pos), "\n"))

        # Subset on overlapping genes
        sample_tCounts_subset <- sample_tCounts[,sample_in_cellTypeData_pos]
        cellTypeData_subset <- cellTypeData[,c(1, cellTypeData_in_sample_pos)]        
        sample_tCounts_subset_order <- order(colnames(sample_tCounts_subset))
        cellTypeData_subset_order <- order(colnames(cellTypeData_subset[2:ncol(cellTypeData_subset)]))
        sample_tCounts_subset_ordered <- sample_tCounts_subset[,sample_tCounts_subset_order]
        cellTypeData_subset_ordered <- cellTypeData_subset[,c(1, 1+cellTypeData_subset_order)]
        #print(identical(colnames(sample_tCounts_subset_ordered), colnames(cellTypeData_subset_ordered[2:ncol(cellTypeData_subset_ordered)])))

        # Subset sample and stemness data to max number of genes
        if(length(sample_in_cellTypeData_pos) > number_of_genes_to_compare){
            #sample_subset_pos <- 1:number_of_genes_to_compare # take the first n genes
            #stem_subset_pos <- 1:(number_of_genes_to_compare + 1)
            sample_subset_pos <- sort(sample(1:ncol(sample_tCounts_subset_ordered), number_of_genes_to_compare)) # use random positions
            #sample_subset_pos_file <- paste0(cellType_out_dir, "sample_subset_pos-", this_sample_name, "-", train_type, "_", model_type, "-", number_of_genes_to_compare, "_", sample_data_type, ".txt")
            #write.table(sample_subset_pos, file = sample_subset_pos_file, col.names = F, row.names = F, sep = "\t", append = FALSE, quote = FALSE)
            stem_subset_pos <- c(1,(sample_subset_pos + 1))
            sample_tCounts_subset_to_use <- sample_tCounts_subset_ordered[,sample_subset_pos]
            cellTypeData_subset_to_use <- cellTypeData_subset_ordered[,stem_subset_pos]
        } else {
            sample_tCounts_subset_to_use <- sample_tCounts_subset_ordered
            cellTypeData_subset_to_use <- cellTypeData_subset_ordered
        }
        #cat("Genes match between subsets to use:"); print(identical(colnames(sample_tCounts_subset_to_use), colnames(cellTypeData_subset_to_use[2:ncol(cellTypeData_subset_to_use)])))

        #inTrain <-  caret::createDataPartition(cellTypeData_subset_to_use[,1], p=0.8, list=FALSE)[,1]
        inTrain  <-  caret::createDataPartition(cellTypeData_subset_to_use[,1], p=0.6, list=FALSE)[,1]
        training <- cellTypeData_subset_to_use[inTrain,]
        testing  <-  cellTypeData_subset_to_use[-inTrain,]

        cat("Start of predictions\n")        
        #controlProbs <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 5, classProbs = TRUE)
        #controlProbs <- caret::trainControl(method = "repeatedcv", number = 2, repeats = 1, classProbs = TRUE)
        #controlProbs <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs = TRUE)
        controlProbs <- caret::trainControl(method = "repeatedcv", number = 5, repeats = 2, classProbs = TRUE)

        #cat("Fitting\n")
        if(model_type == "gbm"){
            if(show_warnings == TRUE){
                thisFit <- caret::train(Cell_type ~ ., data = training, method = model_type, trControl = controlProbs, train.fraction = 0.8)
            } else {
                thisFit <- suppressWarnings(caret::train(Cell_type ~ ., data = training, method = model_type, trControl = controlProbs, train.fraction = 0.8))
            }
        }
        else if(model_type == "randomForest"){
            if(show_warnings == TRUE){
                thisFit <- randomForest::randomForest(Cell_type ~ ., data = training)
            } else {
                thisFit <- suppressWarnings(randomForest::randomForest(Cell_type ~ ., data = training))
            }
        } else {
            if(show_warnings == TRUE){
                thisFit <- caret::train(Cell_type ~ ., data = na.omit(training), method = model_type, trControl = controlProbs)
            } else {
                thisFit <- suppressWarnings(caret::train(Cell_type ~ ., data = na.omit(training), method = model_type, trControl = controlProbs))
            }
        }
        #save(thisFit, file = paste0(cellType_out_dir, this_sample_name, "-", expressionData_type, "_", model_type, "_fit-", number_of_genes_to_compare, "_", sample_data_type, ".RData"))

        cat("Predictions: testing \n")       
        thisTestingPredict <- stats::predict(thisFit, newdata = testing)
        #print(thisTestingPredict)
        caret::confusionMatrix(thisTestingPredict, testing$Cell_type)
        
        cat("Probabilities: testing\n")
        thisTestingPredict_prob <- stats::predict(thisFit, newdata = testing, type = "prob")
        #print(thisTestingPredict_prob)

        cat(paste0("Predictions: ", this_sample_name, "\n"))
        thisSamplePredict <- stats::predict(thisFit, newdata = sample_tCounts_subset_to_use)
        #print(thisSamplePredict)
        
        cat(paste0("Probabilities: ", this_sample_name, "\n"))
        thisSamplePredict_prob <- stats::predict(thisFit, newdata = sample_tCounts_subset_to_use, type = "prob")
        rownames(thisSamplePredict_prob) <- rownames(sample_tCounts_subset_to_use)
        #print(thisSamplePredict_prob)

        print_df <- data.frame("Cell_names"=rownames(thisSamplePredict_prob),thisSamplePredict_prob)
        write.table(print_df, file = paste0(cellType_out_dir, "/", this_sample_name, "-", expressionData_type, "_", model_type, "_predictions-", number_of_genes_to_compare, "_", sample_data_type, ".txt"),
            col.names = TRUE, row.names = FALSE, sep = "\t", append = FALSE, quote = FALSE)

        if(model_type == "randomForest"){
            thisSamplePredict_prob <- as.data.frame(thisSamplePredict_prob)
        }

        results_list[[i]] <- thisSamplePredict_prob
        names(results_list)[i] <- this_sample_name
    }

    cat("ML cell type predictions completed\n")

    return(results_list)
}

