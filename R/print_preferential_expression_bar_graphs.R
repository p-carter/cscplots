#' Print preferential expression bar graphs for subset of single-cells.
#'
#' This function prints multiple bar graphs to display preferential expression
#' in a defined subset of cells e.g. a group of cancer stem cells.
#'
#' @param pCSC_group_name Name of CSC group to print preferential expression plots for
#'     e.g. 'pCSC_group_1_SC4_Talazoparib.LB17011', but can be anything.
#' @param pCSC_group_positions Row number of CSCs in the CSC group in the expression matrix
#'     tExprMatr.
#' @param HUGO_groups_sets_of_columnNumbers List of lists, i.e. a list of HUGO gene groups,
#'     with each named gene group containing positions (row numbers) in the expression matrix
#'     tExprMatr of all the gene group members present in the expression matrix.
#' @param tExprMatr Expression matrix describing a cancer single-cell population; column
#'     names should be gene names, row names should be cancer cell names.
#' @param tumour_only_positions Positions (row numbers) in the expression matrix tExprMatr
#'     thought to be tumour cells e.g. when the population contains cells annotated or predicted
#'     as other cell types a subset may be defined, or if all cells are thought to be malignant
#'     cancer cells this may be all cells.
#' @param this_patient_tumour_only_positions Optional: Expression matrix tExprMatr positions
#'     (row numbers) for tumour cells from a specific patient i.e. that the CSCs for the
#'     preferential expression plots are taken from. So, the plot needs to fulfil the criteria
#'     that the CSCs are from only one patient, but the tumour cell population is from multiple
#'     patients.
#' @param outDir Output directory for plots.
#' @param tissue_type Name of tissue type e.g. 'CDX'; can be anything. 
#' @param family_type_name Name of the gene group type; usually 'HUGO'.
#' @param patient_source_name Optional: Name of the patient sample that CSC cells are from;
#'     just name for in the plot title, can be anything.
#' @param frag_count_method Expression measurement type, e.g. 'raw'; but just a label, can be anything.
#' @param no_ribo Optional: Set to 1 to remove ribosomal transcripts from some of the bar charts
#' @param separate_row_for_high_expression_transcripts Optional: Set to 1 to print very high
#'     expression transcripts or gene groups on a separate row in some of the bar charts; helpful
#'     when there is a big difference between the most highly expressed and those on the same
#'     bar chart row.
#' @param max_number_of_graph_rows_stacked_transcripts Optional: Set to manually override the
#'     maximum number of rows shown in the stacked transcript plot. 
#' @param max_number_of_graph_rows_stacked_families Optional: Set to manually override the
#'     maximum number of rows shown in the stacked family plot. 
#' @return List of highest ranked preferentially expressed genes and gene families:
#'     1 = ranked transcripts from single transcript barcharts;
#'     2 = ranked families from family barcharts;
#'     3 = ranked transcripts from family barcharts.
#' @export
print_preferential_expression_bar_graphs <- function(pCSC_group_name,
                                                     pCSC_group_positions,
                                                     HUGO_groups_sets_of_columnNumbers,
                                                     tExprMatr,
                                                     tumour_only_positions,
                                                     this_patient_tumour_only_positions = NULL,
                                                     outDir,
                                                     tissue_type,
                                                     family_type_name,
                                                     patient_source_name,
                                                     frag_count_method,
                                                     no_ribo = NULL,
                                                     separate_row_for_high_expression_transcripts = NULL,
                                                     max_number_of_graph_rows_stacked_transcripts = NULL,
                                                     max_number_of_graph_rows_stacked_families = NULL,
                                                     HUGO_abbreviations_and_fullnames = NULL,
                                                     HUGOFullNameFilesDir = NULL){

    #library("ggplot2")
    #library("ggrepel")
    #library("RColorBrewer")
    #require(scales)
    #library("cowplot")

    if(is.null(no_ribo)){
        no_ribo <- FALSE
    }
    if(is.null(separate_row_for_high_expression_transcripts)){
        separate_row_for_high_expression_transcripts <- TRUE
    }
    if(!is.null(this_patient_tumour_only_positions)){
        if(length(this_patient_tumour_only_positions) == 1){
            if(this_patient_tumour_only_positions == ''){
                this_patient_tumour_only_positions <- NULL
            }
        }
    }
    
    #cat(paste0("Output directory = '", outDir, "preferential_expression/'\n"))
    cat(paste0("Output directory = '", outDir, "'\n"))
    
    ## RANKINGS ##########

    #1) STACKED BARCHARTS FOR SINGLE TRANSCRIPTS - rank top 100 transcripts
    #2) STACKED FAMILY EXPRESSION BARCHART SINGLE BARS - rank top 100 families
    #3)                                                - also rank by each transcript's share of top family rank (top family occurrence only)
    
    ranks1_STACKED_BARCHARTS_FOR_SINGLE_TRANSCRIPTS                                                                        <- NULL
    ranks2_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__FAMILIES                                                        <- NULL
    ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted_no_duplicates_no_zeros <- NULL

    ############################################################
    
    #1# TRIPLET SINGLE TRANSCRIPT BARCHARTS

    printTripletSingleTranscriptBarcharts(tExprMatr,
                                          pCSC_group_positions,
                                          tumour_only_positions,
                                          pCSC_group_name,
                                          outDir,
                                          tissue_type,
                                          family_type_name,
                                          patient_source_name,
                                          no_ribo <- TRUE,
                                          number_of_transcripts_to_show <- 100,
                                          HUGO_abbreviations_and_fullnames,
                                          HUGOFullNameFilesDir)

    ############################################################

    #2# BAR CHART SINGLE TRANSCRIPT NO STACKING

    printBarchartSingleTranscriptNoStacking(tExprMatr,
                                            pCSC_group_positions,
                                            pCSC_group_name,
                                            outDir,
                                            tissue_type,
                                            family_type_name,
                                            patient_source_name,
                                            no_ribo <- TRUE,
                                            number_of_transcripts_to_show <- 100,
                                            HUGO_abbreviations_and_fullnames,
                                            HUGOFullNameFilesDir)
    
    ############################################################

    #3# STACKED BARCHARTS FOR SINGLE TRANSCRIPTS

    ranks1_STACKED_BARCHARTS_FOR_SINGLE_TRANSCRIPTS <- stackedBarchartsForSingleTranscripts(tExprMatr,                           
                                                                                            pCSC_group_positions,                  
                                                                                            tumour_only_positions,                 
                                                                                            this_patient_tumour_only_positions,
                                                                                            pCSC_group_name,
                                                                                            tissue_type,
                                                                                            family_type_name,
                                                                                            patient_source_name,
                                                                                            separate_row_for_high_expression_transcripts,
                                                                                            no_ribo,
                                                                                            max_number_of_graph_rows_stacked_transcripts,
                                                                                            HUGO_abbreviations_and_fullnames,
                                                                                            HUGOFullNameFilesDir)

    ############################################################
    
    #4# STACKED FAMILY EXPRESSION BARCHART SINGLE BARS

    family_ranks <- stackedFamilyExpressionBarchartSingleBars(tExprMatr,
                                                              pCSC_group_positions,
                                                              tumour_only_positions,
                                                              this_patient_tumour_only_positions,
                                                              HUGO_groups_sets_of_columnNumbers,
                                                              pCSC_group_name,                             
                                                              outDir,                                      
                                                              tissue_type,                                 
                                                              family_type_name,                            
                                                              patient_source_name,
                                                              frag_count_method,                           
                                                              separate_row_for_high_expression_transcripts,
                                                              max_number_of_graph_rows_stacked_families)
    
    ranks2_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__FAMILIES <- family_ranks[[1]]
    ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted_no_duplicates_no_zeros <- family_ranks[[2]]

    ############################################################

    ## RANKINGS ##########

    these_ranks <- list()

    these_ranks[[1]]  <- ranks1_STACKED_BARCHARTS_FOR_SINGLE_TRANSCRIPTS
    these_ranks[[2]]  <- ranks2_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__FAMILIES
    these_ranks[[3]]  <- ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted_no_duplicates_no_zeros
    
    names(these_ranks) <- c(
        "ranks1_transcripts_from_single_transcript_barcharts",
        "ranks2_families_from_family_barcharts",
        "ranks3_transcripts_from_family_barcharts"
    )

    #cat("# # # # # # # # # EXPRESSION BAR CHARTS COMPLETED # # # # # # # # #\n")
    cat("Expression bar charts completed\n\n")

    return(these_ranks)
}


#' Print expression bar graph for genes/transcripts most highly expressed in selected
#'     putative CSCs vs tumour cells vs all cells.
#'
#' This function prints a bar graph to display relative expression in the genes most
#'     highly expressed in a selected set of putative CSCs, along with their expression 
#'     in tumour cells and all cells. The order is by mean transcript expression in the
#'     pCSC cells i.e. highest first. If the no_ribo parameter is used, the same bar
#'     graph is printed with all transcripts with "ribosomal" in their full name removed
#'     from the results.
#'
#' @param tExprMatr Expression matrix describing a cancer single-cell population; column
#'     names should be gene names, row names should be cancer cell names.
#' @param pCSC_group_positions Row number of CSCs in the CSC group in the expression matrix
#'     tExprMatr.
#' @param tumour_only_positions Positions (row numbers) in the expression matrix tExprMatr
#'     thought to be tumour cells e.g. when the population contains cells annotated or predicted
#'     as other cell types a subset may be defined, or if all cells are thought to be malignant
#'     cancer cells this may be all cells.
#' @param pCSC_group_name Name of CSC group to print preferential expression plots for
#'     e.g. 'pCSC_group_1_SC4_Talazoparib.LB17011', but can be anything.
#' @param outDir Output directory for plots.
#' @param tissue_type Name of tissue type e.g. 'CDX'; can be anything. 
#' @param family_type_name Name of the gene group type; usually 'HUGO'.
#' @param patient_source_name Optional: Name of the patient sample that CSC cells are from;
#'     just name for in the plot title, can be anything.
#' @param no_ribo Optional: Set to 1 to remove ribosomal transcripts from some of the bar charts
#' @param number_of_transcripts_to_show Optional: Number of transcripts to include.
#' @return Nothing returned.
#' @export
printTripletSingleTranscriptBarcharts <- function(tExprMatr,
                                                  pCSC_group_positions,
                                                  tumour_only_positions,
                                                  pCSC_group_name,
                                                  outDir,
                                                  tissue_type,
                                                  family_type_name,
                                                  patient_source_name,
                                                  no_ribo = TRUE,
                                                  number_of_transcripts_to_show = NULL,
                                                  HUGO_abbreviations_and_fullnames = NULL,
                                                  HUGOFullNameFilesDir = NULL
                                                  ){

    #1# TRIPLET SINGLE TRANSCRIPT BARCHARTS
    cat("Printing: Single transcripts barcharts (triplets)\n")

    #require(scales)

    ####
    
    pCSCs_tExprMatr                    <- tExprMatr[pCSC_group_positions,]
    tumour_only_tExprMatr              <- tExprMatr[tumour_only_positions,]

    pCSCs_tExprMatr_colmeans                  <- colMeans(pCSCs_tExprMatr)
    pCSCs_tExprMatr_colmeans_sorted_rev       <- rev(sort(pCSCs_tExprMatr_colmeans))
    pCSCs_tExprMatr_colmeans_sorted_rev_names <- names(pCSCs_tExprMatr_colmeans_sorted_rev)
    pCSCs_tExprMatr_colmeans_rev_order        <- rev(order(pCSCs_tExprMatr_colmeans))

    # Order the full set columns by the colmeans of the pCSC subset & tumour only subset
    tExprMatr_rev_ordered_by_pCSC_group_colmeans                      <- tExprMatr[,pCSCs_tExprMatr_colmeans_rev_order]
    tumour_only_tExprMatr_rev_ordered_by_pCSC_group_colmeans          <- tumour_only_tExprMatr[,pCSCs_tExprMatr_colmeans_rev_order]
    tExprMatr_colmeans_rev_ordered_by_pCSC_group_colmeans             <- colMeans(tExprMatr_rev_ordered_by_pCSC_group_colmeans)
    tumour_only_tExprMatr_colmeans_rev_ordered_by_pCSC_group_colmeans <- colMeans(tumour_only_tExprMatr_rev_ordered_by_pCSC_group_colmeans)

    pCSCs_tExprMatr_colmeans_sorted_rev_full_names <- getGenesFullNames(pCSCs_tExprMatr_colmeans_sorted_rev_names, HUGO_abbreviations_and_fullnames, HUGOFullNameFilesDir)
    
    ####
    
    # df for printing
    prDF1 <- data.frame(
        transcript_abbreviation        = names(pCSCs_tExprMatr_colmeans_sorted_rev),
        transcript_name                = pCSCs_tExprMatr_colmeans_sorted_rev_full_names,
        pCSC_group_expression_average  = pCSCs_tExprMatr_colmeans_sorted_rev,
        tumour_only_expression_average = tumour_only_tExprMatr_colmeans_rev_ordered_by_pCSC_group_colmeans,
        all_samples_expression_average = tExprMatr_colmeans_rev_ordered_by_pCSC_group_colmeans,            
        stringsAsFactors = FALSE
    )
    rownames(prDF1) <- NULL

    prDF2 <- data.frame(
        cell_group_type    = c(rep(colnames(prDF1)[3], length(prDF1[,3])), rep(colnames(prDF1)[4], length(prDF1[,4])), rep(colnames(prDF1)[5], length(prDF1[,5]))),
        transcript_name    = c(prDF1[,2], prDF1[,2], prDF1[,2]),
        expression_average = c(prDF1[,3], prDF1[,4], prDF1[,5]),
        stringsAsFactors = FALSE
    )

    # families retain order (reorder to group together families for ggplot bar format)
    prDF2o <- prDF2[order(match(prDF2[,2], prDF2[,2])),]
    prDF2o$transcript_name <- factor(prDF2o$transcript_name, levels=unique(prDF2o$transcript_name))
    rownames(prDF2o) <- NULL

    if(is.null(number_of_transcripts_to_show)){
        number_of_transcripts_to_show <- 100
    }
    
    triplet_single_to_show_len <- number_of_transcripts_to_show * 3
    df_to_use <- prDF2o[1:triplet_single_to_show_len,]

    tif_name <- paste0(pCSC_group_name, "_transcript_expression_bargraph-triplets.tif")
    #outFile <- paste0(outDir, "preferential_expression/", tif_name)
    outFile <- paste0(outDir, tif_name)
    x_label <- "Transcript"
    y_label <- "Average Expression"
    plot_title <- paste0(toupper(tissue_type), ": [", family_type_name, "] Transcripts for ", pCSC_group_name, " (", length(pCSC_group_positions), " cells)", " from samples: ", patient_source_name)
    axis_text_size <- 6
    legend_title <- "Cells Group"
    legend_lab1  <- pCSC_group_name
    legend_lab2  <- "Tumour"
    legend_lab3  <- "All cell types"

    p <- ggplot2::ggplot(data=df_to_use, aes(x = transcript_name, y = expression_average, fill = factor(cell_group_type, levels=unique(cell_group_type)))) +
    geom_bar(stat="identity", color="black", position=position_dodge(), size=0.3, width=1) +
    theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),        
    panel.grid.minor = element_blank(),        
    panel.background = element_blank(),        
    panel.spacing = unit(0.1, "cm"),
    axis.text.x = element_text(color = "black", angle = 270, hjust = 0, vjust = 0.5, size = axis_text_size),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0, size = 10)) +
    scale_y_continuous(expand = c(0,0), labels = scales::comma) +
    scale_fill_brewer(palette = "Greens", direction=-1, name = legend_title, labels = c(legend_lab1, legend_lab2, legend_lab3))
    p <- p + labs(title = plot_title, x = x_label, y = y_label)

    ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 150, width = 24, height = 9, units = "in")
    print(p)

    ## no ribo - remove ribosomal transcripts ##################

    if(no_ribo == TRUE){
    
        # same plot as above but with no ribosomal
        prDF2o_copy <- prDF2o
        rows_to_remove_vec <- vector()
        for(i in 1:nrow(prDF2o)){
            #this_row_number <- (i * 3) - 2
            #cat(paste0(i, "\n"))
            this_trans_name <- as.vector(prDF2o[i,2])
            #cat(paste0("this trans name = ", this_trans_name, "\n"))
            if(substr(this_trans_name, 1, 9) == "ribosomal"){
                rows_to_remove_vec <- c(rows_to_remove_vec, i)
            }
        }
        prDF2o_no_ribo <- prDF2o[-c(rows_to_remove_vec),]
        df_to_use <- prDF2o_no_ribo[1:triplet_single_to_show_len,]
        
        # retain order
        df_to_use$transcript_name <- factor(df_to_use$transcript_name, levels=unique(df_to_use$transcript_name))
        rownames(df_to_use) <- NULL
        
        tif_name <- paste0(pCSC_group_name, "_transcript_expression_bargraph-triplets-no_ribo.tif")
        #outFile <- paste0(outDir, "preferential_expression/", tif_name)
        outFile <- paste0(outDir, tif_name)
        
        p <- ggplot2::ggplot(data = df_to_use, aes(x = transcript_name, y = expression_average, fill = factor(cell_group_type, levels=unique(cell_group_type)))) +
        geom_bar(stat = "identity", color = "black", position = position_dodge(), size = 0.3, width = 1) +
        theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),        
        panel.background = element_blank(),        
        panel.spacing = unit(0.1, "cm"),
        axis.text.x = element_text(color = "black", angle = 270, hjust = 0, vjust = 0.5, size = axis_text_size),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0, size = 10)) +
        scale_y_continuous(expand = c(0,0), labels = scales::comma) +
        scale_fill_brewer(palette = "Greens", direction=-1, name = legend_title, labels = c(legend_lab1, legend_lab2, legend_lab3))
        p <- p + labs(title = plot_title, x = x_label, y = y_label)
        
        ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 150, width = 24, height = 9, units = "in")
        print(p)
    }
        
    ## end of no ribo ##########################################

    return()
} 


#' Print expression bar graph for genes/transcripts most highly expressed in selected
#'     putative CSCs.
#'
#' This function prints a bar graph to display expression in the genes most highly expressed
#'     in a selected set of putative CSCs. The order is by mean transcript expression in the
#'     pCSC cells i.e. highest first. If the no_ribo parameter is used, the same bar
#'     graph is printed with all transcripts with "ribosomal" in their full name removed
#'     from the results.
#'
#' @param tExprMatr Expression matrix describing a cancer single-cell population; column
#'     names should be gene names, row names should be cancer cell names.
#' @param pCSC_group_positions Row number of CSCs in the CSC group in the expression matrix
#'     tExprMatr.
#' @param pCSC_group_name Name of CSC group to print preferential expression plots for
#'     e.g. 'pCSC_group_1_SC4_Talazoparib.LB17011', but can be anything.
#' @param outDir Output directory for plots.
#' @param tissue_type Name of tissue type e.g. 'CDX'; can be anything. 
#' @param family_type_name Name of the gene group type; usually 'HUGO'.
#' @param patient_source_name Optional: Name of the patient sample that CSC cells are from;
#'     just name for in the plot title, can be anything.
#' @param no_ribo Optional: Set to 1 to remove ribosomal transcripts from some of the bar charts
#' @param number_of_transcripts_to_show Optional: Number of transcripts to include.
#' @return Nothing returned.
#' @export
printBarchartSingleTranscriptNoStacking <- function(tExprMatr,    
                                                    pCSC_group_positions,
                                                    pCSC_group_name,
                                                    outDir,
                                                    tissue_type,
                                                    family_type_name,
                                                    patient_source_name,
                                                    no_ribo = TRUE,
                                                    number_of_transcripts_to_show = NULL,
                                                    HUGO_abbreviations_and_fullnames = NULL,
                                                    HUGOFullNameFilesDir = NULL
                                                    ){

    #2# BAR CHART SINGLE TRANSCRIPT NO STACKING
    cat("Printing: Single transcripts barcharts (no stacking)\n")

    ####

    pCSCs_tExprMatr                                <- tExprMatr[pCSC_group_positions,]                                                     
    pCSCs_tExprMatr_colmeans                       <- colMeans(pCSCs_tExprMatr)
    pCSCs_tExprMatr_colmeans_sorted_rev            <- rev(sort(pCSCs_tExprMatr_colmeans))
    pCSCs_tExprMatr_colmeans_sorted_rev_names      <- names(pCSCs_tExprMatr_colmeans_sorted_rev)
    pCSCs_tExprMatr_colmeans_sorted_rev_full_names <- getGenesFullNames(pCSCs_tExprMatr_colmeans_sorted_rev_names, HUGO_abbreviations_and_fullnames, HUGOFullNameFilesDir)
    
    ####
    
    prDF3 <- data.frame(
        transcript_abbreviation = names(pCSCs_tExprMatr_colmeans_sorted_rev),
        transcript_name         = pCSCs_tExprMatr_colmeans_sorted_rev_full_names,
        expression_average      = pCSCs_tExprMatr_colmeans_sorted_rev,
        stringsAsFactors = FALSE
    )
    rownames(prDF3) <- NULL
    if(is.null(number_of_transcripts_to_show)){
        number_of_transcripts_to_show <- 100
    }
    df_to_use <- prDF3[1:number_of_transcripts_to_show,]
    # retain order
    df_to_use$transcript_name <- factor(df_to_use$transcript_name, levels = unique(df_to_use$transcript_name))
    rownames(df_to_use) <- NULL

    tif_name <- paste0(pCSC_group_name, "_transcript_expression_bargraph.tif")
    #outFile <- paste0(outDir, "preferential_expression/", tif_name)
    outFile <- paste0(outDir, tif_name)
    x_label <- "Transcript"
    y_label <- "Average Expression"
    plot_title <- paste0(toupper(tissue_type), ": [", family_type_name, "] Transcript expression for ", pCSC_group_name, " (", length(pCSC_group_positions), " cells)", " from samples: ", patient_source_name)
    axis_text_size <- 6

    p <- ggplot2::ggplot(data=df_to_use, aes(x = transcript_name, y = expression_average, fill = expression_average)) +
    geom_bar(stat = "identity", color = "black", position = position_dodge(), size = 0.2, width = 1) +
    theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),        
    panel.grid.minor = element_blank(),        
    panel.background = element_blank(),        
    panel.spacing = unit(0.1, "cm"),
    axis.text.x = element_text(color = "black", angle = 270, hjust = 0, vjust = 0.5, size = axis_text_size),
    legend.position = "none",                                                                                        
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0, size = 10),
    plot.margin = ggplot2::margin(0.1, 0.4, 0.1, 0.1, "cm")) +                                                                  
    scale_y_continuous(expand = c(0,0), labels = scales::comma)  +
    #scale_fill_gradient(low = "pink", high = "red4")
    scale_fill_gradientn(colours = c("#0092FFFF", "#00FF92FF", "#49FF00FF", "#FFDB00FF", "#FF0000FF"))
    p <- p + labs(title = plot_title, x = x_label, y = y_label)                                                      

    ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 150, width = 12, height = 9, units = "in")
    print(p)

    ## no ribo #################################################

    if(no_ribo == TRUE){
        
        prDF3_rows_to_remove_vec <- vector()
        for(i in 1:nrow(prDF3)){
            this_trans_name <- as.vector(prDF3[i,2])
            if(substr(this_trans_name, 1, 9) == "ribosomal"){
                prDF3_rows_to_remove_vec <- c(prDF3_rows_to_remove_vec, i)
            }
        }
        prDF3_no_ribo <- prDF3[-c(prDF3_rows_to_remove_vec),]
        
        df_to_use <- prDF3_no_ribo[1:number_of_transcripts_to_show,]
        # retain order
        df_to_use$transcript_name <- factor(df_to_use$transcript_name, levels=unique(df_to_use$transcript_name))
        rownames(df_to_use) <- NULL
        
        tif_name <- paste0(pCSC_group_name, "_transcript_expression_bargraph-no_ribo.tif")
        #outFile <- paste0(outDir, "preferential_expression/", tif_name)
        outFile <- paste0(outDir, tif_name)
        
        p <- ggplot2::ggplot(data=df_to_use, aes(x = transcript_name, y = expression_average, fill = expression_average)) +
        geom_bar(stat="identity", color="black", position=position_dodge(), size=0.2, width=1) +
        theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank(),        
        panel.background = element_blank(),        
        panel.spacing = unit(0.1, "cm"),
        axis.text.x = element_text(color = "black", angle = 270, hjust = 0, vjust = 0.5, size = axis_text_size),
        legend.position = "none",                                                                                        
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0, size = 10),
        plot.margin = ggplot2::margin(0.1, 0.4, 0.1, 0.1, "cm")) +                                                                  
        scale_y_continuous(expand = c(0,0), labels = scales::comma)  +
        #scale_fill_gradient(low = "pink", high = "red4")
        scale_fill_gradientn(colours = c("#0092FFFF", "#00FF92FF", "#49FF00FF", "#FFDB00FF", "#FF0000FF"))
        p <- p + labs(title = plot_title, x = x_label, y = y_label)                                                      
        
        ggplot2::ggsave(file = outFile, device = "tiff", compression = "lzw", dpi = 150, width = 12, height = 9, units = "in")
        print(p)
    }
        
    ## end of no ribo ##########################################

    return()
}


#' Print expression stacked bar graph for genes/transcripts most preferentially expressed in selected
#'     putative CSCs.
#'
#' This function prints a stacked bar graph to display expression in the genes/transcripts most
#'     preferentially expressed in a selected set of putative CSCs. This may use expression in the
#'     form of read counts, TPM, TPM count fractions, Cells Per Million, etc. Cells Per Million or
#'     TPM count fractions are recommended for tumour cell populations containing a mix of cell types.
#'     Read counts are recommended for homogenous cell populations e.g. PDXs, CDXs. Each bar in the
#'     plot generated is subdivided by the proportion of the average expression that is contributed
#'     by each cell in the pCSC group. The bars are ordered by these average expression scores across
#'     the pCSC group. 
#' 
#' @param tExprMatr Expression matrix describing a cancer single-cell population; column
#'     names should be gene names, row names should be cancer cell names.
#' @param pCSC_group_positions Row number of CSCs in the CSC group in the expression matrix
#'     tExprMatr.
#' @param tumour_only_positions Positions (row numbers) in the expression matrix tExprMatr
#'     thought to be tumour cells e.g. when the population contains cells annotated or predicted
#'     as other cell types a subset may be defined, or if all cells are thought to be malignant
#'     cancer cells this may be all cells.
#' @param this_patient_tumour_only_positions Optional: Expression matrix tExprMatr positions
#'     (row numbers) for tumour cells from a specific patient i.e. that the CSCs for the
#'     preferential expression plots are taken from. So, the plot needs to fulfil the criteria
#'     that the CSCs are from only one patient, but the tumour cell population is from multiple
#'     patients.
#' @param pCSC_group_name Name of CSC group to print preferential expression plots for
#'     e.g. 'pCSC_group_1_SC4_Talazoparib.LB17011', but can be anything.
#' @param tissue_type Name of tissue type e.g. 'CDX'; can be anything. 
#' @param family_type_name Name of the gene group type; usually 'HUGO'.
#' @param patient_source_name Optional: Name of the patient sample that CSC cells are from;
#'     just name for in the plot title, can be anything.
#' @param separate_row_for_high_expression_transcripts Optional: Set to 1 to print very high
#'     expression transcripts or gene groups on a separate row in some of the bar charts; helpful
#'     when there is a big difference between the most highly expressed and those on the same
#'     bar chart row.
#' @param no_ribo Optional: Set to 1 to remove ribosomal transcripts from some of the bar charts
#' @param manually_set_max_number_of_graph_rows Optional: Set to manually override the
#'     maximum number of rows shown in the plot. 
#' @return Highest ranked preferentially expressed genes
#' @export
stackedBarchartsForSingleTranscripts <- function(tExprMatr,                           
                                                 pCSC_group_positions,                  
                                                 tumour_only_positions,                 
                                                 this_patient_tumour_only_positions = NULL,
                                                 pCSC_group_name,
                                                 tissue_type,
                                                 family_type_name,
                                                 patient_source_name,
                                                 separate_row_for_high_expression_transcripts = NULL,
                                                 no_ribo = NULL,
                                                 manually_set_max_number_of_graph_rows = NULL,
                                                 HUGO_abbreviations_and_fullnames = NULL,
                                                 HUGOFullNameFilesDir = NULL
                                                 ){
        
    #3# STACKED BARCHARTS FOR SINGLE TRANSCRIPTS
    cat("Printing: Single transcripts barchart (stacked)\n")

    if(is.null(separate_row_for_high_expression_transcripts)){
        separate_row_for_high_expression_transcripts <- FALSE
    }
    if(is.null(no_ribo)){
        no_ribo <- FALSE
    }
    if(!is.null(this_patient_tumour_only_positions)){
        if(length(this_patient_tumour_only_positions) == 1){
            if(this_patient_tumour_only_positions == ''){
                this_patient_tumour_only_positions <- NULL
            }
        }
    }
    
    #### Prep data

    pCSCs_tExprMatr                           <- tExprMatr[pCSC_group_positions,]

    pCSCs_tExprMatr_colmeans                  <- colMeans(pCSCs_tExprMatr)
    pCSCs_tExprMatr_colmeans_sorted_rev       <- rev(sort(pCSCs_tExprMatr_colmeans))
    pCSCs_tExprMatr_colmeans_sorted_rev_names <- names(pCSCs_tExprMatr_colmeans_sorted_rev)
    pCSCs_tExprMatr_colmeans_rev_order        <- rev(order(pCSCs_tExprMatr_colmeans))

    pCSCs_tExprMatr_colmeans_sorted_rev_full_names <- getGenesFullNames(pCSCs_tExprMatr_colmeans_sorted_rev_names, HUGO_abbreviations_and_fullnames, HUGOFullNameFilesDir)

    pCSCs_tExprMatr_rev_ordered_by_colmeans   <- pCSCs_tExprMatr[,pCSCs_tExprMatr_colmeans_rev_order]

    if(length(this_patient_tumour_only_positions) > 0){
        this_patient_tumour_only_tExprMatr <- tExprMatr[this_patient_tumour_only_positions,]
        this_patient_tumour_only_average_expression_for_each_transcript_vec <- round( ( (colSums(this_patient_tumour_only_tExprMatr)) / (nrow(this_patient_tumour_only_tExprMatr)) ) , 2)
        this_patient_tumour_only_average_expression_for_each_transcript_colmeans_rev_ordered_vec <- this_patient_tumour_only_average_expression_for_each_transcript_vec[pCSCs_tExprMatr_colmeans_rev_order]
    }

    columnNumbers <- 1:ncol(tExprMatr)
    numberOfColumns <- length(columnNumbers)

    all_samples_tumour_only_column_mean <- colMeans(tExprMatr[tumour_only_positions,])
    all_samples_column_mean             <- colMeans(tExprMatr)

    pCSC_group_tExprMatr_col_means              <- colMeans(tExprMatr[pCSC_group_positions, , drop = FALSE])
    pCSC_group_column_mean                      <- pCSC_group_tExprMatr_col_means[columnNumbers]
    
    all_transcripts_pCSC_group_expression_mean  <- round(sum(pCSC_group_column_mean)/numberOfColumns, 2)
    all_transcripts_tumour_only_expression_mean <- round(sum(all_samples_tumour_only_column_mean)/numberOfColumns, 2)
    all_transcripts_all_samples_expression_mean <- round(sum(all_samples_column_mean)/numberOfColumns, 2)
    
    ####
    
    pCSCs_tExprMatr_rev_ordered_by_colmeans_list <- list()
    pCSCs_tExprMatr_rev_ordered_by_colmeans_list_proportions <- list()

    # convert to list with one DF column in each element
    for(i in 1:ncol(pCSCs_tExprMatr_rev_ordered_by_colmeans)){
        pCSCs_tExprMatr_rev_ordered_by_colmeans_list[[(length(pCSCs_tExprMatr_rev_ordered_by_colmeans_list) + 1)]] <- pCSCs_tExprMatr_rev_ordered_by_colmeans[,i]
    }

    prDF5 <- data.frame(
        pCSCs_tExprMatr_rev_ordered_by_colmeans_list,
        stringsAsFactors = FALSE
    )
    colnames(prDF5) <- pCSCs_tExprMatr_colmeans_sorted_rev_full_names
    prDF5_averages <- prDF5 / nrow(prDF5)
    number_of_cells <- nrow(prDF5)

    prDF6 <- data.frame(
        cell                            = rep(rownames(prDF5), ncol(prDF5)),
        transcript_name                 = rep(colnames(prDF5), each=nrow(prDF5)),
        expression                      = as.vector(round(unlist(c(prDF5)), 2)),
        expression_average              = as.vector(round(unlist(c(prDF5_averages)), 2)),
        cell_number                     = rep(paste0("#", 1:nrow(prDF5)), ncol(prDF5)),
        expression_average_sum_for_calc = rep(-1, nrow(prDF5)),
        stringsAsFactors = TRUE
    )

    number_of_transcripts_to_show <- 360
    reorder_seq <- pCSCs_tExprMatr_colmeans_rev_order # RIBO EFFECT? YES, BUT CORRECTED BELOW.

    ## no ribo #################################################
    
    if(no_ribo == TRUE){
        prDF6_pos_to_remove_vec <- vector()
        for(i in 1:nrow(prDF6)){
            this_transcript_name <- as.vector(prDF6[i, 2])
            if(substr(this_transcript_name, 1, 9) == "ribosomal"){
                prDF6_pos_to_remove_vec <- c(prDF6_pos_to_remove_vec, i)
            }
        }
        if(length(prDF6_pos_to_remove_vec) > 0){
            prDF6_no_ribo <- prDF6[-c(prDF6_pos_to_remove_vec),]
        } else {
            prDF6_no_ribo <- prDF6
        }
        prDF6 <- prDF6_no_ribo
    
        pCSCs_tExprMatr_colmeans_rev_order_rows_to_remove_vec <- vector()
        for(i in 1:length(pCSCs_tExprMatr_colmeans_sorted_rev_full_names)){
            this_transcript_name <- as.vector(pCSCs_tExprMatr_colmeans_sorted_rev_full_names[i])
            if(substr(this_transcript_name, 1, 9) == "ribosomal"){
                pCSCs_tExprMatr_colmeans_rev_order_rows_to_remove_vec <- c(pCSCs_tExprMatr_colmeans_rev_order_rows_to_remove_vec, i)
            }
        }
        reorder_seq_no_ribo <- reorder_seq[-c(pCSCs_tExprMatr_colmeans_rev_order_rows_to_remove_vec)]
        reorder_seq <- reorder_seq_no_ribo

        if(length(this_patient_tumour_only_positions) > 0){
            this_patient_tumour_only_average_expression_for_each_transcript_colmeans_rev_ordered_vec_no_ribo <-
                this_patient_tumour_only_average_expression_for_each_transcript_colmeans_rev_ordered_vec[-pCSCs_tExprMatr_colmeans_rev_order_rows_to_remove_vec]
            
            this_patient_tumour_only_average_expression_for_each_transcript_colmeans_rev_ordered_vec <-
                this_patient_tumour_only_average_expression_for_each_transcript_colmeans_rev_ordered_vec_no_ribo
        }
    }
    
    ## end of no ribo ##########################################

    start_row <- 1
    end_row <- (number_of_transcripts_to_show * (number_of_cells))
    prDF6 <- prDF6[start_row:end_row,]

    # add the family expression sums to the 6th column of the df
    transcript_names <- prDF6[,2]
    trans_pos <- match(transcript_names, unique(transcript_names))
    unique_trans_pos <- unique(trans_pos)
    for(i in 1:length(unique_trans_pos)){
        this_trans_pos <- unique_trans_pos[i]
        this_trans_all_pos <- which(trans_pos %in% this_trans_pos)
        this_trans_sum_expression <- sum(prDF6[this_trans_all_pos, 4]) # /(length(this_trans_all_pos))
        prDF6[this_trans_all_pos,6] <- this_trans_sum_expression
    }

    # Add cell numbers to be displayed to the df - if <= 5% of the family sum, don't display
    cell_numbers_to_show_vec <- vector()
    for(i in 1:nrow(prDF6)){
        if(prDF6[i,4] <= ((prDF6[i,6]/100) * 5)){ # 5%
        #if(prDF6[i,4] <= ((prDF6[i,6]/100) * (100 / length(pCSC_group_positions)))){ # even share of 100
            cell_numbers_to_show_vec[i] <- ""
        } else {
            cell_numbers_to_show_vec[i] <- as.character(prDF6[i,5])
        }
    }
    prDF6_with_cells_to_show <- cbind(prDF6, cell_numbers_to_show_vec)
    colnames(prDF6_with_cells_to_show)[7] <- "cell_numbers_to_show"
    prDF6 <- prDF6_with_cells_to_show

    # Retain order - reorder the transcript_name factors so print order is correct
    prDF6o <- prDF6[order(match(prDF6[,2], prDF6[,2])),] # no actual reorder
    prDF6o$transcript_name <- factor(prDF6o$transcript_name, levels=unique(prDF6o$transcript_name))
    rownames(prDF6o) <- NULL

    # Order by cell IDs for printing
    prDF6o2 <- prDF6o[order(match(prDF6o[,1], prDF6o[,1])),]
    prDF6o2$cell <- factor(prDF6o2$cell, levels=unique(prDF6o2$cell))
    rownames(prDF6o2) <- NULL

    # Reorder the cell_number factors so print order is correct    
    prDF6o3 <- prDF6o2[order(match(prDF6o2[,5], prDF6o2[,5])),] # no actual reorder
    prDF6o3$cell_number <- factor(prDF6o3$cell_number, levels=unique(prDF6o3$cell_number))
    rownames(prDF6o3) <- NULL

    ## RANKINGS - Store the result ranks for later analysis
    
    number_of_ranks_to_get <- 100
    if(length(pCSCs_tExprMatr_colmeans_sorted_rev) < number_of_ranks_to_get){ cat("\nFewer genes than ranks\n") }
    ranks1_STACKED_BARCHARTS_FOR_SINGLE_TRANSCRIPTS <- number_of_ranks_to_get:1
    names(ranks1_STACKED_BARCHARTS_FOR_SINGLE_TRANSCRIPTS) <- names(pCSCs_tExprMatr_colmeans_sorted_rev)[1:number_of_ranks_to_get]
    #names(ranks1_STACKED_BARCHARTS_FOR_SINGLE_TRANSCRIPTS) <- prDF6o2[1:number_of_ranks_to_get,2] # full names and includes no ribo if selected

    ##

    df_to_use <- prDF6o3
    prDF6o <- NULL
    prDF6o2 <- NULL 
    prDF6o3 <- NULL 

    ## Format transcript names for printing
    
    transcript_names_to_show <- vector()
    
    for(j in 1:nrow(df_to_use)){
        this_transcript_name <- as.vector(df_to_use[j,2])
        #this_transcript_name <- "averyveryveryverylongword andveryveryveryveryveryvery anotherveryveryveryverylongword"
	split_pos <- 0
        under_pos <- gregexpr(pattern = "_", this_transcript_name)
        space_pos <- gregexpr(pattern = " ", this_transcript_name)
        if(space_pos[[1]][1] == -1){
            this_transcript_name_len <- nchar(this_transcript_name)
            if(this_transcript_name_len >= 30){
	        split_point <- (nchar(this_transcript_name) %/% 2)
                this_transcript_name <- paste0(substr(this_transcript_name, 1, split_point), "-\n", substr(this_transcript_name, split_point+1, nchar(this_transcript_name)))
            }
        } else if(max(nchar(strsplit(this_transcript_name," ")[[1]])) > 30){
            split_pos <- which(nchar(strsplit(this_transcript_name, " ")[[1]]) > 30)
            for(k in 1:length(split_pos)){
                this_split_pos <- split_pos[k]
                word_start <- 0
                word_end <- 0
                split_point <- 0
                if(this_split_pos == 1){
                    word_start <- 1
                    word_end   <- space_pos[[1]][this_split_pos] - 1
                } else if(this_split_pos > length(space_pos[[1]])){
                    word_start <- space_pos[[1]][this_split_pos-1] + 1
                    word_end   <- nchar(this_transcript_name)
                } else {
                    word_start <- space_pos[[1]][this_split_pos-1] + 1
                    word_end   <- space_pos[[1]][this_split_pos] - 1
		}
		split_point <- ((word_start + word_end) %/% 2) + (k - 1)
                this_transcript_name <- paste0(substr(this_transcript_name, 1, split_point), "-\n", substr(this_transcript_name, split_point+1, nchar(this_transcript_name)))
            }
        }
        under_pos <- gregexpr(pattern = "_", this_transcript_name)
        space_pos <- gregexpr(pattern = " ", this_transcript_name)
        if((under_pos[[1]][1] > -1)){
            under_pos_len <- length(under_pos[[1]])
            words_on_a_row <- 4
            for(k in 1:under_pos_len){
                if((k %% words_on_a_row) == 0){
                    this_under_pos <- under_pos[[1]][k]
                    substr(this_transcript_name, this_under_pos, this_under_pos) <- "\n"
                }
            }
            this_transcript_name <- gsub("_", " ", this_transcript_name)
        } else if((space_pos[[1]][1] > -1)){
            space_pos_len <- length(space_pos[[1]])
            words_on_a_row <- 3 # split with \n if more than 3 words on a row
            for(k in 1:space_pos_len){
                if((k %% words_on_a_row) == 0){
                    this_space_pos <- space_pos[[1]][k]
                    substr(this_transcript_name, this_space_pos, this_space_pos) <- "\n"
                }
            }
        }
        transcript_names_to_show[j] <- this_transcript_name
    }

    # Add formatted names to the df
    df_to_use[,8] <- transcript_names_to_show
    colnames(df_to_use)[8] <- "transcript_name_to_show"
    df_to_use$transcript_name_to_show <- factor(df_to_use$transcript_name_to_show, levels = unique(df_to_use$transcript_name_to_show))
    
    prDF6o3_with_modified_transcript_names <- df_to_use

    ## Make a fake df to allow printing horizontal lines for means of all cells, tumour only, this patient tumour only
    
    # Cell column is fake i.e. just one abitrary cell label - be careful using this df
    fake_prDF6o3 <- prDF6o3_with_modified_transcript_names[1:number_of_transcripts_to_show,c(1,2,8,6)]

    all_samples_tumour_only_column_mean_reordered        <- all_samples_tumour_only_column_mean[reorder_seq] # previously: reorder_seq <- pCSCs_tExprMatr_colmeans_rev_order
    all_samples_tumour_only_column_mean_reordered_subset <- all_samples_tumour_only_column_mean_reordered[1:number_of_transcripts_to_show]
    all_samples_column_mean_reordered                    <- all_samples_column_mean[reorder_seq]
    all_samples_column_mean_reordered_subset             <- all_samples_column_mean_reordered[1:number_of_transcripts_to_show]

    fake_prDF6o3_with_means <- cbind(fake_prDF6o3, all_samples_tumour_only_column_mean_reordered_subset, all_samples_column_mean_reordered_subset)
    colnames(fake_prDF6o3_with_means)[5:6] <- c("tumour_only", "all_cells")

    if(length(this_patient_tumour_only_positions) > 0){
        fake_prDF6o3_with_means2 <- cbind(fake_prDF6o3_with_means, this_patient_tumour_only_average_expression_for_each_transcript_colmeans_rev_ordered_vec[1:number_of_transcripts_to_show]) # RIBO EFFECT (BUT FIXED ABOVE)
        colnames(fake_prDF6o3_with_means2)[7] <- c("this_patient_tumour_only")
        fake_prDF6o3_with_means <- fake_prDF6o3_with_means2
        fake_prDF6o3_with_means2 <- NULL
    }
    fake_prDF6o3 <- NULL

    ## PRINT THE PLOT

    #library("RColorBrewer")
    #require(scales)

    graph_length_limit <- 60

    p_list <- list()
    df_number_of_rows <- nrow(df_to_use)
    df_to_use_all <- df_to_use

    ## Separate high expression columns if requested
    #separate_row_for_high_expression_transcripts <- FALSE
    if(separate_row_for_high_expression_transcripts == TRUE){
        prDF6_expression_threshold <- 10000 # threshold to split on
        df_to_use_over_threshold_pos <- which(df_to_use_all[c(1:length(unique(df_to_use_all[,1]))),6] >= prDF6_expression_threshold)
        number_of_prDF6_high_expression_to_separate <- length(df_to_use_over_threshold_pos)
        if(number_of_prDF6_high_expression_to_separate > graph_length_limit){
            cat("More number_of_prDF6_high_expression_to_separate to separate than graph len\n")
        }
        if(number_of_prDF6_high_expression_to_separate > 0){
            number_of_graph_rows <- (((number_of_transcripts_to_show - 1) %/% graph_length_limit) + 2)
        } else {
            number_of_graph_rows <- (((number_of_transcripts_to_show - 1) %/% graph_length_limit) + 1)
        }
    } else {
        number_of_graph_rows <- (((number_of_transcripts_to_show - 1) %/% graph_length_limit) + 1)
    }

    if(!is.null(manually_set_max_number_of_graph_rows)){
        number_of_graph_rows <- manually_set_max_number_of_graph_rows
    }
    
    ## no ribo #################################################

    if(no_ribo == TRUE){
        tif_name <- paste0(pCSC_group_name, "_transcript_expression_bargraph-stacked-no_ribo.tif")
    } else {
        tif_name <- paste0(pCSC_group_name, "_transcript_expression_bargraph-stacked.tif")
    }

    ## end of no ribo ##########################################

    #outFile <- paste0(outDir, "preferential_expression/", tif_name)
    outFile <- paste0(outDir, tif_name)

    #number_of_graph_rows <- 1

    for(i in 1:number_of_graph_rows){
        #cat(paste0("=> Printing row ", i, "\n"))

        # Get the df row numbers for all transcripts on this row for all cells and subset the df with them
        if((separate_row_for_high_expression_transcripts == 1) && (number_of_prDF6_high_expression_to_separate > 0)){
            if(i == 1){
                df_first_transcript_name_pos <- 1
                df_last_transcript_name_pos  <- number_of_prDF6_high_expression_to_separate
            } else if (i == 2){
                df_first_transcript_name_pos <- number_of_prDF6_high_expression_to_separate + 1
                df_last_transcript_name_pos  <- graph_length_limit
            } else {
                df_first_transcript_name_pos <- ((graph_length_limit * ((i-1)-1)) + 1)
                df_last_transcript_name_pos  <- (graph_length_limit * (i-1))
            }
        } else {
            df_first_transcript_name_pos <- ((graph_length_limit * (i-1)) + 1)
            df_last_transcript_name_pos  <- (graph_length_limit * i)
        }
        if(df_last_transcript_name_pos > number_of_transcripts_to_show){
            df_last_transcript_name_pos <- number_of_transcripts_to_show
        }
        these_transcripts_pos <- c(df_first_transcript_name_pos:df_last_transcript_name_pos)
        these_transcripts_pos_converted_vec <- vector()
        for(j in 1:length(these_transcripts_pos)){
            this_transcript_pos <- these_transcripts_pos[j]
            this_transcript_pos_converted_vec <- vector()
            for(k in 1:number_of_cells){
                this_transcript_pos_converted_vec[k] <- ((k - 1) * number_of_transcripts_to_show) + this_transcript_pos
            }
            these_transcripts_pos_converted_vec <- c(these_transcripts_pos_converted_vec, this_transcript_pos_converted_vec)
        }
        these_transcripts_pos_converted <- sort(these_transcripts_pos_converted_vec)

        df_to_use_this_row_subset <- df_to_use_all[these_transcripts_pos_converted,]
        df_to_use <- df_to_use_this_row_subset #; df_to_use_this_row_subset <- NULL
        
        ##
        
        # Get all the positions in this df subset for the transcripts on this row and find the max out of the sums for all cells for each transcript i.e. max row height
        df_to_use_transcript_index <- match(df_to_use[,2], unique(df_to_use[,2]))
        unique_df_to_use_transcript_index <- unique(df_to_use_transcript_index)

        these_transcripts_sums_vec <- vector()
        for(j in 1:length(unique_df_to_use_transcript_index)){
            #j <- 1
            this_unique_df_to_use_transcript_index <- unique_df_to_use_transcript_index[j]    
            these_transcripts_sums_vec[j] <- sum(df_to_use[which(df_to_use_transcript_index %in% this_unique_df_to_use_transcript_index),4])
        }
        maxHeightThisRow <- max(these_transcripts_sums_vec)
        
        ##
        
        # Subset the fake df (contains means for horizontal lines) for the transcripts in this row
        these_fake_pos <- df_first_transcript_name_pos:df_last_transcript_name_pos
        this_fake_prDF6o3_with_means <- fake_prDF6o3_with_means[these_fake_pos,]

        ##
        
        # Print stuff
        
        x_label <- "Transcript"
        y_label <- frag_count_method
        plot_title <- ""
        if(i == 1){
            plot_title <- paste0(toupper(tissue_type), ": [", family_type_name, "] Transcript expression for ", pCSC_group_name, " (", length(pCSC_group_positions), " cells)", " from samples: ", patient_source_name)
        }
        axis_text_size <- 6
    
        getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel2"))
        
        if(i == 1){    
            p <- ggplot2::ggplot(data=df_to_use, aes(x = transcript_name_to_show, y = expression_average, fill = cell)) +
            geom_bar(stat="identity", color="black", size=0.1, width=1) +
            theme(
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),        
            panel.grid.minor = element_blank(),        
            panel.background = element_blank(),        
            panel.spacing = unit(0.1, "cm"),
            axis.text.x = element_text(color = "black", angle = 270, hjust = 0, vjust = 0.5, size = axis_text_size),
            axis.text.y = element_text(color = "black"),
            legend.text = element_text(size = 8),  
            legend.box = "horizontal",                                              
            legend.position = c(0.92, 0.85), 
            plot.title = element_text(hjust = 0, size = 10),
            plot.margin = ggplot2::margin(0.1, 0.4, 0.1, 0.1, "cm"))
        } else {
            p <- ggplot2::ggplot(data = df_to_use, aes(x = transcript_name_to_show, y = expression_average, fill = cell)) +
            geom_bar(stat = "identity", color = "black", size = 0.1, width = 1) +
            theme(
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),        
            panel.grid.minor = element_blank(),        
            panel.background = element_blank(),        
            panel.spacing = unit(0.1, "cm"),
            axis.text.x = element_text(color = "black", angle = 270, hjust = 0, vjust = 0.5, size = axis_text_size),
            axis.text.y = element_text(color = "black"),
            legend.position = "none",                                                                                        
            plot.title = element_text(hjust = 0, size = 10),
            plot.margin = ggplot2::margin(0.1, 0.4, 0.1, 0.1, "cm"))
        }
        p <- p +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0), labels = scales::comma) +
        #scale_fill_manual(values = getPalette(number_of_cells), guide = FALSE) +
        scale_fill_manual(values = getPalette(number_of_cells), guide = "none") +
        geom_hline(aes(yintercept = all_transcripts_pCSC_group_expression_mean, color = "pCSC group: all transcripts mean", linetype = "pCSC group: all transcripts mean"), size = 0.5, alpha = 0.7) + # black dotted
        #geom_hline(aes(yintercept = all_transcripts_all_samples_expression_mean), color = "blue", linetype = "dotted", size = 0.5, alpha = 0.7) # black dotted
        # Cell number label alternatives
        #geom_text(aes(label = cell_numbers_to_show, x=transcript_name_to_show), color = "black", position = position_stack(vjust = 0.5), size = 2, vjust = -0.5, hjust = 0.5) + # GOOD
        #geom_text(aes(label = cell_numbers_to_show, x=transcript_name_to_show), color = "black", position = position_stack(vjust = 0.5), size = 2, hjust = 0.5) + # GOOD
        geom_text(aes(label = cell_numbers_to_show, x=transcript_name_to_show), color = "black", position = position_stack(vjust = 0.5), size = 1.6, hjust = 0.5) # GOOD
        #options(ggrepel.max.overlaps = Inf)
        #ggrepel::geom_text_repel(aes(label = cell_numbers_to_show, x = transcript_name_to_show), color = "black", position = position_stack(vjust = 0.5), size = 2, 
        #                segment.color = 'black', segment.size = 0.1, force = 0.5, direction = "y") + #, nudge_x = -0.5) + # GOOD    
        
        #if(length(this_patient_tumour_only_positions) > 0){
        #if(length(which(this_patient_tumour_only_positions != tumour_only_positions)) > 0){
        if((length(this_patient_tumour_only_positions) > 0) && (!identical(sort(this_patient_tumour_only_positions), sort(tumour_only_positions)))){ # if this_patient_tumour_only_positions are not the same as tumour_only_positions
            p <- p +
                     geom_hline(aes(yintercept = all_transcripts_tumour_only_expression_mean, color = "All tumour cells (all patients): all transcripts mean", linetype = "All tumour cells (all patients): all transcripts mean"),
                                size = 0.5, alpha = 0.7) + # red dotted
                     geom_line(data = this_fake_prDF6o3_with_means, aes(x = transcript_name_to_show, y = tumour_only, color = "Tumour cells only (all patients) average", linetype = "Tumour cells only (all patients) average"), group = 1, size = 0.5, alpha = 0.5)  +
                     geom_line(data = this_fake_prDF6o3_with_means, aes(x = transcript_name_to_show, y = all_cells, color = "All cells (all patients) average", linetype = "All cells (all patients) average"), group = 1, size = 0.5, alpha = 0.5) +
    
                     geom_line(data = this_fake_prDF6o3_with_means, aes(x = transcript_name_to_show, y = this_patient_tumour_only, color = "pCSC group source patient all tumour cells average", linetype = "pCSC group source patient all tumour cells average"),
                               group = 1, size = 0.5, alpha = 0.5) + 
                     scale_color_manual(breaks = c("pCSC group source patient all tumour cells average", "Tumour cells only (all patients) average", "All cells (all patients) average", "pCSC group: all transcripts mean",
                                        "All tumour cells (all patients): all transcripts mean"),
                                        values = c("pCSC group source patient all tumour cells average"="darkgreen", "Tumour cells only (all patients) average"="red", "All cells (all patients) average"="blue", "pCSC group: all transcripts mean"="black",
                                        "All tumour cells (all patients): all transcripts mean"="red")) +
                     scale_linetype_manual(breaks = c("pCSC group source patient all tumour cells average", "Tumour cells only (all patients) average", "All cells (all patients) average", "pCSC group: all transcripts mean",
                                           "All tumour cells (all patients): all transcripts mean"),
                                           values = c("pCSC group source patient all tumour cells average"="solid", "Tumour cells only (all patients) average"="solid", "All cells (all patients) average"="solid", "pCSC group: all transcripts mean"="dotted",
                                           "All tumour cells (all patients): all transcripts mean"="dotted"))
        } else {
            p <- p +
                     geom_hline(aes(yintercept = all_transcripts_tumour_only_expression_mean, color = "All tumour cells: all transcripts mean", linetype = "All tumour cells: all transcripts mean"),
                                size = 0.5, alpha = 0.7) + # red dotted
                     geom_line(data = this_fake_prDF6o3_with_means, aes(x = transcript_name_to_show, y = tumour_only, color = "Tumour cells only average", linetype = "Tumour cells only average"), group = 1, size = 0.5, alpha = 0.5)  +
                     geom_line(data = this_fake_prDF6o3_with_means, aes(x = transcript_name_to_show, y = all_cells, color = "All cells average", linetype = "All cells average"), group = 1, size = 0.5, alpha = 0.5) +

                     scale_color_manual(breaks = c("Tumour cells only average", "All cells average", "pCSC group: all transcripts mean", "All tumour cells: all transcripts mean"),
                                        values = c("Tumour cells only average"="red", "All cells average"="blue", "pCSC group: all transcripts mean"="black", "All tumour cells: all transcripts mean"="red")) +
                     scale_linetype_manual(breaks = c("Tumour cells only average", "All cells average", "pCSC group: all transcripts mean", "All tumour cells: all transcripts mean"),
                                           values = c("Tumour cells only average"="solid", "All cells average"="solid", "pCSC group: all transcripts mean"="dotted", "All tumour cells: all transcripts mean"="dotted"))
        }
        p <- p + labs(title = plot_title, x = x_label, y = y_label, colour = "Cell Group Expression Averages", linetype = "Cell Group Expression Averages", size = "", fill = "", shape = "")

        p_list[[i]] <- p        
    }

    #library("cowplot")

    graph_row_lengths <- list()
    for(i in 1:length(p_list)){
	graph_row_lengths[[i]] <- length(unique(p_list[[i]][[1]][[2]]))
    }
    max_row_len <- max(unlist(graph_row_lengths))
    p_draw <- cowplot::ggdraw()
    p_list_len <- length(p_list)
    first_row_len <- max(unlist(graph_row_lengths))
    for(i in 1:length(p_list)){
        if(graph_row_lengths[[i]] < max_row_len){
            p_draw <- p_draw + cowplot::draw_plot(p_list[[i]], x = 0, y = (1 - ((1/p_list_len) * i)), width = (((1/max_row_len) * (graph_row_lengths[[i]])) + 0.05), height = (1/p_list_len)) 
        } else {
            p_draw <- p_draw + cowplot::draw_plot(p_list[[i]], x = 0, y = (1 - ((1/p_list_len) * i)), width = ((1/max_row_len) * (graph_row_lengths[[i]])), height = (1/p_list_len))
        }
    }

    dpiVal <- 150
    thisTifWidth <- ((max(unlist(graph_row_lengths))) / 3)  + 2
    if(thisTifWidth < 4){  thisTifWidth <- 4  } 
    if(thisTifWidth > 49){  thisTifWidth <- 49  }
    thisTifHeight <- (7 * number_of_graph_rows)
    #outFile <- paste0(outDir, "preferential_expression/", tif_name)
    outFile <- paste0(outDir, tif_name)
    if(thisTifHeight > 120){
        thisTifHeight <- 120
        dpiVal <- 100
    }

    tiff(filename = outFile, compression = "lzw", res = dpiVal, width = thisTifWidth, height = thisTifHeight, units = "in", pointsize = 12)
    print(p_draw)
    dev.off()    
    #ggplot2::ggsave(p_draw, file = outFile, device = "tiff", compression = "lzw", dpi = dpiVal, width = thisTifWidth, height = thisTifHeight, units = "in")
    #ggplot2::ggsave(p_draw, file = outFile, device = "jpg", compression = "lzw", dpi = dpiVal, width = thisTifWidth, height = thisTifHeight, units = "in")
    print(p_draw)

    return(ranks1_STACKED_BARCHARTS_FOR_SINGLE_TRANSCRIPTS)
}


#' Print expression stacked bar graph for HUGO gene groups most preferentially expressed in
#'     selected putative CSCs. 
#'
#' This function prints a stacked bar graph to display expression in HUGO gene groups most
#'     preferentially expressed in a selected set of putative CSCs. This may use expression in the
#'     form of read counts, TPM, TPM count fractions, Cells Per Million, etc. Cells Per Million or
#'     TPM count fractions are recommended for tumour cell populations containing a mix of cell types.
#'     Read counts are recommended for homogenous cell populations e.g. PDXs, CDXs. The height of
#'     each bar is the average expression value for the members of a HUGO gene group in the cells
#'     of the selected pCSC group. Each bar is subdivided by the proportion of the gene group
#'     average expression contributed by each gene; subdivisions are labelled with gene names, and
#'     where the subdivision is too small to fit the gene name in the graphic they are labelled
#'     nearby. Bars are ordered by the HUGO gene group average expression value for the pCSC group
#'     minus the gene group average expression in all tumour cells from all patients; the difference
#'     between the pCSC group average expression and all tumour cells average expression is shown by
#'     a pink dashed line. Other horizontal lines: the pCSC group’s source patient average expression
#'     value for all tumour cells = dashed green line; average expression value for tumour cells from
#'     all patients = red dashed line; all cells from all patients (e.g. including immune cells)
#'     average expression value = blue dashed line; and the average expression value for all
#'     transcripts in the pCSC group = black dotted line. N.B. not all genes officially defined in
#'     a HUGO gene group are present in the scRNA-Seq dataset, depending on the experimental platform
#'     used.
#'
#' @param tExprMatr Expression matrix describing a cancer single-cell population; column
#'     names should be gene names, row names should be cancer cell names.
#' @param pCSC_group_positions Row number of CSCs in the CSC group in the expression matrix
#'     tExprMatr.
#' @param tumour_only_positions Positions (row numbers) in the expression matrix tExprMatr
#'     thought to be tumour cells e.g. when the population contains cells annotated or predicted
#'     as other cell types a subset may be defined, or if all cells are thought to be malignant
#'     cancer cells this may be all cells.
#' @param this_patient_tumour_only_positions Optional: Expression matrix tExprMatr positions
#'     (row numbers) for tumour cells from a specific patient i.e. that the CSCs for the
#'     preferential expression plots are taken from. So, the plot needs to fulfil the criteria
#'     that the CSCs are from only one patient, but the tumour cell population is from multiple
#'     patients.
#' @param HUGO_groups_sets_of_columnNumbers List of lists, i.e. a list of HUGO gene groups,
#'     with each named gene group containing positions (row numbers) in the expression matrix
#'     tExprMatr of all the gene group members present in the expression matrix.
#' @param pCSC_group_name Name of CSC group to print preferential expression plots for
#'     e.g. 'pCSC_group_1_SC4_Talazoparib.LB17011', but can be anything.
#' @param outDir Output directory for plots.
#' @param tissue_type Name of tissue type e.g. 'CDX'; can be anything. 
#' @param family_type_name Name of the gene group type; usually 'HUGO'.
#' @param patient_source_name Optional: Name of the patient sample that CSC cells are from;
#'     just name for in the plot title, can be anything.
#' @param frag_count_method Expression measurement type, e.g. 'raw'; but just a label, can be anything.
#' @param separate_row_for_high_expression_transcripts Optional: Set to 1 to print very high
#'     expression transcripts or gene groups on a separate row in some of the bar charts; helpful
#'     when there is a big difference between the most highly expressed and those on the same
#'     bar chart row.
#' @param manually_set_max_number_of_graph_rows Optional: Set to manually override the
#'     maximum number of rows shown in the stacked family plot. 
#' @return List of highest ranked preferentially expressed genes and gene families:
#'     - ranked families from family barcharts;
#'     - ranked transcripts from family barcharts.
#' @export
stackedFamilyExpressionBarchartSingleBars <- function(tExprMatr,                                            
                                                      pCSC_group_positions,
                                                      tumour_only_positions,
                                                      this_patient_tumour_only_positions,
                                                      HUGO_groups_sets_of_columnNumbers,
                                                      pCSC_group_name,                             
                                                      outDir,                                      
                                                      tissue_type,                                 
                                                      family_type_name,                            
                                                      patient_source_name,
                                                      frag_count_method,                           
                                                      separate_row_for_high_expression_transcripts = NULL,
                                                      manually_set_max_number_of_graph_rows = NULL){
    
    #4# STACKED FAMILY EXPRESSION BARCHART SINGLE BARS
    cat("Printing: Family transcripts barchart (stacked)\n")
    
    if(is.null(separate_row_for_high_expression_transcripts)){
        separate_row_for_high_expression_transcripts <- FALSE
    }
    if(!is.null(this_patient_tumour_only_positions)){
        if(length(this_patient_tumour_only_positions) == 1){
            if(this_patient_tumour_only_positions == ''){
                this_patient_tumour_only_positions <- NULL
            }
        }
    }

    #### Prep data

    columnNumbers <- 1:ncol(tExprMatr)
    numberOfColumns <- length(columnNumbers)

    all_transcript_names <- colnames(tExprMatr)
    HUGO_column_groups_names <- names(HUGO_groups_sets_of_columnNumbers)

    pCSCs_tExprMatr                    <- tExprMatr[pCSC_group_positions,]
    tumour_only_tExprMatr              <- tExprMatr[tumour_only_positions,]
    if(length(this_patient_tumour_only_positions) > 0){
        this_patient_tumour_only_tExprMatr <- tExprMatr[this_patient_tumour_only_positions,]
    }
    
    all_samples_tumour_only_column_mean <- colMeans(tExprMatr[tumour_only_positions,])
    all_samples_column_mean             <- colMeans(tExprMatr)
    
    pCSC_group_tExprMatr_col_means              <- colMeans(tExprMatr[pCSC_group_positions, , drop = FALSE])
    pCSC_group_column_mean                      <- pCSC_group_tExprMatr_col_means[columnNumbers]
    
    all_transcripts_pCSC_group_expression_mean  <- round(sum(pCSC_group_column_mean)/numberOfColumns, 2)
    all_transcripts_tumour_only_expression_mean <- round(sum(all_samples_tumour_only_column_mean)/numberOfColumns, 2)
    all_transcripts_all_samples_expression_mean <- round(sum(all_samples_column_mean)/numberOfColumns, 2)

    ####
    
    ## Prep Data for Stacked (Single-Cell) Barcharts
    ## -> collect expression level & frequency info for HUGO gene groups & store in DFs

    pCSC_group_average_expression_across_each_family                         <- vector()
    pCSC_group_average_expression_each_family_transcript_list                <- list()
    pCSC_group_average_expression_each_family_transcript_percentage_max_list <- list()

    tumour_only_total_expression_across_each_family   <- vector()
    tumour_only_average_expression_across_each_family <- vector()
    
    #this_patient_tumour_only_total_expression_across_each_family   <- vector()
    this_patient_tumour_only_average_expression_across_each_family <- vector()
    
    all_samples_all_cells_total_expression_across_each_family   <- vector()
    all_samples_all_cells_average_expression_across_each_family <- vector()

    # Collect info for each HUGO gene group (for gene family expression DF, pCS group subset, tumour only subset, this patient tumour only):
    #   - average_expression_across_this_family                  = average expression for all matrix cells of family expression DF and subsets of it i.e. pCSC cells subset, tumour only subset, this patient tumour only 
    #   - average_expression_this_family_transcripts             = average expression over each matrix column (gene) for family expression DF and subsets of it
    #   - average_expression_this_family_transcripts_percentages = percentage (for each transcript) of sum of average_expression_this_family_transcripts for this family
    #   - total_expression_across_this_family                    = sum of family expression DF

    for(k in 1:length(HUGO_groups_sets_of_columnNumbers)){ #  each transcript family i.e. HUGO gene group
        if(length(HUGO_groups_sets_of_columnNumbers[[k]]) == 1){
            cat("Too few columns\n")
        }

        ## numberOfCols == number of transcripts in this gene group
        numberOfCols <- length(HUGO_groups_sets_of_columnNumbers[[k]])
        
        pCSC_group_average_expression_across_this_family <- round(((sum(pCSCs_tExprMatr[, HUGO_groups_sets_of_columnNumbers[[k]]])) / (length(pCSC_group_positions) * numberOfCols)), 2)
        pCSC_group_average_expression_this_family_transcripts <- vector()
        pCSC_group_average_expression_this_family_transcripts <- colSums(pCSCs_tExprMatr[, HUGO_groups_sets_of_columnNumbers[[k]]]) / nrow(pCSCs_tExprMatr)
        pCSC_group_average_expression_this_family_transcripts_sum <- sum(pCSC_group_average_expression_this_family_transcripts)
        if(pCSC_group_average_expression_this_family_transcripts_sum == 0){
            pCSC_group_average_expression_this_family_transcripts_percentages_max <- 0
        } else {
            pCSC_group_average_expression_this_family_transcripts_percentages <- 100 * (pCSC_group_average_expression_this_family_transcripts / pCSC_group_average_expression_this_family_transcripts_sum)
            pCSC_group_average_expression_this_family_transcripts_percentages_max <- max(pCSC_group_average_expression_this_family_transcripts_percentages)
        }
        pCSC_group_average_expression_each_family_transcript_percentage_max_list[k] <- pCSC_group_average_expression_this_family_transcripts_percentages_max
        names(pCSC_group_average_expression_each_family_transcript_percentage_max_list)[k] <- names(HUGO_groups_sets_of_columnNumbers)[k]
        
        tumour_only_total_expression_across_this_family   <- sum(tumour_only_tExprMatr[, HUGO_groups_sets_of_columnNumbers[[k]]])
        tumour_only_average_expression_across_this_family <- round(((sum(tumour_only_tExprMatr[, HUGO_groups_sets_of_columnNumbers[[k]]])) / (nrow(tumour_only_tExprMatr) * numberOfCols)), 2)
        
        if(length(this_patient_tumour_only_positions) > 0){
            this_patient_tumour_only_total_expression_across_this_family   <- sum(this_patient_tumour_only_tExprMatr[, HUGO_groups_sets_of_columnNumbers[[k]]])
            this_patient_tumour_only_average_expression_across_this_family <- round(((sum(this_patient_tumour_only_tExprMatr[, HUGO_groups_sets_of_columnNumbers[[k]]])) / (nrow(this_patient_tumour_only_tExprMatr) * numberOfCols)), 2)
        }        

        all_cells_total_expression_across_this_family   <-         sum(tExprMatr[, HUGO_groups_sets_of_columnNumbers[[k]]])
        all_cells_average_expression_across_this_family <- round(((sum(tExprMatr[, HUGO_groups_sets_of_columnNumbers[[k]]])) / (nrow(tExprMatr) * numberOfCols)), 2)
    
        pCSC_group_average_expression_across_each_family[[k]]        <- pCSC_group_average_expression_across_this_family
        names(pCSC_group_average_expression_across_each_family[[k]]) <- HUGO_column_groups_names[k]

        pCSC_group_average_expression_each_family_transcript_list[[k]]      <- pCSC_group_average_expression_this_family_transcripts
        names(pCSC_group_average_expression_each_family_transcript_list)[k] <- HUGO_column_groups_names[k]

        tumour_only_total_expression_across_each_family[[k]]   <- tumour_only_total_expression_across_this_family
        tumour_only_average_expression_across_each_family[[k]] <- tumour_only_average_expression_across_this_family
        names(tumour_only_total_expression_across_each_family[[k]])   <- HUGO_column_groups_names[k]
        names(tumour_only_average_expression_across_each_family[[k]]) <- HUGO_column_groups_names[k]

        if(length(this_patient_tumour_only_positions) > 0){
            #this_patient_tumour_only_total_expression_across_each_family[[k]]          <- this_patient_tumour_only_total_expression_across_this_family
            this_patient_tumour_only_average_expression_across_each_family[[k]]        <- this_patient_tumour_only_average_expression_across_this_family
            #names(this_patient_tumour_only_total_expression_across_each_family[[k]])   <- HUGO_column_groups_names[k]
            names(this_patient_tumour_only_average_expression_across_each_family[[k]]) <- HUGO_column_groups_names[k]
        }

        all_samples_all_cells_total_expression_across_each_family[[k]]   <- all_cells_total_expression_across_this_family
        all_samples_all_cells_average_expression_across_each_family[[k]] <- all_cells_average_expression_across_this_family
        names(all_samples_all_cells_total_expression_across_each_family[[k]])   <- HUGO_column_groups_names[k]
        names(all_samples_all_cells_average_expression_across_each_family[[k]]) <- HUGO_column_groups_names[k]    
    }
            
    pCSC_group_all_cells_list_total_expression_across_each_family   <- list()
    pCSC_group_all_cells_list_average_expression_across_each_family <- list()

    ## Collect expression info for gene groups for individual cells:
    ##   - total family expression
    ##   - average family expression
    for(l in 1:length(pCSC_group_positions)){
        #l <- 1
        this_pCSC_group_position <- pCSC_group_positions[l]

        this_cell_total_expression_across_each_family   <- vector()
        this_cell_average_expression_across_each_family <- vector()

        for(m in 1:length(HUGO_groups_sets_of_columnNumbers)){ # each transcript family
            #m <- 1
            numberOfCols <- length(HUGO_groups_sets_of_columnNumbers[[m]])

            this_cell_total_expression_across_this_family    <- NULL
            this_cell_average_expression_across_this_family  <- NULL

            this_cell_total_expression_across_this_family    <- sum(tExprMatr[this_pCSC_group_position, HUGO_groups_sets_of_columnNumbers[[m]]])
            this_cell_average_expression_across_this_family  <- round(((this_cell_total_expression_across_this_family) / numberOfCols), 2)

            this_cell_total_expression_across_each_family[[m]]   <- this_cell_total_expression_across_this_family  
            this_cell_average_expression_across_each_family[[m]] <- this_cell_average_expression_across_this_family
            names(this_cell_total_expression_across_each_family[[m]])   <- HUGO_column_groups_names[m]
            names(this_cell_average_expression_across_each_family[[m]]) <- HUGO_column_groups_names[m]
        }
    
        pCSC_group_all_cells_list_total_expression_across_each_family[[l]]   <- this_cell_total_expression_across_each_family
        pCSC_group_all_cells_list_average_expression_across_each_family[[l]] <- this_cell_average_expression_across_each_family
    }

    # Store info in DFs for use in plots
    prDF4 <- data.frame(
        pCSC_group_expression_average                             = pCSC_group_average_expression_across_each_family,
        tumour_only_expression_average                            = tumour_only_average_expression_across_each_family,
        all_samples_expression_average                            = all_samples_all_cells_average_expression_across_each_family,
        pCSC_group_expression_average_each_transcript             = I(pCSC_group_average_expression_each_family_transcript_list),
        pCSC_group_max_transcript_percentage_expression_in_family = as.vector(unlist(pCSC_group_average_expression_each_family_transcript_percentage_max_list)),
        pCSC_group_expression_average_per_cell                    = pCSC_group_all_cells_list_average_expression_across_each_family,
        stringsAsFactors = FALSE
    )    

    rownames(prDF4) <- HUGO_column_groups_names
    colnames(prDF4) <- c(
                       "pCSC_group_expression_average",
                       "tumour_only_expression_average",
    	               "all_samples_expression_average",
                       "pCSC_group_expression_average_each_transcript",    
                       "pCSC_group_max_transcript_percentage_expression_in_family",
    	               paste0("pCSC_group_expression_average_cell_", rownames(tExprMatr[pCSC_group_positions,]))
                      )    

    if(length(this_patient_tumour_only_positions) > 0){
        prDF4_with_patient_tumour_info <- data.frame(
            pCSC_group_expression_average                             = pCSC_group_average_expression_across_each_family,
            tumour_only_expression_average                            = tumour_only_average_expression_across_each_family,
            all_samples_expression_average                            = all_samples_all_cells_average_expression_across_each_family,
            pCSC_group_expression_average_each_transcript             = I(pCSC_group_average_expression_each_family_transcript_list),
            pCSC_group_max_transcript_percentage_expression_in_family = as.vector(unlist(pCSC_group_average_expression_each_family_transcript_percentage_max_list)),
            this_patient_tumour_only_expression_average               = this_patient_tumour_only_average_expression_across_each_family,
            pCSC_group_expression_average_per_cell                    = pCSC_group_all_cells_list_average_expression_across_each_family,   
            stringsAsFactors = FALSE
        )

        rownames(prDF4_with_patient_tumour_info) <- HUGO_column_groups_names
        colnames(prDF4_with_patient_tumour_info) <- c(
            "pCSC_group_expression_average",
            "tumour_only_expression_average",
            "all_samples_expression_average",
            "pCSC_group_expression_average_each_transcript",    
            "pCSC_group_max_transcript_percentage_expression_in_family",
            "this_patient_tumour_only_expression_average",
            paste0("pCSC_group_expression_average_cell_", rownames(tExprMatr[pCSC_group_positions,]))
        )			  
    }

    expression_rev_order <- rev(order(pCSC_group_average_expression_across_each_family - tumour_only_average_expression_across_each_family))
    prDF4ro2  <- prDF4[expression_rev_order,]
    prDF4_expression_diff_ordered <- prDF4ro2 # prDF4 ordered by decreasing difference between pCSC_group_average_expression_across_each_family and tumour_only_average_expression_across_each_family)
    
    if(length(this_patient_tumour_only_positions) > 0){
        prDF4ro2_w_patient_t_info  <- prDF4_with_patient_tumour_info[expression_rev_order,]    
        prDF4_with_patient_tumour_info_expression_diff_ordered <- prDF4ro2_w_patient_t_info#; prDF4ro2_w_patient_t_info <- NULL
    }

    ###########################################
    
    ## Make a df with family expression info ordered by pCSC_group_average_expression_across_each_family - tumour_only_average_expression_across_each_family
    ## but remove families where all transcription is from one gene
    
    #cat("Printing: Gene family expression barchart (stacked)\n")

    prDF4 <- prDF4_expression_diff_ordered # Can order it differently here if needed

    # Remove groups where one transcript is 100% i.e. the whole family's expression
    prDF4_to_remove_vec <- vector()
    for(i in 1:nrow(prDF4)){
        if(prDF4$pCSC_group_max_transcript_percentage_expression_in_family[i] == 100){
        #if(prDF4$pCSC_group_max_transcript_percentage_expression_in_family[i] > 99){
            prDF4_to_remove_vec <- c(prDF4_to_remove_vec, i)
        }
    }
    
    if(length(prDF4_to_remove_vec) > 0){
        prDF4_no_removals <- prDF4
        prDF4_minus_over_transcription_threshold <- prDF4[-c(prDF4_to_remove_vec),]
        prDF4 <- prDF4_minus_over_transcription_threshold
        prDF4_minus_over_transcription_threshold <- NULL

        if(length(this_patient_tumour_only_positions) > 0){
            ## diff_ordered is from rev(order(pCSC_group_average_expression_across_each_family - tumour_only_average_expression_across_each_family))
            prDF4_with_patient_tumour_info_expression_diff_ordered_no_removals <- prDF4_with_patient_tumour_info_expression_diff_ordered
            prDF4_with_patient_tumour_info_expression_diff_ordered_minus_over_transcription_threshold <- prDF4_with_patient_tumour_info_expression_diff_ordered[-c(prDF4_to_remove_vec),]
            prDF4_with_patient_tumour_info_expression_diff_ordered <- prDF4_with_patient_tumour_info_expression_diff_ordered_minus_over_transcription_threshold
        }    
    }

    # Collect the transcription info for each gene family and put in a df
    
    pCSC_group_all_transcripts_family_name_vec                   <- vector()
    pCSC_group_all_transcripts_names_vec                         <- vector()
    pCSC_group_all_transcripts_names_and_family_vec              <- vector()
    pCSC_group_all_transcripts_expression_vec                    <- vector()
    pCSC_group_all_transcripts_expression_average_proportion_vec <- vector()

    pCSC_group_all_transcripts_expression_info <- prDF4$pCSC_group_expression_average_each_transcript # prDF4[,7]
    # Sort the transcripts within the families
    pCSC_group_all_transcripts_expression_info_sorted <- lapply(pCSC_group_all_transcripts_expression_info, sort, decreasing=F)

    k <- 1
    for(i in 1:length(pCSC_group_all_transcripts_expression_info_sorted)){
        this_family_len <- length(pCSC_group_all_transcripts_expression_info_sorted[[i]])
        for(j in 1:this_family_len){
            pCSC_group_all_transcripts_family_name_vec[k]                   <- names(pCSC_group_all_transcripts_expression_info_sorted)[[i]]
            pCSC_group_all_transcripts_names_vec[k]                         <- names(pCSC_group_all_transcripts_expression_info_sorted[[i]][j])
            pCSC_group_all_transcripts_names_and_family_vec[k]              <- paste0(names(pCSC_group_all_transcripts_expression_info_sorted[[i]][j]), "[", names(pCSC_group_all_transcripts_expression_info_sorted)[[i]], "]") 
            pCSC_group_all_transcripts_expression_vec[k]                    <- pCSC_group_all_transcripts_expression_info_sorted[[i]][j]
            pCSC_group_all_transcripts_expression_average_proportion_vec[k] <- (pCSC_group_all_transcripts_expression_info_sorted[[i]][j] / this_family_len) # NOT REALLY AVERAGE PROPORTION; JUST EXPRESSION DIVIDED BY FAMILY SIZE
            k <- (k + 1)
        }
    }
    
    number_of_families    <- length(unique(pCSC_group_all_transcripts_family_name_vec))
    number_of_transcripts <- length(pCSC_group_all_transcripts_family_name_vec)

    prDF7 <- data.frame(
        cell_group_type               = c(rep("pCSC group", length(pCSC_group_all_transcripts_family_name_vec))),
        family                        = pCSC_group_all_transcripts_family_name_vec,
        transcript                    = pCSC_group_all_transcripts_names_vec,
        transcript_and_family         = pCSC_group_all_transcripts_names_and_family_vec,
        expression                    = pCSC_group_all_transcripts_expression_vec,
        expression_average_proportion = pCSC_group_all_transcripts_expression_average_proportion_vec,
        stringsAsFactors = FALSE
    )
    df_to_use <- prDF7

    ## RANKINGS

    #2) STACKED FAMILY EXPRESSION BARCHART SINGLE BARS - rank top 100 families also
    #3)                                                - rank top 100 but only for those >5% of the family
    
    number_of_ranks_to_get <- 100
    expression_percentage_threshold <- 0.01

    ## Simple ranks of all families
    ranks2_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__FAMILIES        <- number_of_ranks_to_get:1
    names(ranks2_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__FAMILIES) <- rownames(prDF4)[1:number_of_ranks_to_get]
    
    ## Rank genes using their occurrences in the top families
    
    ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION <- vector()
    options(scipen = 999)
    
    unique_family_names <- unique(prDF7[,2])[1:number_of_ranks_to_get]
    unique_family_have_match_pos <- which(prDF7[,2] %in% unique_family_names)
    these_family_pos <- vector()
    for(i in 1:length(unique_family_have_match_pos)){
        this_pos <- unique_family_have_match_pos[i]
        these_family_pos[i] <- match(prDF7[this_pos,2], unique_family_names)
    }
    for(i in 1:length(these_family_pos)){
        this_family_number <- these_family_pos[i]
        this_family_pos <- which(these_family_pos %in% this_family_number)
        this_family_expression_sum <- sum(prDF7[this_family_pos,5])
    
        if(this_family_expression_sum > 0){
            this_percentage <- (prDF7[i,5] / this_family_expression_sum)
        } else {
            this_percentage <- 0
        }
    
        ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION[i] <- ((100 - this_family_number) + 1) * this_percentage # So family 1: each transcript's share of 100; family 2: share of 99, family 3, share of 98 ..
        names(ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION)[i] <- prDF7[i,3]
    }
    
    # Tidy up the gene ranks to only take the larger rank value
    ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted <- rev(sort(ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION))
    ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted_no_duplicates <-
        ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted[match(unique(names(ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted)),
        names(ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted))]
    ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted_no_duplicates_no_zeros <-
        ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted_no_duplicates[-c(which(ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted_no_duplicates == 0))]
    
    # Top genes using family ranks:
    #cat("\nTop genes using family ranks:\n")
    #print(ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted_no_duplicates_no_zeros[1:10])
    #cat('\n')

    ############################################################

    ## Get family transcript info, format transcript names for printing, & add to df

    family_names_to_show <- vector()
    transcript_info_to_show <- vector()

    for(j in 1:nrow(df_to_use)){
        this_family_name <- as.vector(df_to_use[j,2])
        this_family_size <- length(which(df_to_use[,2] %in% this_family_name))

        # Get the gene symbols and values for each family
        family_start_pos <- unique(match(df_to_use[,2], df_to_use[,2]))
        if(j %in% family_start_pos){ # only do for the first position for each family
            this_family_name_hugo_pos <- which(names(HUGO_groups_sets_of_columnNumbers) %in% as.vector(df_to_use[j,2]))
            this_family_size <- length(unlist(HUGO_groups_sets_of_columnNumbers[this_family_name_hugo_pos]))
            this_family_gene_symbols <- all_transcript_names[unlist(HUGO_groups_sets_of_columnNumbers[this_family_name_hugo_pos])]
            this_family_expression_means_each_transcript_average_across_pCSC_group <- as.vector(round(colMeans(pCSCs_tExprMatr[,unlist(HUGO_groups_sets_of_columnNumbers[this_family_name_hugo_pos])]), 1))
            these_gene_symbols_and_vals <- NULL
            this_order <- order(this_family_expression_means_each_transcript_average_across_pCSC_group)
            for(k in length(this_family_gene_symbols):1){ # Reverse order
                these_gene_symbols_and_vals <- paste0(these_gene_symbols_and_vals, this_family_gene_symbols[this_order[k]], "=", this_family_expression_means_each_transcript_average_across_pCSC_group[this_order[k]], "\n")
            }
            transcript_info_to_show[j] <- these_gene_symbols_and_vals
        } else {
            transcript_info_to_show[j] <- ""
        }
        ## this_family_name <- "averyveryveryverylongword andveryveryveryveryveryvery anotherveryveryveryverylongword"
	split_pos <- 0
        under_pos <- gregexpr(pattern = "_", this_family_name)
        space_pos <- gregexpr(pattern = " ", this_family_name)
        if(space_pos[[1]][1] == -1){
            this_family_name_len <- nchar(this_family_name)
            if(this_family_name_len >= 20){
	        split_point <- (nchar(this_family_name) %/% 2)
                this_family_name <- paste0(substr(this_family_name, 1, split_point), "-\n", substr(this_family_name, split_point+1, nchar(this_family_name)))
            }
        } else if(max(nchar(strsplit(this_family_name," ")[[1]])) > 20){
            split_pos <- which(nchar(strsplit(this_family_name, " ")[[1]]) > 20)
            for(k in 1:length(split_pos)){
                this_split_pos <- split_pos[k]
                word_start <- 0
                word_end <- 0
                split_point <- 0
                if(this_split_pos == 1){
                    word_start <- 1
                    word_end   <- space_pos[[1]][this_split_pos] - 1
                } else if(this_split_pos > length(space_pos[[1]])){
                    word_start <- space_pos[[1]][this_split_pos-1] + 1
                    word_end   <- nchar(this_family_name)
                } else {
                    word_start <- space_pos[[1]][this_split_pos-1] + 1
                    word_end   <- space_pos[[1]][this_split_pos] - 1
		}
		split_point <- ((word_start + word_end) %/% 2) + (k - 1)
                this_family_name <- paste0(substr(this_family_name, 1, split_point), "-\n", substr(this_family_name, split_point+1, nchar(this_family_name)))
            }
        }
        under_pos <- gregexpr(pattern = "_", this_family_name)
        space_pos <- gregexpr(pattern = " ", this_family_name)
        if((under_pos[[1]][1] > -1)){
            under_pos_len <- length(under_pos[[1]])
            words_on_a_row <- 4
            for(k in 1:under_pos_len){
                if((k %% words_on_a_row) == 0){
                    this_under_pos <- under_pos[[1]][k]
                    substr(this_family_name, this_under_pos, this_under_pos) <- "\n"
                }
            }
            this_family_name <- gsub("_", " ", this_family_name)
        } else if((space_pos[[1]][1] > -1)){
            space_pos_len <- length(space_pos[[1]])
            words_on_a_row <- 2
            for(k in 1:space_pos_len){
                if((k %% words_on_a_row) == 0){
                    this_space_pos <- space_pos[[1]][k]
                    substr(this_family_name, this_space_pos, this_space_pos) <- "\n"
                }
            }
        }
        this_family_name <- paste0(this_family_name, "\n[", this_family_size, " genes]")
        family_names_to_show[j] <- this_family_name
    }

    df_to_use[,2] <- family_names_to_show
    df_to_use[,7] <- as.factor(transcript_info_to_show)
    colnames(df_to_use)[7] <- "transcript_info_to_show"

    prDF8 <- df_to_use 

    # Convert families to factor for printing
    prDF8o <- prDF8[order(match(prDF8[,2], prDF8[,2])),]
    prDF8o$family <- factor(prDF8o$family, levels=unique(prDF8o$family))
    rownames(prDF8o) <- NULL

    # Convert transcript_and_family to factor for printing
    prDF8o2 <- prDF8o[order(match(prDF8o[,2], prDF8o[,2])),]
    prDF8o2$transcript_and_family <- factor(prDF8o2$transcript_and_family, levels=unique(prDF8o2$transcript_and_family))
    rownames(prDF8o2) <- NULL

    df_to_use <- prDF8o2

    ## Make a fake df for help in printing

    # prDF9 Has subset of data (one transcript per family) in column 2 to allow printing but non-sensical -> be careful using prDF9 column2, as only one gene name (first) out of the family
    prDF9_fake_col2_for_printing <- df_to_use[match(unique(df_to_use[,2]), df_to_use[,2]),c(2,4)] # take the first row for each family
    prDF9_fake_col2_for_printing[,3] <- prDF4$tumour_only_expression_average
    prDF9_fake_col2_for_printing[,4] <- prDF4$all_samples_expression_average
    prDF9_fake_col2_for_printing[,5] <- (prDF4$pCSC_group_expression_average - prDF4$tumour_only_expression_average)
    colnames(prDF9_fake_col2_for_printing)[3:5] <- c("tumour_only_expression", "all_samples_expression", "pCSC_group_tumour_difference_expression")

    if(length(this_patient_tumour_only_positions) > 0){
        prDF9_fake_col2_for_printing[,6] <- prDF4_with_patient_tumour_info_expression_diff_ordered$this_patient_tumour_only_expression_average
        prDF9_fake_col2_for_printing[,7] <- (prDF4_with_patient_tumour_info_expression_diff_ordered$pCSC_group_expression_average - prDF4_with_patient_tumour_info_expression_diff_ordered$this_patient_tumour_only_expression_average)
        colnames(prDF9_fake_col2_for_printing)[6:7] <- c("this_tumour_only_expression", "pCSC_group_this_tumour_difference_expression")
    }

    ## Find where graph should end i.e. pCSC_group_versus_tumour_threshold_reached
    
    prDF10 <- data.frame(
        cell_group_type    = c(rep("pCSC_group_expression_average", 1+(nrow(prDF4))), rep("tumour_only_expression_average", 1+(nrow(prDF4))), rep("all_samples_expression_average", 1+(nrow(prDF4)))), 
        family             = c("All", rownames(prDF4), "All", rownames(prDF4), "All", rownames(prDF4)),
        expression_average = c(all_transcripts_pCSC_group_expression_mean, prDF4$pCSC_group_expression_average, all_transcripts_tumour_only_expression_mean, prDF4$tumour_only_expression_average,
                               all_transcripts_all_samples_expression_mean, prDF4$all_samples_expression_average), 
        stringsAsFactors = FALSE
    )

    # Reorder to group together families & convert to factor
    prDF10o <- prDF10[order(match(prDF10[,2], prDF10[,2])),]
    prDF10o$family <- factor(prDF10o$family, levels=unique(prDF10o$family))
    rownames(prDF10o) <- NULL

    graph_length_limit <- 40 # shouldn't it be 60?
    display_threshold <- 1.5
    max_number_of_transcript_labels_to_show <- 25
    max_rows_when_threshold_fails <- 6
    pCSC_group_versus_tumour_threshold_reached <- 0

    if(prDF10o[1,3] <= prDF10o[2,3]){ # If all_transcripts_pCSC_group_expression_mean is <= all_transcripts_tumour_only_expression_mean (colmeans for respective subsets then averaged over all genes)
        cat(paste0("pCSC_group_expression_average is >= tumour_only_expression_average so will do ", max_rows_when_threshold_fails, " rows\n"))
        pCSC_group_versus_tumour_threshold_reached <- graph_length_limit * max_rows_when_threshold_fails
        if(length(unique(prDF10o[,2])) < pCSC_group_versus_tumour_threshold_reached){ cat("ERROR, not enough families\n") }
    } else {
        # Check if threshold reached i.e. (this family pCSC_group_expression_average - this family tumour_only_expression_average) < (1.5 * (all_transcripts_pCSC_group_expression_mean - all_transcripts_tumour_only_expression_mean))
        #all_pCSC_group_minus_tumour_all_expression_vals_list <- list()
        for(i in 1:((length(prDF10o[,3])/3) - 1)){
            #i <- 1
            thisIndex <- (i * 3) + 1
            if((prDF10o[thisIndex,3] - prDF10o[(thisIndex + 1),3]) < ((prDF10o[1,3] - prDF10o[2,3]) * display_threshold)){
                pCSC_group_versus_tumour_threshold_reached <- (i - 1)
                break
            }
            #all_pCSC_group_minus_tumour_all_expression_vals_list[[i]] <- (prDF10o[thisIndex,3] - prDF10o[(thisIndex + 1),3])
        }
    }
    
    ##

    prDF8o2_family_first_pos <- match(unique(prDF8o2[,2]), prDF8o2[,2])
    prDF8o2_family_last_pos <- ((length(prDF8o2[,2]) - (match(unique(prDF8o2[,2]), rev(prDF8o2[,2])))) + 1)

    # Subset df for length of graph i.e. where threshold reached
    if(pCSC_group_versus_tumour_threshold_reached > 0){
        prDF8o2_pCSC_group_versus_tumour_threshold_reached_pos <- prDF8o2_family_last_pos[pCSC_group_versus_tumour_threshold_reached]
        prDF8o2_subset <- prDF8o2[1:prDF8o2_pCSC_group_versus_tumour_threshold_reached_pos,]
        number_of_graph_rows <- ((pCSC_group_versus_tumour_threshold_reached - 1) %/% graph_length_limit) + 1

	df_to_use <- prDF8o2_subset
    } else {
        number_of_column_group_names <- length(HUGO_column_groups_names)
        number_of_sets_of_columnNumbers <- length(HUGO_groups_sets_of_columnNumbers)
        if(number_of_column_group_names != number_of_sets_of_columnNumbers){
            cat("Mismatch in the number of groups")
        } else {
            number_of_groups <- number_of_column_group_names
        }
        number_of_graph_rows <- ((number_of_groups - 1) %/% graph_length_limit) + 1
    }
    effective_number_of_graph_rows <- (number_of_graph_rows)

    df_to_use_all <- df_to_use

    ## PRINT THE PLOT

    p_list <- list()
    df_number_of_rows <- nrow(df_to_use)
    max_graph_rows <- 6
    #max_graph_rows <- 1
    if(!is.null(manually_set_max_number_of_graph_rows)){
        max_graph_rows <- manually_set_max_number_of_graph_rows
    }

    # Separate high expression columns
    #separate_row_for_high_expression_transcripts <- FALSE
    if(separate_row_for_high_expression_transcripts == TRUE){
        prDF8o2_expression_threshold <- 1000
        these_family_pos <- match(df_to_use_all[,2], unique(df_to_use_all[,2]))
        unique_these_family_pos <- unique(these_family_pos)
        these_family_transcripts_sums_vec <- vector()
        for(j in 1:length(unique_these_family_pos)){
            #j <- 1
            this_unique_these_family_pos <- unique_these_family_pos[j]    
            these_family_transcripts_sums_vec[j] <- sum(df_to_use_all[which(these_family_pos %in% this_unique_these_family_pos),6]) # sum of expression_average_proportion i.e. family average
        }
        number_of_families_over_prDF8o2_threshold_pos <- which(these_family_transcripts_sums_vec >= prDF8o2_expression_threshold) # bars to separate
        number_of_prDF8o2_families_high_expression_to_separate <- length(number_of_families_over_prDF8o2_threshold_pos)
        if(number_of_prDF8o2_families_high_expression_to_separate > graph_length_limit){
            cat("More number_of_prDF8o2_families_high_expression_to_separate than graph_length_limit\n")
        }
        if(number_of_prDF8o2_families_high_expression_to_separate > 0){
            number_of_graph_rows <- (number_of_graph_rows + 1)
        }
    }
    
    ##
    
    # Subset df if too many rows
    if(effective_number_of_graph_rows > max_graph_rows){
        #cat(paste0("Too many graph rows (", effective_number_of_graph_rows, ") so changing to ", max_graph_rows, "\n"))
        number_of_graph_rows <- max_graph_rows
        df_to_use_max_subset <- df_to_use[1:(prDF8o2_family_last_pos[(graph_length_limit * max_graph_rows)]),]
        df_to_use_all <- df_to_use_max_subset
    }

    for(i in 1:number_of_graph_rows){
        #cat(paste0("=> Printing row ", i, "\n"))

        # Get the df row numbers for all families on this row for all cells and subset the df with them
        if((separate_row_for_high_expression_transcripts == TRUE) && (number_of_prDF8o2_families_high_expression_to_separate > 0)){
            if(i == 1){
                df_first_row <- prDF8o2_family_first_pos[1]
                df_last_row  <- prDF8o2_family_last_pos[number_of_prDF8o2_families_high_expression_to_separate]
            } else if (i == 2){
                df_first_row <- prDF8o2_family_first_pos[(number_of_prDF8o2_families_high_expression_to_separate + 1)]
                df_last_row  <- prDF8o2_family_last_pos[graph_length_limit]
            } else {
                df_first_row <- prDF8o2_family_first_pos[((graph_length_limit * ((i-1)-1)) + 1)]
                df_last_row  <- prDF8o2_family_last_pos[(graph_length_limit * (i-1))]
            }
        } else {
            df_first_row <- prDF8o2_family_first_pos[((graph_length_limit * (i-1)) + 1)]
            df_last_row  <- prDF8o2_family_last_pos[(graph_length_limit * i)]
        }
        if(df_last_row > df_number_of_rows){
            df_last_row <- df_number_of_rows
        }
        if(df_first_row > df_number_of_rows){
            cat("ERROR IN THE ROW NUMBER CALCULATIONS (TOO MANY GRAPH ROWS?)\n")
        }

        these_rows <- c(df_first_row:df_last_row)
        df_to_use_this_row_subset <- df_to_use_all[these_rows,]
        df_to_use <- df_to_use_this_row_subset #; df_to_use_this_row_subset <- NULL

        ##
        
        # Get the average expression for these families
        these_families_pos <- match(df_to_use[,2], unique(df_to_use[,2]))
        unique_these_families_pos <- unique(these_families_pos)

        family_transcription_sums_vec <- vector()
        for(j in 1:length(unique_these_families_pos)){
            #j <- 1
            this_unique_these_families_pos <- unique_these_families_pos[j]    
            family_transcription_sums_vec[j] <- sum(df_to_use[which(these_families_pos %in% this_unique_these_families_pos),6])
        }
        maxHeightThisRow <- max(family_transcription_sums_vec)

        # Subset the fake df for this row
        these_fake_pos <- match(df_to_use$family, unique(df_to_use_all$family))
        this_prDF9_fake_col2_for_printing <- prDF9_fake_col2_for_printing[these_fake_pos,]

        # Add to df family averages and remaining average expression as progress down family
        family_expression_average_vec <- vector()
        family_expression_average_vec <- as.numeric(rep("0", nrow(df_to_use)))
        family_expression_average_remaining_vec <- vector()
        family_expression_average_remaining_vec <- as.numeric(rep("0", nrow(df_to_use)))
        family_pos <- match(df_to_use[,2], df_to_use[,2])
        unique_family_pos <- unique(family_pos)
        for(j in 1:length(unique_family_pos)){
            #j <- 1
            this_family_pos <- unique_family_pos[j]
            this_family_all_pos <- which(family_pos %in% this_family_pos)
            this_family_average_expression <- sum(df_to_use[this_family_all_pos, 6])
            family_expression_average_vec[this_family_all_pos] <- this_family_average_expression
            for(k in 1:length(this_family_all_pos)){
                this_family_remaining_average_expression <- sum(df_to_use[this_family_all_pos[k:length(this_family_all_pos)], 6])
                family_expression_average_remaining_vec[this_family_all_pos[k]] <- this_family_remaining_average_expression
            }
        }
        df_to_use_with_family_averages <- cbind(df_to_use, family_expression_average_vec)
        df_to_use <- df_to_use_with_family_averages
        colnames(df_to_use)[[ncol(df_to_use)]] <- "family_expression_average"
        df_to_use_with_family_remaining_averages <- cbind(df_to_use, family_expression_average_remaining_vec)
        df_to_use <- df_to_use_with_family_remaining_averages
        colnames(df_to_use)[[ncol(df_to_use)]] <- "family_expression_remaining_average_for_stacked_print"
        #df_to_use_with_family_averages <- NULL; df_to_use_with_family_remaining_averages <- NULL

        ##
        
        #library("RColorBrewer")
        #require(scales)

        tif_name <- paste0(pCSC_group_name, "_families_expression_bargraph-stacked-", length(HUGO_column_groups_names), ".tif")
        #outFile <- paste0(outDir, "preferential_expression/", tif_name)
        outFile <- paste0(outDir, tif_name)
        x_label <- "Family"
        y_label <- frag_count_method
        plot_title <- ""
        if(i == 1){
            plot_title <- paste0(toupper(tissue_type), ": [", family_type_name, "] Family expression for ", pCSC_group_name, " (", length(pCSC_group_positions), " cells)", " from samples: ", patient_source_name)
        }
        axis_text_size <- 6
    
        number_of_all_families   <- length(unique(df_to_use_all$family))
        getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel2")) # GOOD
        palette_cols_all <- getPalette(number_of_all_families) # palette for all rows
        palette_vec <- vector()
        palette_vec <- palette_cols_all[match(df_to_use$family, unique(df_to_use_all$family))] # palette for family transcripts on this row

        # Create and truncate family transcript labels (set to "") if over max_number_of_transcript_labels_to_show
        transcript_labels_vec <- paste0(df_to_use[,3], "(", round((df_to_use[,6]/df_to_use[,8]), 5) * 100, "%)") # TO USE %
        df_to_use <- cbind(df_to_use, transcript_labels_vec, stringsAsFactors = FALSE)
        colnames(df_to_use)[ncol(df_to_use)] <- "transcript_info_to_show_truncated"
        
        unique_family_pos <- match(df_to_use[, 2], unique(df_to_use[, 2]))
        family_sizes <- rev(sort(table(unique_family_pos)))

        # Decreasing through values until more than max number of transcripts - set to "" after that, and also for previous 0s if they were found
        for(l in 1:length(family_sizes)){
            if(family_sizes[[l]] > max_number_of_transcript_labels_to_show){
                zero_reached <- 0
                oversized_family_unique_index <-  as.integer(names(family_sizes)[l]) # name is family pos number e.g. 1 from 1 1 1 1 2 2 2 3 3 3 3
                oversized_family_transcript_pos <- which(df_to_use[,2] %in% unique(df_to_use[, 2])[oversized_family_unique_index]) # all the pos of this family e.g. 145 146 147 148
        	oversized_family_transcript_pos_rev <- rev(oversized_family_transcript_pos)
        	for(m in 1:length(oversized_family_transcript_pos_rev)){ 
                    this_transcript_expression_proportion <- df_to_use$expression_average_proportion[(oversized_family_transcript_pos_rev[m])]
                    if((this_transcript_expression_proportion == 0) && (zero_reached == 0)){
                        zero_reached = m
                    }
                    if(m > max_number_of_transcript_labels_to_show){
                        if((m == max_number_of_transcript_labels_to_show + 1) && (zero_reached > 0)){ # set to "" from first zero reached to m
                            #cat(paste0("Reset at ", m, "\n"))
                            df_to_use$transcript_info_to_show_truncated[(oversized_family_transcript_pos_rev[zero_reached:m])] <- ""
                            zero_reached <- 0
                        } else { # or otherwise current position only
                            df_to_use$transcript_info_to_show_truncated[(oversized_family_transcript_pos_rev[m])] <- ""
                        }
                    }                    
                }
            }
        }
        
        ##

        if(((i == 1) && (separate_row_for_high_expression_transcripts == FALSE)) ||
           ((i == 1) && (number_of_prDF8o2_families_high_expression_to_separate > 2)) ||
           ((i == 2) && ((separate_row_for_high_expression_transcripts == TRUE) && (number_of_prDF8o2_families_high_expression_to_separate <= 2))) # Not sure why??
           ){ # include legend box
            p <- ggplot2::ggplot(data=df_to_use, aes(x = family, y = expression_average_proportion, fill = transcript_and_family)) +
            geom_bar(stat="identity", color="black", size=0.1, width=1) +
            theme(
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),        
            panel.grid.minor = element_blank(),        
            panel.background = element_blank(),        
            panel.spacing = unit(0.1, "cm"),
            axis.text.x = element_text(color = "black", angle = 270, hjust = 0, vjust = 0.5, size = axis_text_size),
            axis.text.y = element_text(color = "black"),
            legend.text = element_text(size = 6),  
            legend.box = "horizontal",                                              
            legend.position = c(0.9, 0.85), 
            plot.title = element_text(hjust = 0, size = 10),
            plot.margin = ggplot2::margin(0.1, 0.4, 0.1, 0.1, "cm"))
        } else { # no legend box
            p <- ggplot2::ggplot(data = df_to_use, aes(x = family, y = expression_average_proportion, fill = transcript_and_family)) +
            geom_bar(stat = "identity", color = "black", size = 0.1, width = 1) +
            theme(
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),        
            panel.grid.minor = element_blank(),        
            panel.background = element_blank(),        
            panel.spacing = unit(0.1, "cm"),
            axis.text.x = element_text(color = "black", angle = 270, hjust = 0, vjust = 0.5, size = axis_text_size),
            axis.text.y = element_text(color = "black"),
            legend.position = "none",
            plot.title = element_text(hjust = 0, size = 10),
            plot.margin = ggplot2::margin(0.1, 0.4, 0.1, 0.1, "cm"))
        }
        p <- p +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) + 
        #scale_fill_manual(values = palette_vec, guide = FALSE) +
        scale_fill_manual(values = palette_vec, guide = "none") +

        geom_hline(aes(yintercept = all_transcripts_pCSC_group_expression_mean, color = "pCSC group: all transcripts mean", linetype = "pCSC group: all transcripts mean"), size = 0.5, alpha = 0.7) # black dotted

        p <- p + geom_text_repel(aes(label = transcript_info_to_show_truncated, x = family), position = position_stack(vjust = 0.5), hjust = 0.5, color = "black", size = 2, segment.color = 'black', segment.size = 0.1, force = 0.8, direction = "y")
        # Alternative text formats
        #p <- p + geom_text(aes(label = transcript_info_to_show_truncated, x = family), color = "black", position = position_stack(vjust = 0.5), size = 2, hjust = 0.5)
        #geom_text(aes(label = transcript, x = family), color = "black", position = position_stack(vjust = 0.5), size = 2, hjust = 0.5)
        #geom_text(aes(label = cell_numbers_to_show, x=transcript_name_to_show), color = "black", position = position_stack(vjust = 0.5), size = 2, hjust = 0.5)
        #geom_text(aes(label = transcript_info_to_show, x=family, y = 2), color = "black", position = position_nudge(x = -0.4), size = 1.4, vjust = 0, hjust = 0)
        #geom_text_repel(aes(label = transcript_info_to_show_truncated, x = family, y = (family_expression_remaining_average_for_stacked_print-expression_average_proportion)),
        #                    color = "black", size = 2, segment.color = 'black', segment.size = 0.1, force = 1)
        #geom_text(aes(label = transcript_info_to_show, x=family, y = 2), color = "black", position = position_nudge(x = -0.4), size = 1.4, vjust = 0, hjust = 0)

        #if(length(this_patient_tumour_only_positions) > 0){ # if this_patient_tumour_only_positions are not the same as tumour_only_positions
        #if(length(which(this_patient_tumour_only_positions != tumour_only_positions)) > 0){ # if this_patient_tumour_only_positions are not the same as tumour_only_positions
        #if(!identical(sort(this_patient_tumour_only_positions), sort(tumour_only_positions))){ # if this_patient_tumour_only_positions are not the same as tumour_only_positions
        if((length(this_patient_tumour_only_positions) > 0) && (!identical(sort(this_patient_tumour_only_positions), sort(tumour_only_positions)))){ # if this_patient_tumour_only_positions are not the same as tumour_only_positions
            p <- p +
                     geom_line(data = this_prDF9_fake_col2_for_printing, aes(x = family, y = this_tumour_only_expression, color = "pCSC group source patient all tumour cells average",
                                                                             linetype = "pCSC group source patient all tumour cells average"), group = 1, size = 0.5, alpha = 0.5) + # darrkgreen dashed

                     geom_line(data = this_prDF9_fake_col2_for_printing, aes(x = family, y = tumour_only_expression, color = "Tumour cells only (all patients) average", linetype = "Tumour cells only (all patients) average"), group = 1, size = 0.5, alpha = 0.5) + # red dashed
                     geom_line(data = this_prDF9_fake_col2_for_printing, aes(x = family, y = all_samples_expression, color = "All cells (all patients) average", linetype = "All cells (all patients) average"), group = 1, size = 0.5, alpha = 0.5) + # blue dashed
                     geom_line(data = this_prDF9_fake_col2_for_printing, aes(x = family, y = pCSC_group_tumour_difference_expression, color = "pCSC group cells average expression minus\ntumour cells (all patients) average expression",
                                                                             linetype = "pCSC group cells average expression minus\ntumour cells (all patients) average expression"), group = 1, size = 0.5, alpha = 0.7) + # deeppink1 

                     scale_color_manual(breaks = c("pCSC group cells average expression minus\ntumour cells (all patients) average expression", "pCSC group source patient all tumour cells average", "Tumour cells only (all patients) average",
                                                   "All cells (all patients) average", "pCSC group: all transcripts mean"),
                                        values = c("pCSC group cells average expression minus\ntumour cells (all patients) average expression"="deeppink1", "pCSC group source patient all tumour cells average"="darkgreen",
                                                   "Tumour cells only (all patients) average"="red", "All cells (all patients) average"="blue", "pCSC group: all transcripts mean"="black")) +
                     scale_linetype_manual(breaks = c("pCSC group cells average expression minus\ntumour cells (all patients) average expression", "pCSC group source patient all tumour cells average", "Tumour cells only (all patients) average",
                                                      "All cells (all patients) average", "pCSC group: all transcripts mean"),
                                           values = c("pCSC group cells average expression minus\ntumour cells (all patients) average expression"="twodash", "pCSC group source patient all tumour cells average"="dashed",
                                                      "Tumour cells only (all patients) average"="dashed", "All cells (all patients) average"="dashed", "pCSC group: all transcripts mean"="dotted")) + 
                     labs(title = plot_title, x = x_label, y = y_label, colour = "Cell Group Expression Averages", linetype = "Cell Group Expression Averages", size = "", fill = "", shape = "")
        } else {
            p <- p +
                     geom_line(data = this_prDF9_fake_col2_for_printing, aes(x = family, y = tumour_only_expression, color = "Tumour cells only average", linetype = "Tumour cells only average"), group = 1, size = 0.5, alpha = 0.5) + # red dashed
                     geom_line(data = this_prDF9_fake_col2_for_printing, aes(x = family, y = all_samples_expression, color = "All cells average", linetype = "All cells average"), group = 1, size = 0.5, alpha = 0.5) + # blue dashed
                     geom_line(data = this_prDF9_fake_col2_for_printing, aes(x = family, y = pCSC_group_tumour_difference_expression, color = "pCSC group cells average expression minus\ntumour cells average expression",
                                                                             linetype = "pCSC group cells average expression minus\ntumour cells average expression"), group = 1, size = 0.5, alpha = 0.7) + # deeppink1 

                     scale_color_manual(breaks = c("pCSC group cells average expression minus\ntumour cells average expression", "Tumour cells only average", "All cells average", "pCSC group: all transcripts mean"),
                                        values = c("pCSC group cells average expression minus\ntumour cells average expression"="deeppink1", "Tumour cells only average"="red",
                                                   "All cells average"="blue", "pCSC group: all transcripts mean"="black")) +
                     scale_linetype_manual(breaks = c("pCSC group cells average expression minus\ntumour cells average expression", "Tumour cells only average", "All cells average", "pCSC group: all transcripts mean"),
                                           values = c("pCSC group cells average expression minus\ntumour cells average expression"="twodash", "Tumour cells only average"="dashed",
                                                      "All cells average"="dashed", "pCSC group: all transcripts mean"="dotted")) + 
                     labs(title = plot_title, x = x_label, y = y_label, colour = "Cell Group Expression Averages", linetype = "Cell Group Expression Averages", size = "", fill = "", shape = "")
        }
        p_list[[i]] <- p
        #ggplot2::ggsave(p, file = outFile, device = "tiff", compression = "lzw", dpi = 150, width = 49, height = 9, units = "in") # For printing each row i.e. dev testing        
    }

    #library("cowplot")

    graph_row_lengths <- list()
    for(i in 1:length(p_list)){
	graph_row_lengths[[i]] <- length(unique(p_list[[i]][[1]][[2]]))
    }

    p_draw <- cowplot::ggdraw()
    p_list_len <- length(p_list)
    max_row_len <- max(unlist(graph_row_lengths))
    
    for(i in 1:length(p_list)){
        if(graph_row_lengths[[i]] < max_row_len){
            p_draw <- p_draw + cowplot::draw_plot(p_list[[i]], x = 0, y = (1 - ((1/p_list_len) * i)), width = (((1/max_row_len) * (graph_row_lengths[[i]])) + 0.05), height = (1/p_list_len)) # not perfect but OK for the moment.
        } else {
            p_draw <- p_draw + cowplot::draw_plot(p_list[[i]], x = 0, y = (1 - ((1/p_list_len) * i)), width = ((1/max_row_len) * (graph_row_lengths[[i]])), height = (1/p_list_len))
        }
    }

    if(graph_length_limit <= 40){
        width_scale <- 1.3
    } else {
        width_scale <- 1
    }

    dpiVal <- 150
    thisTifWidth <- (max(unlist(graph_row_lengths))/width_scale)  + 2
    if(thisTifWidth < 4){  thisTifWidth <- 4  } 
    if(thisTifWidth > 49){  thisTifWidth <- 49  }
    thisTifHeight <- (7 * number_of_graph_rows)

    if(thisTifHeight > 49){  thisTifHeight <- 49  }

    #outFile <- paste0(outDir, "preferential_expression/", tif_name)
    outFile <- paste0(outDir, tif_name)
    if(thisTifHeight > 120){
        thisTifHeight <- 120
        dpiVal <- 100
    }
    #cat(paste0("tif height = ", thisTifHeight, "\n"))

    tiff(filename = outFile, compression = "lzw", res = dpiVal, width = thisTifWidth, height = thisTifHeight, units = "in", pointsize = 12)
    print(p_draw)
    dev.off()    
    #ggplot2::ggsave(p_draw, file = outFile, device = "tiff", compression = "lzw", dpi = dpiVal, width = thisTifWidth, height = thisTifHeight, units = "in")
    #ggplot2::ggsave(p_draw, file = outFile, device = "jpg", compression = "lzw", dpi = dpiVal, width = thisTifWidth, height = thisTifHeight, units = "in")
    print(p_draw)

    family_ranks_to_use_list <- list()
    family_ranks_to_use_list[[1]] <- ranks2_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__FAMILIES
    family_ranks_to_use_list[[2]] <- ranks3_STACKED_FAMILY_EXPRESSION_BARCHART_SINGLE_BARS__TRANSCRIPTS_PERCENTEXPRESSION_rev_sorted_no_duplicates_no_zeros
    
    return(family_ranks_to_use_list)
}

