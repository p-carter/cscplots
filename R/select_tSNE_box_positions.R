#' Select cells from a t-SNE that are located within a defined rectangle.
#'
#' Identify and return the positions of cells that are located in the t-SNE within
#'     the rectangle defined by parameters Y1_borders and Y2_borders. Can be used,
#'     for example, to find the cells that make up a cluster of cells in the t-SNE.
#'
#' @param tSNE_info Data frame where first 2 columns are Y1 and Y2 results from a
#'     t-SNE generated using a cancer poulation expression matrix. Rows are cancer
#'     cells. Additional columns may contain additional info, such as scores of
#'     the cells or clinical info about the cells.
#' @param Y1_borders Opposite sides of the rectangle to locate cells within it, in
#'     Y1 positions in a vector, e.g. Y1_borders = c(3,5).
#' @param Y2_borders Opposite sides of the rectangle in Y2 positions in a
#'     vector, e.g. Y2_borders = c(2,4).
#' @return Positions of the cells in the DF rows, if any, that are located in the
#'     rectangle defined by Y1_borders and Y2_borders.
#' @export
select_tSNE_box_positions <- function(tSNE_info, Y1_borders, Y2_borders){

    if((length(Y1_borders) == 2) && (length(Y2_borders) == 2)){
        Y1_border1 = Y1_borders[1]
        Y1_border2 = Y1_borders[2]
        Y2_border1 = Y2_borders[1]
        Y2_border2 = Y2_borders[2]
        if(Y1_border1 > Y1_border2){
            tmpVal <- Y1_border1
            Y1_border1 <- Y1_border2
            Y1_border2 <- tmpVal
        }
        if(Y2_border1 > Y2_border2){
            tmpVal <- Y2_border1
            Y2_border1 <- Y2_border2
            Y2_border2 <- tmpVal
        }

        box_pos1 <- which(tSNE_info$Y1 > Y1_border1)
        box_pos2 <- which(tSNE_info$Y1 < Y1_border2)
        box_pos3 <- which(tSNE_info$Y2 > Y2_border1)
        box_pos4 <- which(tSNE_info$Y2 < Y2_border2)
        
        box_subset1 <- box_pos1[which(box_pos1 %in% box_pos2)]
        box_subset2 <- box_pos3[which(box_pos3 %in% box_pos4)]
        
        box_contents_pos <- box_subset1[which(box_subset1 %in% box_subset2)]
        
        return(box_contents_pos)
    }

    return(NULL)
}
