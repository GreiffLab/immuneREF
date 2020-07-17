#' Condense single IMMUNENET LAYERS
#'
#' @param list_single_layers List of calculated immuneREF similarity layers
#' @param weights Weights for each layer
#' @param method Method for condensing (so far only "standard" is available)
#' @return output matrix
#' @examples
#' \dontrun{
#' condense_layers(list_single_layers, weights, method="spearman")
#' }

condense_layers <- function(list_single_layers,weights,method="standard"){
  #condenses single layers depending on weights

  if(method=="standard"){
    wt <- weights

    wt_norm<-wt/sum(wt)
    #  Make a 3D array from list of matrices
    arr<-array(unlist(list_single_layers),c(nrow(list_single_layers[[1]]),nrow(list_single_layers[[1]]),length(list_single_layers)))

    #  Get mean of third dimension
    output_matrix<-apply(arr,1:2, weighted.mean,wt_norm)
    rownames(output_matrix)<-rownames(list_single_layers[[1]])
    colnames(output_matrix)<-rownames(output_matrix)
  }else if(method=="majority_vote"){
    cat("majority vote not implemented yet \n")
    output_matrix <- NA
  }else{
    cat("not implemented yet \n")
    output_matrix <- NA
  }

  return(output_matrix)

}
