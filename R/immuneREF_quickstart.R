#' Quickstart function to calculate immuneREF layers from list of reperotires
#'
#' @param repertoire_list List of AIRR-compliant repertoires.
#' @param convergence_layer Determines whether Convergence layer is calculated using Overlap or Immunosignatures
#' @return List of immuneREF similarity layers
#' @examples
#' \dontrun{
#' immuneREF_quickstart(tutorial_repertoires)
#' }


immuneREF_quickstart <- function(repertoire_list,convergence_layer="overlap"){

  if(convergence_layer=="overlap"){
    #calculate overlap layer as a backup layer (in place of immunosignature layer)
    overlap_layer<-repertoire_overlap(repertoire_list,basis="CDR3_aa")
  }else{
    overlap_layer <- ""
  }
  #for naming: size of repertoires (as they are simulated to 1000 no need to subsample)
  subsample_size<-10000

  #for naming: save length of repertoire and names
  repertoire_lengths <- sapply(1:length(repertoire_list),function(x) nrow(repertoire_list[[x]]))
  repertoire_names <- sapply(1:length(repertoire_list),function(x) as.character(unique(repertoire_list[[x]]$name_repertoire)))


  #calculate all features for each repertoire using the calc_characteristics function
  repertoires_analyzed<-list()
  for(i in 1:length(repertoire_list)){
    repertoires_analyzed[[repertoire_names[i]]]<-calc_characteristics(
      repertoire_df=repertoire_list[[i]],
      species=strsplit(repertoire_names[i],"_")[[1]][2],
      receptor=strsplit(repertoire_names[i],"_")[[1]][3],
      chain=strsplit(repertoire_names[i],"_")[[1]][4],
      identifier_rep=repertoire_names[i])
  }

  #calculate similarities across layers
  cat("Calculate Repertoire similarities \n")

  list_single_layers<-calculate_similarities(
    repertoires_analyzed=repertoires_analyzed,
    overlap=overlap_layer)

  cat("Calculate Condensed network \n")
  #Calculate condensed network (here equal weights for each layer)
  cormat <- condense_layers(list_single_layers,
    weights = c(1,1,1,1,1,1),
    method = "standard")

  list_all_layers <- list_single_layers
  list_all_layers[["Condensed"]] <- cormat

  return(list_all_layers)


}
