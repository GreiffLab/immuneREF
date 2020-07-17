#' Calculates Similarity Layers BETWEEN REPERTOIRES
#'
#' @param repertoires_analyzed immuneREF analyzed repertoires
#' @param overlap_layer overlap_layer if calculated_previously and convergence method is 'overlap'
#' @param vdj_vj_weights How VDJ usage should be treated for similarity calculation (V,D,J,VJusage on/off)
#' @param convergence Which convergence measure should be used overlap or immunosignatures (Default: overlap)
#' @param cor_methods vector of lengths 6 determining which correlation method should be used for each layer.
#' @return list_single_layers
#' @examples
#' \dontrun{
#' calculate_similarities(repertoires_analyzed,vdj_vj_weights=c(1,0,1,0),convergence="overlap",overlap_layer="",cor_methods=c("Diversity"="pearson","AAfreq"='pearson',"Architecture"="pearson","Convergence"='pearson',"VDJ_usage"="pearson","k-mers"='pearson'))
#' }

calculate_similarities <- function(repertoires_analyzed,overlap_layer="",vdj_vj_weights=c(1,1,1,0),convergence="overlap", cor_methods=c("Diversity"="pearson","AAfreq"='pearson',"Architecture"="pearson","Convergence"='pearson',"VDJ_usage"="pearson","k-mers"='pearson')){
  list_single_layers<-list()
  list_single_layers[["Diversity"]]<-make_cormat(repertoires_analyzed,weights_overall=c(1,0,0,0,0,0),correlation_method=cor_methods["Diversity"])
  list_single_layers[["AAfreq"]]<-make_cormat(repertoires_analyzed,weights_overall=c(0,1,0,0,0,0),correlation_method=cor_methods["AAfreq"])
  list_single_layers[["Architecture"]]<-make_cormat(repertoires_analyzed,weights_overall=c(0,0,1,0,0,0),correlation_method=cor_methods["Architecture"])
  if(convergence=="overlap" & class(overlap_layer)!="character"){
    list_single_layers[["Convergence"]]<-overlap_layer
  }else{
    list_single_layers[["Convergence"]]<-make_cormat(repertoires_analyzed,weights_overall=c(0,0,0,1,0,0),correlation_method=cor_methods["Convergence"])
  }
  list_single_layers[["VDJ_usage"]]<-make_cormat(repertoires_analyzed,weights_VDJ=vdj_vj_weights,weights_overall=c(0,0,0,0,1,0),correlation_method=cor_methods["VDJ_usage"])
  list_single_layers[["k-mers"]]<-make_cormat(repertoires_analyzed,weights_overall=c(0,0,0,0,0,1),correlation_method=cor_methods["k-mers"])

  return(list_single_layers)
}
