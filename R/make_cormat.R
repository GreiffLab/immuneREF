#' Calculates similarities between processed repertoires
#'
#' @param repertoires_analyzed Output of feature analysis
#' @param weights_VDJ V,D,J,VJ weights for similarity calculation
#' @param weights_overall Layer weights for similarity calculation
#' @param correlation_method Method for correlation calculations
#' @param arch_method How architecture layer should be calculated (Default="all")
#' @return Matrix containing similarity of each repertoire pair
#' @examples
#' \dontrun{
#' make_cormat(repertoires_analyzed,weights_VDJ=c(1,1,1,0),
#' weights_overall=c(1,1,1,1,1,1),correlation_method='pearson',
#' arch_method="all")
#' }


make_cormat<-function(repertoires_analyzed,weights_VDJ=c(1,0,1,0),weights_overall=c(1,1,1,1,1,1),correlation_method='pearson',arch_method="all"){
  #calculate all correlations in repertoires_analyzed list
  #and output correlation matrix for a given set of feature weights and VDJ weights.
  mat_combinations<-t(combn(c(1:length(repertoires_analyzed)),2))
  #mat_combinations<-expand.grid(c(1:length(repertoires_analyzed)),c(1:length(repertoires_analyzed)))

  names_samples<-sapply(1:length(repertoires_analyzed),function(x) repertoires_analyzed[[x]][[9]])

  cor_matrix<-diag(1,length(names_samples),length(names_samples))
  dimnames(cor_matrix)<-list(names_samples,names_samples)

  #fill in correlation matrix with correlations between samples
  for(i in 1:nrow(mat_combinations)){
    sample_1<-names_samples[mat_combinations[i,1]]
    sample_2<-names_samples[mat_combinations[i,2]]

    cor_matrix[sample_1,sample_2]<-.comparing_nodes(rep_1=repertoires_analyzed[[sample_1]],rep_2=repertoires_analyzed[[sample_2]],weights_VDJ=weights_VDJ,weights_overall=weights_overall,cor_method=correlation_method,arch_method=arch_method)
    cor_matrix[sample_2,sample_1]<-cor_matrix[sample_1,sample_2]
  }

  return(cor_matrix)

}
