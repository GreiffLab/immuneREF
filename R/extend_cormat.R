# Extends existing correlation matrix with new data
#
# @param repertoires_analyzed x
# @param existing_cormat x
# @param weights_VDJ x
# @param weights_overall x
# @param correlation_method x
# @return Correlation matrix
# @examples
# \dontrun{
# extend_cormat(tutorial_repertoires, existing_cormat,
# weights_VDJ, weights_overall, correlation_method)
# }


.extend_cormat<-function(repertoires_analyzed,existing_cormat,weights_VDJ,weights_overall,correlation_method='pearson'){
  #extends an existing correlation matrix with novel repertoires added circumventing calculation of full cor matrix.
  #and outputs new extended correlation matrix for a given set of feature weights and VDJ weights.
  mat_combinations<-t(combn(c(1:length(repertoires_analyzed)),2))
  #mat_combinations<-expand.grid(c(1:length(repertoires_analyzed)),c(1:length(repertoires_analyzed)))

  names_samples<-sapply(1:length(repertoires_analyzed),function(x) repertoires_analyzed[[x]][[9]])

  cor_matrix<-diag(1,length(names_samples),length(names_samples))
  dimnames(cor_matrix)<-list(names_samples,names_samples)

  existing_samples <- rownames(existing_cormat)
  new_samples <- names_samples[!(names_samples %in% rownames(existing_cormat))]

  #fill in correlation matrix with correlations between samples
  for(i in 1:nrow(mat_combinations)){
    sample_1<-names_samples[mat_combinations[i,1]]
    sample_2<-names_samples[mat_combinations[i,2]]
    if((sample_1 %in% new_samples) | (sample_2 %in% new_samples)){
      cor_matrix[sample_1,sample_2]<-.comparing_nodes(rep_1=repertoires_analyzed[[sample_1]],rep_2=repertoires_analyzed[[sample_2]],weights_VDJ=weights_VDJ,weights_overall=weights_overall,cor_method=correlation_method)
    } else{
      cor_matrix[sample_1,sample_2] <- existing_cormat[sample_1,sample_2]
    }
    cor_matrix[sample_2,sample_1] <- cor_matrix[sample_1,sample_2]

  }

  return(cor_matrix)

}
