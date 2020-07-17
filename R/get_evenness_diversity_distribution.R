# Calculates Evenness and Diversity measure of repertoire
# params repertoire_df as obtained from load_mixcr_output
# reqs vegan package for renyi function)
# returns list of evenness and diversity profiles of repertoire (in that order)
#



#clonal expansion
.get_evenness_diversity_distribution<-function(repertoire_df){
  #calculates evenness and diversity of input repertoire
  #returns list. first entry evenness, second entry: diversity

  #steps for div parameter
  alphas <- seq(0,10,0.1)

  #evenness,diversity of unique nt seqs (DOUBLE CHECK DIVERSITY)
  profiles_evenness_matrix  <- sapply(1:length(alphas), function(x) exp(vegan::renyi(repertoire_df[, "freqs"]/sum(repertoire_df[, "freqs"]), scales = alphas[x])))/nrow(repertoire_df)
  profiles_diversity_matrix  <- sapply(1:length(alphas), function(x) vegan::renyi(repertoire_df[, "freqs"]/sum(repertoire_df[, "freqs"]), scales = alphas[x]))

  list_results<-list()
  list_results[[1]]<-profiles_evenness_matrix
  list_results[[2]]<-profiles_diversity_matrix
  names(list_results)<-c("evenness","diversity")

  list_results
}
