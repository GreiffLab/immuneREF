.adj_matrix_calc<-function(sequence_vec){
  #create levenshtein distance matrix for you sequences for downstream graph analysis
  ld_matrix <- stringdist::stringdistmatrix(sequence_vec,sequence_vec,method="lv")
  rownames(ld_matrix)<-sequence_vec
  colnames(ld_matrix)<-sequence_vec

  #turn ld_matrix into adjacency matrix
  adj_matrix<-ld_matrix

  adj_matrix

}
