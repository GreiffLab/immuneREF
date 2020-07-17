
.graph_similarity<-function(sequence_vec,similarity=1){
  #create levenshtein distance matrix for you sequences for downstream graph analysis

  ld_matrix <- stringdist::stringdistmatrix(sequence_vec,sequence_vec)
  rownames(ld_matrix)<-sequence_vec
  colnames(ld_matrix)<-sequence_vec

  #turn ld_matrix into adjacency matrix
  adj_matrix<-ld_matrix
  adj_matrix[adj_matrix > similarity]<-0
  adj_matrix[adj_matrix==similarity]<-1

  #construct graph
  graph_similarity<-igraph::graph_from_adjacency_matrix(adj_matrix,mode=c('undirected'))

  graph_similarity

}
