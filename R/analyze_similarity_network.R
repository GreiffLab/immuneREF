#' Analyzes immuneNET similarity network
#'
#' @param similarity_network immuneREF output layer
#' @return List of network features
#' @examples
#' \dontrun{
#' analyze_similarity_network(similarity_network)
#' }

#simple analysis of network features

analyze_similarity_network <- function(similarity_network){
  #Calculate Network Features

  #create graph object from adjaceny matrix (similarity score matrix)
  similarity_graph <- igraph::graph_from_adjacency_matrix(
    similarity_network,
    mode = c("undirected"),
    weighted = TRUE)

  #evaluate components of graph
  components_of_graph <- igraph::components(similarity_graph)

  #calculate 12 networks features using the igraph package
  network_features <- list()
  network_features[[1]]<-igraph::degree(similarity_graph) #degrees of nodes
  network_features[[2]]<-2*igraph::ecount(similarity_graph)/igraph::gorder(similarity_graph) #avg degree
  network_features[[3]]<-igraph::degree.distribution(similarity_graph, cumulative=T, mode="all") #degree distribution
  network_features[[4]]<-igraph::mean_distance(similarity_graph, directed=F) #mean distance
  network_features[[5]]<-igraph::diameter(similarity_graph, directed=F) #diameter of graph
  network_features[[6]]<-mean(igraph::hub_score(similarity_graph, weights=NA)$vector,na.rm=TRUE) #mean hub scores across all nodes
  network_features[[7]]<-igraph::betweenness(similarity_graph, directed=F, weights=NA) #betweenness of nodes
  network_features[[8]]<-igraph::transitivity(similarity_graph,type=c("global")) #global clustering coeff
  network_features[[9]]<-igraph::transitivity(similarity_graph,type=c("local")) #local clustering coeff
  network_features[[10]]<-components_of_graph$no #total number of clusters
  network_features[[11]]<-base::max(components_of_graph$csize)/base::nrow(similarity_network) #percent of sequences that are in main cluster
  network_features[[12]]<-c(Nodes=igraph::gorder(similarity_graph),Edges=igraph::ecount(similarity_graph)) #number of nodes and edges

  #name each list entry
  names(network_features)<-c("degrees",
                             "avg_degree",
                             "degree_distribution",
                             "mean_dist",
                             "diameter",
                             "hub_score",
                             "betweenness",
                             "transitivity_global",
                             "transitivity_local",
                             "Nb_clusters",
                             "main_cluster_perc",
                             "Nodes_Edges")

  return(network_features)
}
