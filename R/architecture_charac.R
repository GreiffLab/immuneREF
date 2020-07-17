
.architecture_charac<-function(repertoire_df,max_ld,option="AA"){
  #creates a variety of measures related to repertoire architecture based on Miho2017.
  #builds graph from levenshtein distance matrix and calculates key indicators using igraph package

  if(option=="AA"){
    CDR3s_char<-as.character(repertoire_df$junction_aa)
  } else if(option=="nt"){
    CDR3s_char<-toupper(as.character(repertoire_df$junction))
  }

  #check if more than 10000 if yes: subsample to 10000
  if(length(CDR3s_char)>10000){
    CDR3s <- sample(CDR3s_char,10000)
  }else{
    CDR3s <- CDR3s_char
  }

  #take unique CDR3s (faster network calculation, no self reference)
  #Important: take unique after subsampling to 10'000 otherwise there will be a bias against shorter seqs
  CDR3s <- unique(CDR3s)

  #create graph object from adjaceny matrix (similarity score matrix)
  adj_matrix<-.adj_matrix_calc(CDR3s)

  #for each levenshtein distance from 1:max_ld calcualte network features
  list_results_graph_per_ld <- list()
  for(i in 1:max_ld){
    #cat("architecture LD = ",i,"\n")
    adj_matrix_ld <- adj_matrix #make new var to change

    similarity <- i

    #set all distances that are smaller or equal to threshold to 1
    #then 0 out all that are not 1 (all above similarity threshold)
    adj_matrix_ld[adj_matrix_ld <= similarity]<-1
    adj_matrix_ld[adj_matrix_ld!=1]<-0

    #create graph object from adj matrix of current edit distance
    CDR3s_graph<-igraph::graph_from_adjacency_matrix(adj_matrix_ld,mode=c('undirected'),diag = FALSE)

    #evaluate components of graph
    components_of_graph <- igraph::components(CDR3s_graph)

    #calculate 9 networks features using the igraph package
    list_results_graph<-list()
    list_results_graph[[1]]<-c(Nodes=igraph::gorder(CDR3s_graph),Edges=igraph::ecount(CDR3s_graph)) #number of nodes and edges
    list_results_graph[[2]]<-2*igraph::ecount(CDR3s_graph)/igraph::gorder(CDR3s_graph) #avg degree
    list_results_graph[[3]]<-igraph::mean_distance(CDR3s_graph, directed=F) #mean distance
    list_results_graph[[4]]<-igraph::diameter(CDR3s_graph, directed=F) #diameter of graph
    list_results_graph[[5]]<-mean(igraph::hub_score(CDR3s_graph, weights=NA)$vector,na.rm=TRUE) #mean hub scores across all nodes
    list_results_graph[[6]]<-igraph::transitivity(CDR3s_graph,type=c("global")) #global clustering coeff
    list_results_graph[[7]]<-components_of_graph$no #total number of clusters
    list_results_graph[[8]]<-max(components_of_graph$csize)/length(CDR3s) #percent of sequences that are in main cluster
    list_results_graph[[9]]<-table(components_of_graph$csize) #components/cluster and number of members
    list_results_graph[[10]]<-igraph::degree_distribution(CDR3s_graph) #components/cluster and number of members
    list_results_graph[[11]]<-igraph::degree_distribution(CDR3s_graph,cumulative=TRUE) #components/cluster and number of members

    names(list_results_graph)<-c("Nodes_Edges","avg_degree",'mean_dist','diameter',"mean_hub_score","transitivity_global",'nb_components',"main_cluster_perc","component_size_and_members","degree_distribution","degree_distribution_cumulative")

    list_results_graph_per_ld[[i]]<-list_results_graph
  }

  list_results_graph_per_ld

}
