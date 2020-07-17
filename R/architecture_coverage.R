
.architecture_coverage<-function(repertoire_df,max_ld=4,option="AA"){
  #creates a variety of measures related to repertoire architecture based on Miho2017.
  #builds graph from levenshtein distance matrix and calculates key indicators using igraph package

  #set variable CDR3s depending on option (either nucleotide or amino acid)
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

  #call adj_matrix calculation
  adj_matrix<-.adj_matrix_calc(CDR3s)

  #for each levenshtein distance in 1:max_ld calculate network and components
  no_of_clusters_list <- list()
  perc_in_main_list <- list()
  for(i in 1:max_ld){
    cat("architecture LD = ",i,"\n")
    adj_matrix_ld <- adj_matrix #make new var to change

    similarity <- i

    adj_matrix_ld[adj_matrix_ld<=similarity]<-1
    adj_matrix_ld[adj_matrix_ld>similarity]<-0

    #create graph object from adj matrix
    CDR3s_graph<-igraph::graph_from_adjacency_matrix(adj_matrix_ld,mode=c('undirected'))

    #find components and evaluate number and size(percentage of seqs) in each component
    components_of_graph <- igraph::components(CDR3s_graph)
    no_of_clusters_list[[i]] <- components_of_graph$no #total number of clusters
    perc_in_main_list[[i]] <- max(components_of_graph$csize)/length(CDR3s) #percent of sequences that are in main cluster
  }

  #save in list --> perc in main cluster across edit distances and number of components at each edit distance
  arch_coverage <- list(perc_in_main=unlist(perc_in_main_list), no_of_clusters=unlist(no_of_clusters_list))

  return(arch_coverage)

}
