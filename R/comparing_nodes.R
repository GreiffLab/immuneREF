

.comparing_nodes<-function(rep_1,rep_2,weights_VDJ,weights_overall,cor_method,arch_method){
  #compares each repertoire characteristic btw 2 repertoires
  #based on weighted characteristics and weights for relative V,D,J importance
  #correlation is calculated based on method input. default: pearson.

  #comparing two nodes:
  node_correlations<-list()
  #evenness #pearson otherwise 1
  if(weights_overall[1]!=0){
    node_correlations[[1]]<-abs(cor(rep_1[[1]][[1]],rep_2[[1]][[1]],method=cor_method))
  }else{
    node_correlations[[1]]<-0
  }
  #AA freq (fixed for all lengths) ks: ks<-ks.test(A,B)$statistic
  if(weights_overall[2]!=0){
    node_correlations[[2]]<-mean(sapply(1:length(rep_1[[2]]), function(x)
      mean(sapply(1:length(rep_1[[2]][[x]]),function(y)
        abs(cor(as.vector(rep_1[[2]][[x]][[y]]),as.vector(rep_2[[2]][[x]][[y]]),method=cor_method,use="everything"))),na.rm=TRUE)),na.rm=TRUE)
  }else{
    node_correlations[[2]]<-0
  }
  #rep architecture #for now: only including avg degree.
  if(weights_overall[3]!=0){
    if(arch_method=="all"){
      list_arch<-list()
      shortest <- min(length(rep_1[[3]][[1]][["degree_distribution_cumulative"]]), length(rep_2[[3]][[1]][["degree_distribution_cumulative"]]))
      deg_dist_1 <- tail(rep_1[[3]][[1]][["degree_distribution_cumulative"]], shortest)
      deg_dist_2 <- tail(rep_2[[3]][[1]][["degree_distribution_cumulative"]], shortest)
      list_arch[[1]]<-cor(deg_dist_1,deg_dist_2,method="pearson") #correlate degree distributions
      #list_arch[[1]]<-1-(abs(rep_1[[3]][[1]][["avg_degree"]]-rep_2[[3]][[1]][["avg_degree"]])) #difference btw avg degree
      list_arch[[2]]<-1-(abs(rep_1[[3]][[1]][["mean_hub_score"]]-rep_2[[3]][[1]][["mean_hub_score"]])) #difference btw mean hub score
      #list_arch[[3]]<-1-(abs(rep_1[[3]][[1]][["transitivity_global"]]-rep_2[[3]][[1]][["transitivity_global"]])) #difference btw global_clust
      list_arch[[3]]<-1-(abs(rep_1[[3]][[1]][["nb_components"]]-(rep_2[[3]][[1]][["nb_components"]]))/max(rep_1[[3]][[1]][["Nodes_Edges"]][["Nodes"]],rep_2[[3]][[1]][["Nodes_Edges"]][["Nodes"]])) #difference btw nb of unconnected clusters
      list_arch[[4]]<-1-(abs(rep_1[[3]][[1]][["main_cluster_perc"]]-rep_2[[3]][[1]][["main_cluster_perc"]])) #difference btw perc of seqs in main cluster
      list_arch_vec<-unlist(list_arch)
      if(TRUE %in% (list_arch_vec>1)){
        cat("Warning: Repertoire architecture contains invalid entry")
      }
      node_correlations[[3]]<-mean(list_arch_vec,na.rm=TRUE)
    }else if(arch_method=="coverage"){
      #extract vector of perc in max cluster
      rep_1_percmax <- unlist(sapply(1:length(rep_1[[3]]),function(x) rep_1[[3]][[x]][["main_cluster_perc"]]))
      rep_2_percmax <- unlist(sapply(1:length(rep_2[[3]]),function(x) rep_2[[3]][[x]][["main_cluster_perc"]]))

      node_correlations[[3]] <- abs(cor(rep_1_percmax,rep_2_percmax,method='pearson'))
    }else if(arch_method=="mix"){
      #extract vector of perc in max cluster
      rep_1_percmax <- unlist(sapply(1:length(rep_1[[3]]),function(x) rep_1[[3]][[x]][["main_cluster_perc"]]))
      rep_2_percmax <- unlist(sapply(1:length(rep_2[[3]]),function(x) rep_2[[3]][[x]][["main_cluster_perc"]]))

      list_arch<-list()
      shortest <- min(length(rep_1[[3]][[1]][["degree_distribution_cumulative"]]), length(rep_2[[3]][[1]][["degree_distribution_cumulative"]]))
      deg_dist_1 <- tail(rep_1[[3]][[1]][["degree_distribution_cumulative"]], shortest)
      deg_dist_2 <- tail(rep_2[[3]][[1]][["degree_distribution_cumulative"]], shortest)
      list_arch[[1]]<-cor(deg_dist_1,deg_dist_2,method="pearson") #correlate degree distributions
      #list_arch[[1]]<-1-(abs(rep_1[[3]][[1]][["avg_degree"]]-rep_2[[3]][[1]][["avg_degree"]])) #difference btw avg degree
      list_arch[[2]]<-1-(abs(rep_1[[3]][[1]][["mean_hub_score"]]-rep_2[[3]][[1]][["mean_hub_score"]])) #difference btw mean hub score
      #list_arch[[3]]<-1-(abs(rep_1[[3]][[1]][["transitivity_global"]]-rep_2[[3]][[1]][["transitivity_global"]])) #difference btw global_clust
      list_arch[[3]]<-1-(abs(rep_1[[3]][[1]][["nb_components"]]-(rep_2[[3]][[1]][["nb_components"]]))/max(rep_1[[3]][[1]][["Nodes_Edges"]][["Nodes"]],rep_2[[3]][[1]][["Nodes_Edges"]][["Nodes"]])) #difference btw nb of unconnected clusters
      list_arch[[4]]<-abs(cor(rep_1_percmax,rep_2_percmax,method='pearson'))
      list_arch_vec<-unlist(list_arch)
      if(TRUE %in% (list_arch_vec>1)){
        cat("Warning: Repertoire architecture contains invalid entry")
      }
      node_correlations[[3]]<-mean(list_arch_vec,na.rm=TRUE)
    }
  }else{
    node_correlations[[3]]<-0
  }

  #rep immunosignature_presence #Extende for use with more than one model (not yet done in main functions. 20180112)
  #if one of two entries is 0 -->  0 (no pos classified seqs in one rep: very different reps wrt to the class),
  #otherwise 1-difference between the two (when some pos classified)
  if(weights_overall[4]!=0){
    node_correlations[[4]]<-mean(sapply(1:length(rep_1[[4]]), function(z)
      if((rep_1[[4]][[z]]==0) | (rep_2[[4]][[z]]==0)){
        0
      }else{
        1-(abs(rep_1[[4]][[z]]-rep_2[[4]][[z]]))
      }),na.rm=TRUE)
  }else{
    node_correlations[[4]]<-0
  }

  #rep vdj diversity #check if both same species. only compare VDJ if it is the case otherwise VDJ usage does not play a role as initial genes differ
  if(weights_overall[5]!=0){
    weights_V_D_J<-weights_VDJ
    if((rep_1[[7]]!=rep_2[[7]])|(rep_1[[8]]!=rep_2[[8]])){ #check if species are same and VDJ weights nonzero

      node_correlations[[5]]<-0

    } else if(sum(weights_V_D_J)==0){

      node_correlations[[5]]<-1

    } else{

      weights_V_D_J_norm<-weights_V_D_J/sum(weights_V_D_J)

      cor_V<-abs(cor(rep_1[[5]][[1]],rep_2[[5]][[1]],method=cor_method))
      cor_D<-abs(cor(rep_1[[5]][[2]],rep_2[[5]][[2]],method=cor_method))
      cor_J<-abs(cor(rep_1[[5]][[3]],rep_2[[5]][[3]],method=cor_method))
      cor_VJ<-abs(cor(rep_1[[5]][[4]],rep_2[[5]][[4]],method=cor_method))
      #calc weighted mean
      node_correlations[[5]]<-weighted.mean(c(cor_V,cor_D,cor_J,cor_VJ),weights_V_D_J_norm)
    }
  }else{
    node_correlations[[5]]<-0
  }

  #rep kmer_occurrence nt
  if(weights_overall[6]!=0){

    node_correlations[[6]]<-abs(cor(rep_1[[6]]$freq_counts_overall,rep_2[[6]]$freq_counts_overall,method=cor_method))

  }else{
    node_correlations[[6]]<-0
  }

  #weighted mean correlation across all categories
  wt <- weights_overall
  wt_norm<-wt/sum(wt)
  weighted_mean_correlation<-weighted.mean(unlist(node_correlations),wt_norm)

  if(is.na(weighted_mean_correlation)){
    weighted_mean_correlation<-0
  }

  return(weighted_mean_correlation)
}
