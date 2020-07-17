#' Function to recalculate single immuneREF layer
#'
#' @param analyzed_list prepared and named repertoires_analyzed list
#' @param characteristic layer of interest
#' @param repertoire_df repertoire that needs to recalculated
#' @param models immunosignature models
#' @param species species of repertoire
#' @param receptor receptor of repertoire
#' @param chain chain of repertoire
#' @param identifier_rep name of repertoire
#' @param vdj_list reference germline gene list
#' @param dictionary dictionary for kmer calculation
#' @param architecture_depth maximal LD for network calculation (a network is calculated for each distance 1:architecture_depth)
#' @param aa_range range of sequence lengths for which amino acid frequencies should be evaluated
#' @return List of repertoire layers
#' @examples
#' \dontrun{
#' recalc_single_characteristic(analyzed_list,characteristic,
#' repertoire_df, svm, "mm", "igh", "repertoire_1",
#' list_germline_genes_allele_01, kmer_dict)
#' }


#if a single characteristic is faulty this function allows for reevaluation of a single characteristic (no need to run full analysis)
recalc_single_characteristic<-function(analyzed_list,characteristic,repertoire_df,models,species,receptor,chain,identifier_rep,vdj_list,dictionary,architecture_depth=1,aa_range=c(8:20)){
  #call all repertoire analysis fcts and
  #save all results in a characteristics list.

  if(characteristic=="Evenness"){
    cat("1 start evenness \n")
    analyzed_list[[1]]<-.get_evenness_diversity_distribution(repertoire_df)
  }
  if(characteristic=="AA_frequency"){
    cat("2 start AA freq \n")
    analyzed_list[[2]]<-.get_AA_freq_distribution(repertoire_df,aa_length_range=aa_range)
  }
  if(characteristic=="Architecture"){
    cat("3 start arch \n")
    analyzed_list[[3]]<-.architecture_charac(repertoire_df,max_ld=architecture_depth,option="AA")
  }
  if(characteristic=="Immunosignatures"){
    cat("4 start immunosig \n")
    analyzed_list[[4]]<-.immunosignature_presence(repertoire_df,models)
  }
  if(characteristic=="VDJ_usage"){
    cat("5 start vdj div \n")
    analyzed_list[[5]]<-.vdj_diversity(repertoire_df,vdj_list,species,receptor,chain)
  }
  if(characteristic=="kmer_occurrence"){
    cat("6 start kmer nt \n")
    k_mers_in_dict <- as.integer(names(dictionary[["nt"]]))
    gap_sizes <- dictionary[["gap_nt"]]
    #    analyzed_list[[6]]<-kmer_occurrence(repertoire_df,"nt",k_mers_in_dict,gap_sizes,dictionary_counts)
    analyzed_list[[6]]<-.kmer_occurrence_opt_k_m(repertoire_df,"nt",k_mers_in_dict,gap_sizes,dictionary)
  }
  if(characteristic=="architecture_coverage"){
    cat("3 start arch_coverage \n")
    analyzed_list[[3]]<-.architecture_coverage(repertoire_df,max_ld=architecture_depth,option="AA")
  }

  names(analyzed_list)<-c("diversity","AA_freq","average_degree","immunosignatures","vdj_diversity","kmers_nt","species","receptor","repertoire","number_of_seqs")

  analyzed_list

}
