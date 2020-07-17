#' CALCULATES IMMUNEREF LAYER CHARACTERISTICS
#'
#' @param repertoire_df Input repertoire (AIRR standard column naming)
#' @param models List of ML models compatible with kebabs::predict() function
#' @param species Species information of input repertoire ("hs" or "mm")
#' @param receptor Receptor information of input repertoire ("ig" or "tr")
#' @param chain Chain information of input repertoire ("h","l","b" or "a")
#' @param identifier_rep ID-name for input repertoire
#' @param vdj_list Reference list for germline gene naming (provided dataset: list_germline_genes)
#' @param dictionary Dictionary for kmer-occurrence search (provided: dictionary_counts)
#' @param architecture_depth maximal LD for network calculation (a network is calculated for each distance 1:architecture_depth)
#' @param aa_range range of sequence lengths for which amino acid frequencies should be evaluated
#' @return List of repertoire layers
#' @examples
#' \dontrun{
#' calc_characteristics(repertoire_df,models=list(svm_public_mm=svm_model_public),species,receptor,chain,identifier_rep,vdj_list=list_germline_genes,dictionary=dictionary_counts)
#' }

calc_characteristics<-function(repertoire_df,models=list(svm_public_mm=svm_model_public),species,receptor,chain,identifier_rep,vdj_list=list_germline_genes,dictionary=dictionary_counts,architecture_depth=1,aa_range=c(8:20)){
  #call all repertoire analysis fcts and
  #save all results in a characteristics list.

  #calcualte each repertoire features and output progress message in btw.
  cat("Current repertoire: ",identifier_rep, "\n")
  cat("Begin calculating repertoire features \n")
  repertoire_characteristics<-list()
  cat("1/6 Calculating Diversity \n")
  repertoire_characteristics[[1]]<-.get_evenness_diversity_distribution(repertoire_df)
  cat("2/6 Calculating Amino acid frequency \n")
  repertoire_characteristics[[2]]<-.get_AA_freq_distribution(repertoire_df, aa_length_range=aa_range)
  cat("3/6 Calculating Repertoire architecture \n")
  repertoire_characteristics[[3]]<-.architecture_charac(repertoire_df,max_ld=architecture_depth,option="AA")
  #repertoire_characteristics[[3]]<-.architecture_coverage(repertoire_df,max_ld=6,option="AA")
  cat("4/6 Calculating Convergence \n")
  repertoire_characteristics[[4]]<-.immunosignature_presence(repertoire_df,models)
  cat("5/6 Calculating Germline gene usage \n")
  repertoire_characteristics[[5]]<-.vdj_diversity(repertoire_df,vdj_list,species,receptor,chain)
  cat("6/6 Calculating K-mer occurrence \n")
  #extract parameters k and m from dictionary
  k_mers_in_dict <- as.integer(names(dictionary[["nt"]]))
  gap_sizes <- dictionary[["gap_nt"]]
  repertoire_characteristics[[6]]<-.kmer_occurrence_opt_k_m(repertoire_df,"nt",k_mers_in_dict,gap_sizes,dictionary)
  #save additional information on repertoire (to make it easily accessible in downstream analysis)
  repertoire_characteristics[[7]]<-species #species for comparison step downstream
  repertoire_characteristics[[8]]<-paste(receptor,chain,sep="") #species for comparison step downstream
  repertoire_characteristics[[9]]<-identifier_rep
  repertoire_characteristics[[10]]<-length(repertoire_df$junction)

  names(repertoire_characteristics)<-c("diversity","AA_freq","average_degree","immunosignatures","vdj_diversity","kmers_nt","species","receptor","repertoire","number_of_seqs")
  cat("Repertoire features calculated\n")

  return(repertoire_characteristics)

}
