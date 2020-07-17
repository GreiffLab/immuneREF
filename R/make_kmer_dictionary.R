#' Preparatory step for kmer analysis makes kmer dictionary
#'
#' @param k_nt length of nucleotide pattern
#' @param k_aa length of amino acid pattern
#' @param gap_size_nt gap size for nt pattern
#' @param gap_size_aa gap size for aa pattern
#' @return Dictionary of gapped kmers
#' @examples
#' \dontrun{
#' make_kmer_dictionary(k_nt,k_aa,
#' gap_size_nt=c(0,1,2,3,4,5,6,7,8,9,10),
#' gap_size_aa=c(0,1,2,3))
#' }

# Preparatory step for kmer analysis makes kmer dictionary
# params k_nt
# params k_aa
# returns
#

make_kmer_dictionary<-function(k_nt,k_aa,gap_size_nt=c(0,1,2,3,4,5,6,7,8,9,10),gap_size_aa=c(0,1,2,3)){
  #makes a dictionary for kmer combo counts depending on kmer sizes used in analysis
  #only needs to be run once.

  cat('Preparing nt dictionary ...\n')

  k_prep_nt<-k_nt#c(1,2,3,4)
  dictionaries_nt<-list()
  for(i in 1:length(k_prep_nt)){
    k_length<-k_prep_nt[[i]]

    dictionary<-c("A","C","G","T")

    combn_it<-k_length
    all_combs<-dictionary
    while(combn_it>1){
      combn_dict<-base::merge(dictionary,all_combs)
      all_combs<-apply(combn_dict, 1, paste, collapse="")
      combn_it<-combn_it-1
    }#

    all_pairs_of_combs<-apply(base::merge(all_combs,all_combs),1,paste,collapse=",")
    nt_pair_dict_gen<-rep(0,length(all_pairs_of_combs))
    names(nt_pair_dict_gen)<-all_pairs_of_combs#

    dictionaries_nt[[i]]<-nt_pair_dict_gen
  }
  names(dictionaries_nt)<-as.character(k_nt)	#c("1","2","3","4")

  cat('nt dictionary Done \n')
  cat('Preparing aa dictionary ...\n')


  k_prep_aa<-k_aa #c(1,2)
  dictionaries_aa<-list()
  for(i in 1:length(k_prep_aa)){
    k_length<-k_prep_aa[[i]]

    dictionary<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*")

    combn_it<-k_length
    all_combs<-dictionary
    while(combn_it>1){
      combn_dict<-base::merge(dictionary,all_combs)
      all_combs<-apply(combn_dict, 1, paste, collapse="")
      combn_it<-combn_it-1
    }#

    all_pairs_of_combs<-apply(base::merge(all_combs,all_combs),1,paste,collapse=",")
    aa_pair_dict_gen<-rep(0,length(all_pairs_of_combs))
    names(aa_pair_dict_gen)<-all_pairs_of_combs#

    dictionaries_aa[[i]]<-aa_pair_dict_gen

  }
  names(dictionaries_aa)<-as.character(k_aa)

  dictionary_counts<-list()
  dictionary_counts[[1]]<-dictionaries_nt
  dictionary_counts[[2]]<-dictionaries_aa
  dictionary_counts[[3]] <- gap_size_nt
  dictionary_counts[[4]] <- gap_size_aa

  names(dictionary_counts)<-c("nt","AA","gap_nt","gap_aa")
  cat(' Done \n')

  return(dictionary_counts)
}
