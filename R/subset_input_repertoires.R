#' Subsets list of input repertoires to desired nb of sequences.
#'
#' @param list_repertoires list of input repertoires
#' @param subset_size maximal number of sequences per repertoire (default 10'000)
#' @param random If TRUE subsetting is done by random subsampling, else top x clones are picked (by frequency)
#' @return List of subsampled repertoires features
#' @examples
#' \dontrun{
#' subset_input_repertoires(tutorial_repertoires,100, FALSE)
#' }

subset_input_repertoires <- function(list_repertoires,subset_size=10000,random=FALSE){

  for(i in 1:length(list_repertoires)){
    curr_repertoire<-list_repertoires[[i]]

    #check if needs to be subsampled
    if(nrow(curr_repertoire) > subset_size){
      if(random==FALSE){

        #order sequences in repertoire by frequencies
        curr_repertoire <- curr_repertoire[order(-curr_repertoire$freqs),]
        #pick top x sequences
        curr_repertoire <- curr_repertoire[c(1:subset_size), ]

      }else if(random == TRUE){
        #randomly choose x rows of sequences
        curr_repertoire <- curr_repertoire[sample(c(1:nrow(curr_repertoire)), subset_size), ]        
      }

    }
  
    list_repertoires[[i]] <- curr_repertoire
  }

  return(list_repertoires)

}
