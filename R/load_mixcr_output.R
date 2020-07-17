#' Loads output from MiXCR
#'
#' @param filename Address of file to be loaded
#' @param AA_length_of_seqs if !=0 will only include sequences of one specified length

#' @return Processed repertoires
#' @examples
#' \dontrun{
#' load_mixcr_output(filename,AA_length_of_seqs)
#' }

# Loads regular output from MiXCR.
# params filename: address of file
# params AA_length_of_seqs: if !=0 will only include sequences of one specified length
# returns repertoire dataframe ready for immuneNET processing
#

load_mixcr_output<-function(filename,AA_length_of_seqs=0){
  #reads in txt file (mixcr output) and returns df
  #if a length is specified returns only sequences of the specified lenght
  #otherwise will return full df
  repertoire_df_all_lengths<-read.table(filename,header=TRUE,dec=".",sep="\t")

  if(AA_length_of_seqs==0){
    repertoire_df_return<-repertoire_df_all_lengths
  } else{
    repertoire_df_return<-repertoire_df_all_lengths[nchar(as.character(repertoire_df_all_lengths$AA..Seq.CDR3))==AA_length_of_seqs & nchar(as.character(repertoire_df_all_lengths$N..Seq..CDR3))==(AA_length_of_seqs*3),]
  }

  cat("ADD NEW NAMES!!!")

  return(repertoire_df_return)
}
