# Calculates positional AA frequency distribution
# params repertoire_df as obtained from load_mixcr_output
# returns list of lists of positional AA frequency of CDR3s for each length btw 6 and 28 AAs (list index1: length, index2:position)
#

.get_AA_freq_distribution<-function(repertoire_df,aa_length_range){
  #for lenghts 6:28 returns AA frequency list
  #each list contains a sublist for each length
  #which in turn contains a list AA frequencies per position.

  AA_list <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*")

  unique_lengths<-aa_length_range #c(8:20)

  list_AA_freqs_per_position<-list()
  for(i in 1:length(unique_lengths)){
    curr_length<-unique_lengths[[i]]
    curr_seqs<-as.character(repertoire_df[nchar(as.character(repertoire_df[,"junction_aa"]))==curr_length,"junction_aa"])

    list_AA_freqs_per_position_temp<-list()
    #check whether there are seqs of this length if yeas calc positional frequences
    if(length(curr_seqs)>0){
      for(j in 1:curr_length){
        pos_ <- sapply(1:length(curr_seqs), function(x) unlist(strsplit(curr_seqs[x],""))[j])
        pos_freqs <- table(factor(pos_,levels = AA_list))/length(pos_)
        list_AA_freqs_per_position_temp[[j]] <- pos_freqs
      }
    } else{ #otherwise set all to 0
      for(j in 1:curr_length){
        pos_ <- rep(0,length(AA_list))
        pos_freqs<-table(factor(pos_,levels=AA_list))
        list_AA_freqs_per_position_temp[[j]] <- pos_freqs
      }
    }
    names(list_AA_freqs_per_position_temp)<-c(1:curr_length)

    list_AA_freqs_per_position[[i]]<-list_AA_freqs_per_position_temp

  }
  names(list_AA_freqs_per_position)<-unique_lengths

  #return. first index: length, second index position.
  list_AA_freqs_per_position

}
