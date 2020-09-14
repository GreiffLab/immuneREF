
.immunosignature_presence<-function(repertoire_df,models){
  #runs a previously trained ML model on the data and return freq of 'positively'
  #classified sequences
  ##for_example number of predicted public_clones based on svm model from our JI publication Greiff2017

  perc_classification_list<-list()
  for(i in 1:length(models)){
    model<-models[[i]]
    test_set<-Biostrings::DNAStringSet(repertoire_df$junction)
    pred <- kebabs::predict(model, test_set)
    perc_class<-table(pred)[[1]]/sum(table(pred))
    #percPrivate in case of publicprivate model
    perc_classification_list[[i]]<-perc_class
  }
  perc_classification_list
}
