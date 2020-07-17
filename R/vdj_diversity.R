# Calculates VDJ diversity
# params repertoire_df as obtained from load_mixcr_output
# params list_vdjgenes: reference for vdj names contains BCR, TCR data from human and mice from IMGT. IMPORTANT: VDJ names in input need to be same as in this vdj list
# params species: defines what species data is from 'hs' (human) or 'mm' (mouse)
# params receptor: defines whether data is from BCR: 'ig' or TCR: 'tr' data
# params chain: defines heavy,light for ig either 'h','k' or 'l'. for tr: either 'b' or 'a'
# params repertoire_df as obtained from load_mixcr_output
# returns list with 3 entries (V,D,J) for which each
#

.vdj_diversity<-function(repertoire_df,list_vdjgenes,species,receptor,chain){
  #returns list with three entries containing frequencies of each V,D and J usage.
  #first checks which species, repertoires are from and loads corresponding vdj list
  #then extracts main info about v,d,j from repertoire (f.e. IGHV1-86)

  #load appropriate VDJ information
  vdj_genes<-list_vdjgenes[[species]][[receptor]][[chain]]

  #extract essential V,D,J name portions and table it
  v_genes_uncut<-as.character(repertoire_df$v_call)
  #v_genes_present<-gsub("*[*].*","",v_genes_uncut)
  v_genes_present<-v_genes_uncut
  v_genes<-table(factor(v_genes_present,levels=unique(c(as.character(vdj_genes$V$gene),""))))
  v_genes_freq<-v_genes/sum(v_genes)

  if(paste(receptor,chain,sep="") %in% c("igh","trb")){
    d_genes_uncut<-as.character(repertoire_df$d_call)
    #d_genes_present<-gsub("*[*].*","",d_genes_uncut)
    d_genes_present<-d_genes_uncut

    if("" %in% vdj_genes$D$gene){
      d_genes<-table(factor(d_genes_present,levels=unique(c(as.character(vdj_genes$D$gene)))))
    }else{
      d_genes<-table(factor(d_genes_present,levels=unique(c(as.character(vdj_genes$D$gene),""))))
    }
    #d_genes<-table(factor(d_genes_present,levels=c(as.character(vdj_genes$D$gene),"")))
    d_genes_freq<-d_genes/sum(d_genes)
  }else{
    d_genes_freq<-NA
  }

  j_genes_uncut<-as.character(repertoire_df$j_call)
  #if(receptor=="ig"){
  #	j_genes_present<-gsub("*[-].*","",j_genes_uncut)
  #}else{
  #	j_genes_present<-j_genes_uncut
  #}
  j_genes_present<-j_genes_uncut
  j_genes<-table(factor(j_genes_present,levels=unique(c(as.character(vdj_genes$J$gene),""))))
  j_genes_freq<-j_genes/sum(j_genes)

  #save all vj -ocmbos
  v_j_combos_present <- paste(v_genes_present,j_genes_present,sep="_")
  df_vj <- expand.grid(V = as.character(vdj_genes$V$gene), J = as.character(vdj_genes$J$gene))
  level_vj <- paste(df_vj$V,df_vj$J,sep="_")
  vj_genes<-table(factor(v_j_combos_present,levels=unique(c(level_vj,""))))
  vj_genes_freq <-vj_genes/sum(vj_genes)

  list_vdj_diversity<-list()
  list_vdj_diversity[[1]]<-v_genes_freq
  list_vdj_diversity[[2]]<-d_genes_freq
  list_vdj_diversity[[3]]<-j_genes_freq
  list_vdj_diversity[[4]]<-vj_genes_freq

  names(list_vdj_diversity)<-c("v","d","j","vj")

  list_vdj_diversity
}
