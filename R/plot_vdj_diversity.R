#RETURNS VDJ frequency plots either for one repertoire
# or for a list of repertoires (need to plot using wrapping)

#.plot_vdj_diversity(repertoires_analyzed[[1]],"single")
#.plot_vdj_diversity(repertoires_analyzed[1:5],"list")

.plot_vdj_diversity<-function(repertoire,chain="H",mode="single"){


  if(mode=="single"){
    #make ready for plot by ordering V,D,J
    vdj_occurrence<-repertoire[[5]][1:3]
    names(vdj_occurrence)<-c("V","D","J")

    names(vdj_occurrence[["V"]])[names(vdj_occurrence[["V"]])==""]<-"V_gene_unkown"
    names(vdj_occurrence[["D"]])[names(vdj_occurrence[["D"]])==""]<-"D_gene_unkown"
    names(vdj_occurrence[["J"]])[names(vdj_occurrence[["J"]])==""]<-"J_gene_unkown"

    list_order_genes<-list()
    list_order_genes[["V"]]<-names(sort(-vdj_occurrence[["V"]]))
    list_order_genes[["D"]]<-names(sort(-vdj_occurrence[["D"]]))
    list_order_genes[["J"]]<-names(sort(-vdj_occurrence[["J"]]))


    vdj_occurrence_plot_df <- reshape2::melt(vdj_occurrence)
    names(vdj_occurrence_plot_df)<-c("Var1","Freq","gene")

    vdj_occurrence_plot_df$gene <- factor(vdj_occurrence_plot_df$gene,levels=c("V","D","J"))

    #for each of V,D,J create plot
    vdj_occurrence_plot_list<-list()
    for(i in 1:length(unique(vdj_occurrence_plot_df$gene))){
      curr_vdj_occurrence_plot_df <- vdj_occurrence_plot_df[vdj_occurrence_plot_df$gene==unique(vdj_occurrence_plot_df$gene)[i],]

      curr_vdj_occurrence_plot_df$Var1<-factor(curr_vdj_occurrence_plot_df$Var1,levels=list_order_genes[[unique(vdj_occurrence_plot_df$gene)[i]]])
      vdj_occurrence_plot_list[[i]]<-curr_vdj_occurrence_plot_df
    }

    vdj_occurrence_plot_all<-data.table::setDF(data.table::rbindlist(vdj_occurrence_plot_list))
    #add variable definitions pre plotting.
    Var1 <- Freq <- NULL

    #name_plot <- paste(folder,"/vdj_occurrence_",unique(vdj_occurrence_plot_df$gene)[i],"_",vdj_occurrence$repertoire,"_axis.pdf",sep="")
    vdj_occurrence_plot_out <- ggplot2::ggplot(data=vdj_occurrence_plot_all, ggplot2::aes(Var1, Freq))+
      ggplot2::geom_bar(stat="identity",position="dodge") +
      ggplot2::facet_wrap(.~gene, scales="free")+
      ggplot2::theme_bw()+
      ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
      .theme.akbar()+
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=2, angle = 90))

    list_plots_vdj<-vdj_occurrence_plot_out

  }else if(mode=="list"){

    list_vdj_rep<-list()
    for(z in 1:length(repertoire)){

      repertoire_name<-names(repertoire)[z]

      vdj_occurrence<-repertoire[[z]][[5]][c(1,3)]
      names(vdj_occurrence)<-c("V","J")

      names(vdj_occurrence[["V"]])[names(vdj_occurrence[["V"]])==""]<-"V_gene_unkown"
      #names(vdj_occurrence[["D"]])[names(vdj_occurrence[["D"]])==""]<-"D_gene_unkown"
      names(vdj_occurrence[["J"]])[names(vdj_occurrence[["J"]])==""]<-"J_gene_unkown"

      list_order_genes<-list()
      list_order_genes[["V"]]<-names(sort(-vdj_occurrence[["V"]]))
      #list_order_genes[["D"]]<-names(sort(-vdj_occurrence[["D"]]))
      list_order_genes[["J"]]<-names(sort(-vdj_occurrence[["J"]]))


      vdj_occurrence_plot_df <- reshape2::melt(vdj_occurrence)
      names(vdj_occurrence_plot_df)<-c("Var1","Freq","gene")

      vdj_occurrence_plot_df$gene <- factor(vdj_occurrence_plot_df$gene,levels=c("V","J"))

      list_vdj_rep[[z]]<-vdj_occurrence_plot_df
    }

    vdj_maxmin<-merge(list_vdj_rep[[1]],list_vdj_rep[[2]],by="Var1")
    vdj_maxmin$gene.y<-NULL
    names(vdj_maxmin)<-c("gene_name","freq_A","gene","freq_B")
    vdj_maxmin<-vdj_maxmin[vdj_maxmin$freq_A>=0.01 | vdj_maxmin$freq_B>=0.01,]


    vdj_maxmin_A<-vdj_maxmin[,c("gene_name","gene","freq_A")]
    vdj_maxmin_B<-vdj_maxmin[,c("gene_name","gene","freq_B")]
    names(vdj_maxmin_A)<-c("gene_name","gene","freq")
    vdj_maxmin_A$repertoire<-"max"
    names(vdj_maxmin_B)<-c("gene_name","gene","freq")
    vdj_maxmin_B$repertoire<-"min"

    vdj_maxmin_AB<-rbind(vdj_maxmin_A,vdj_maxmin_B)
    names(vdj_maxmin_AB)<-c("Gene","VJ","Freq","Category")

    vdj_maxmin_v<-vdj_maxmin_AB[vdj_maxmin_AB$VJ=="V",]
    vdj_maxmin_j<-vdj_maxmin_AB[vdj_maxmin_AB$VJ=="J",]

    vdj_occurrence_plot_comp_J <- ggplot2::ggplot(data=vdj_maxmin_j, ggplot2::aes(Gene, Freq,fill=Category))+
      ggplot2::geom_bar(stat="identity",position="dodge",colour='black',show.legend=FALSE) +
      ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("red","royalblue3"))+#3rd pos))+
      ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, angle = 90,vjust=0.5),axis.text.y = ggplot2::element_text(size=20),axis.title.y = ggplot2::element_text(size=20))

    vdj_occurrence_plot_comp_V <- ggplot2::ggplot(data=vdj_maxmin_v, ggplot2::aes(Gene, Freq,fill=Category))+
      ggplot2::geom_bar(stat="identity",position="dodge",colour='black',show.legend=FALSE) +
      ggplot2::theme_bw()+
      ggplot2::scale_fill_manual(values=c("red","royalblue3"))+#3rd pos))+
      ggplot2::labs(x = "", y = "Frequency (%)", fill = "", colour = "Repertoire") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, angle = 90,vjust=0.5),axis.text.y = ggplot2::element_text(size=20),axis.title.y = ggplot2::element_text(size=20))

    list_plots_vdj<-list()
    list_plots_vdj[["V"]]<-vdj_occurrence_plot_comp_V
    list_plots_vdj[["J"]]<-vdj_occurrence_plot_comp_J

  }

  return(list_plots_vdj)
}






