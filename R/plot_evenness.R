#plots evenness curves. either for single repertoire or all repertoires in a list

##plot evenness for single repertoire or list
#plot_evenness(repertoires_analyzed[[1]],mode="single",var="evenness")
#plot_evenness(repertoires_analyzed[1:5],mode="list",var="evenness")
#plot_evenness(repertoires_analyzed[[1]],mode="single",var="diversity")
#plot_evenness(repertoires_analyzed[1:5],mode="list",var="diversity")

.plot_evenness<-function(repertoire,mode="single",var="evenness"){
  if(var == "diversity"){
    variable_to_plot<-"diversity"
    ytitle<-"Diversity"
  }else{
    variable_to_plot<-"evenness"
    ytitle<-"Evenness"
  }
  if(mode=="single"){

    evenness_df<-reshape2::melt(repertoire[[1]][[variable_to_plot]])
    evenness_df$x<-as.numeric(row.names(evenness_df))
    evenness_df<-evenness_df[evenness_df$x %in% seq(0,10,1),]

    evenness_plot <- ggplot2::ggplot(evenness_df, ggplot2::aes(x=x,y=value))+
      ggplot2::geom_line(size=1.5,alpha=0.6,show.legend=FALSE)+
      ggplot2::geom_point(show.legend=TRUE)+
      ggplot2::theme_bw()+
      ggplot2::labs(x = "alpha", y = ytitle,,colour="") +
      ggplot2::scale_x_continuous(breaks = seq(0,10,1))+
      ggplot2::scale_fill_manual(values = c("black"))#+
      #.theme.akbar()
  }else if(mode=="list"){
    list_evenness_df<-list()
    for(i in 1:length(repertoire)){
      evenness_df<-reshape2::melt(repertoire[[i]][[1]][[variable_to_plot]])
      evenness_df$x<-as.numeric(row.names(evenness_df))
      evenness_df$name<-repertoire[[i]][[9]]
      evenness_df<-evenness_df[evenness_df$x %in% seq(0,10,1),]
      list_evenness_df[[i]]<-evenness_df
    }
    evenness_plot_df<-data.table::setDF(data.table::rbindlist(list_evenness_df))

    evenness_plot_df$name <- factor(evenness_plot_df$name,levels=names(repertoire))

    evenness_plot <- ggplot2::ggplot(evenness_plot_df, ggplot2::aes(x=x,y=value,colour=name))+
      ggplot2::geom_line(size=1.5,alpha=0.6,show.legend=FALSE)+
      ggplot2::geom_point(show.legend=TRUE)+
      ggplot2::theme_bw()+
      ggplot2::scale_x_continuous(breaks = seq(0,10,1))+
      ggplot2::scale_colour_manual(values = c("red","royalblue3"),
                                 guide = ggplot2::guide_legend(
                                  direction = "horizontal",
                                  label.position = "bottom",
                                  override.aes = list(size=6),
                                  label.theme = ggplot2::element_text(size=15,angle = 0)))+
      ggplot2::labs(x = "Alpha", y = ytitle,colour="")+
      ggplot2::theme(legend.position="bottom")+
      ggplot2::theme(axis.title.x = ggplot2::element_text(size=20),axis.title.y = ggplot2::element_text(size=20),axis.text.x = ggplot2::element_text(size=15, angle = 0),axis.text.y = ggplot2::element_text(size=20, angle = 0))
      #ggplot2::scale_fill_manual(values = c("black"))+
      #.theme.akbar()
  }

  return(evenness_plot)

}






