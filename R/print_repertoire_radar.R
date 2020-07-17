#' Prints violin plots of similarity values per layer, across categories
#'
#' @param list_similarity_matrices list of immuneREF layers
#' @param to_compare List with entries 'roi' (repertoires of interest), 'roi names', 'ref',' 'plot_names' and 'colors'
#' @param path_figure Path where plots should be stored
#' @param name_plot Name of plot
#' @return TRUE
#' @examples
#' \dontrun{
#' print_repertoire_radar(list_similarity_matrices,to_compare,path_figure,name_plot)
#' }
#
print_repertoire_radar<-function(list_similarity_matrices,to_compare,path_figure,name_plot){

  similarity_values<-list()
  for(z in 1:length(to_compare[["roi"]])){
    similarity_value_tmp<-list()
    for(i in 1:length(list_similarity_matrices)){
      cormat<-list_similarity_matrices[[i]]

      similarity_value_tmp[[names(list_similarity_matrices)[i]]] <-cormat[to_compare[["roi"]][z],to_compare[["ref"]]]
    }

    #similarity_value_tmp_df<-reshape2::melt(t(similarity_value_tmp))
    similarity_value_tmp_df<-as.data.frame(do.call(cbind,similarity_value_tmp))

    similarity_value_tmp_df$roi<-to_compare[["roi"]][z]
    similarity_value_tmp_df$ref<-to_compare[["ref"]]

    similarity_values[[to_compare[["roi"]][z]]]<-similarity_value_tmp_df
  }

  similarity_values_df<-do.call(rbind,similarity_values)
  rownames(similarity_values_df)<-c(1:nrow(similarity_values_df))

  similarity_values_df$Cohort<-to_compare[["roi_names"]]

  color_vec_fill<-to_compare[["colors"]]

  color_list <- list()
  for(i in 1:length(color_vec_fill)){
    color_list[[i]] <- .darken_color(color_vec_fill[i],factor=1.4)
  }
  color_vec_frame<- unlist(color_list)


  radar_plot<-ggiraphExtra::ggRadar(data=similarity_values_df,ggplot2::aes(x=c(Diversity,AAfreq,Architecture,Convergence,VDJ_usage,k-mers),group=Cohort,colour=Cohort,fill=Cohort),
    rescale=FALSE,
    scales="free",
    size = 1,
    alpha=0.4,
    interactive=FALSE,
    use.label=FALSE,
    legend.position="right")+#,
   # axis.label.size = 1, grid.label.size = 1, legend.text.size = 10)+
  ggplot2::theme_bw()+
  ggplot2::scale_fill_manual(values = color_vec_fill)+
  ggplot2::scale_colour_manual(values = color_vec_frame)
  #theme_minimal() +
  #theme(text = element_blank(), # custom font size
  #  axis.text.y = element_blank())

  grDevices::pdf(file.path(path_figure,paste("radar_plot_",name_plot,".pdf",sep="")), width = 11,  height = 8,onefile=FALSE)
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
  base::print(radar_plot, vp=vplayout(1,1))   # plot for row 1, column 1
  grDevices::dev.off()

}
