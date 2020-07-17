#' Prints violin plots of similarity values per layer, across categories
#'
#' @param list_similarity_matrices list of immuneREF layers
#' @param categories_list List with entries 'categories' and 'color' and 'subset' for plot structure
#' @param path_figure Path where plots should be stored
#' @return TRUE
#' @examples
#' \dontrun{
#' print_global_similarity(list_similarity_matrices,categories_list,path_figure="figures")
#' }


print_global_similarity<-function(list_similarity_matrices,categories_list,path_figure="figures"){

  categories_df<-categories_list[["categories"]]
  color_categories<-categories_list[["color"]]
  subsets <- categories_list[["subset"]]

  for(i in 1:length(list_similarity_matrices)){

    cormat<-list_similarity_matrices[[i]]
    current_layer<-names(list_similarity_matrices)[i]

    #reduce to upper triangle of matrix and reshape into dataframe
    cormat_overall<-reshape2::melt(cormat[upper.tri(cormat)])
    cormat_overall$condition<-"All"

    unique_subsets<-unique(categories_df[,subsets])

    list_of_subsets<-list()
    for(z in 1:length(unique_subsets)){

      rows_current_subset<-categories_df[categories_df[[subsets]] == unique_subsets[z],"sample_id"]

      curr_matrix<-cormat[rows_current_subset,rows_current_subset]

      curr_matrix_upper<-reshape2::melt(curr_matrix[upper.tri(curr_matrix)])
      curr_matrix_upper$condition<-unique_subsets[z]

      list_of_subsets[[z]]<-curr_matrix_upper
    }

    list_of_subsets[[z+1]]<-cormat_overall
    bound_subsets<-do.call(rbind,list_of_subsets)

    #summary_stats<-summarySE(bound_subsets,measurevar = "value",groupvars=c("name"))

    plot_global_similarity <- ggplot2::ggplot(bound_subsets, ggplot2::aes(x=condition,y=value,fill=condition))+                    # basic graphical object
      ggplot2::geom_violin()+
      ggplot2::geom_boxplot(width=0.1, outlier.size = 0.25)+
      ggplot2::scale_fill_manual(values = color_categories)+
      ggplot2::labs(x="",y="Similarity score",fill="Condition")+
      ggplot2::theme_bw()+
      ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=ggplot2::rel(1), hjust=0), plot.margin = grid::unit(c(30, 2, 2, 2),"points"), axis.text.x = ggplot2::element_text(size=12, angle = 0),  axis.text.y = ggplot2::element_text(size=14), axis.title= ggplot2::element_text(size=14, vjust = 0.4)) +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank(), panel.border = ggplot2::element_rect(colour = "black"))+
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"), strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.8)),  panel.grid.major.x = ggplot2::element_blank())

      grDevices::pdf(file.path(path_figure,paste("global_similarity_",current_layer,".pdf",sep="")), width = 7,  height = 5.5,onefile=FALSE)
      grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
      vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
      base::print(plot_global_similarity, vp=vplayout(1,1))   # plot for row 1, column 1
      grDevices::dev.off()
  }

  return(TRUE)
}

