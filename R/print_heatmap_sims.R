#' Prints heatmap of all immuneREF layers in input list
#'
#' @param list_similarity_matrices list of immuneREF layers
#' @param annotation_list List with entries 'categories' and 'colors' for heatmap annotation
#' @param path_figure Path where plots should be stored
#' @return TRUE
#' @examples
#' \dontrun{
#' print_heatmap_sims(list_similarity_matrices,annotation_list,path_figure="figures")
#' }

#
print_heatmap_sims<-function(list_similarity_matrices,annotation_list,path_figure="figures"){

  for(i in 1:length(list_similarity_matrices)){

    cormat<-list_similarity_matrices[[i]]
    layer_name <- names(list_similarity_matrices)[i]
  
    #make categories dataframe  
  
    ha = ComplexHeatmap::HeatmapAnnotation(df=annotation_list[["categories"]],
                           col = annotation_list[["colors"]],
                           show_annotation_name=TRUE,
                           annotation_name_side = "right",
                           height = grid::unit(10, "cm"),
                           gap = grid::unit(0.03, "cm"),
                           simple_anno_size = grid::unit(1.5, "cm")
     )
  
  
    ht_map<-ComplexHeatmap::Heatmap(cormat,
                    name = "Similarity", 
                    col = circlize::colorRamp2(c(0, 0.5,1), c("blue","yellow","red")),
                    top_annotation = ha,             
                    row_names_gp = grid::gpar(fontsize = 0),
                    column_names_gp = grid::gpar(fontsize = 0)
                    )
    
    grDevices::pdf(file.path(path_figure, paste("heatmap_",layer_name,".pdf",sep="")), width = 12.5,  height = 10,onefile=FALSE)
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,1), height=1,width=1)) # 3 rows, 1 columns
    vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
    ComplexHeatmap::draw(ht_map,merge_legend = FALSE, heatmap_legend_side = "right", annotation_legend_side = "right")
    grDevices::dev.off()
  }
  return(TRUE)

}

