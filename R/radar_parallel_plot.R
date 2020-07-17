#RETURNS Spider and parallel plot
# or for a list of repertoires (need to plot using wrapping)

.radar_parallel_plot<-function(list_similarity_matrices,to_compare,name_plot="sim"){

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

  similarity_values_df$roi<-to_compare[["roi_names"]]

  radar_plot<-ggiraphExtra::ggRadar(data=similarity_values_df,aes(x=c(diversity,aa_freq,architecture,overlap,vdj_diversity,kmers),group=roi),
                                    rescale=FALSE,
                                    scales="free",
                                    size = 1,
                                    alpha=0.2,
                                    interactive=FALSE,
                                    legend.position="right")+
    theme_bw()

  parallel_plot<-GGally::ggparcoord(similarity_values_df,
                                    columns=1:6,
                                    groupColumn = "roi",
                                    alphaLines = 0.8)+
    theme_bw()


  pdf(paste("figures/radar_parallel/reference_",name_plot,"_radar_parallel.pdf",sep=""), width = 14,  height = 6,onefile=FALSE)
  pushViewport(viewport(layout=grid.layout(1,2), height=1,width=1)) # 3 rows, 1 columns
  vplayout<-function(x,y)viewport(layout.pos.row=x,layout.pos.col=y)
  print(radar_plot, vp=vplayout(1,1))   # plot for row 1, column 1
  print(parallel_plot, vp=vplayout(1,2))   # plot for row 1, column 1
  dev.off()


}
