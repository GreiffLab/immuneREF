##Plots component size and number of occurrences of component of each size
# Options: single or list (either one repertoire or many)
#         which variable should be plotted: for now only components (later add across LD)
#         log_plot option: decides if plot should be plotted with loglog axes.

#.plot_arch(repertoires_analyzed[[1]],mode="single",plot_var="components",log_plot=TRUE)
#.plot_arch(repertoires_analyzed[1:5],mode="list",plot_var="components",log_plot=TRUE)


.plot_arch<-function(repertoire,mode="single",plot_var="components",log_plot=TRUE){

  if(plot_var=="components"){
    variable_to_plot<-"component_size_and_members"

    x_var_name<-"Component size (number of nodes)"
    y_var_name<-"Number of occurrences"

  }else{
    cat("Error: Plotting for this variable is not yet available\n")
    return()
  }

  if(mode=="single"){

    components_members<-reshape2::melt(repertoire[[3]][[1]][[variable_to_plot]])
    names(components_members)<-c("component_size","occurrences")

    if(log_plot==TRUE){
      components_members$component_size<-log(components_members$component_size)
      components_members$occurrences<-log(components_members$occurrences)

      x_var_name<-paste("log(",x_var_name,")",sep="")
      y_var_name<-paste("log(",y_var_name,")",sep="")

    }

    arch_plot <- ggplot2::ggplot(components_members, ggplot2::aes(x=component_size,y=occurrences))+
      ggplot2::geom_line()+
      ggplot2::theme_bw()+
      ggplot2::labs(x = x_var_name, y = y_var_name) +
      ggplot2::scale_fill_manual(values = c("black"))+
      .theme.akbar()

  }else if(mode=="list"){
    list_arch_df<-list()
    for(i in 1:length(repertoire)){

      components_members<-reshape2::melt(repertoire[[i]][[3]][[1]][[variable_to_plot]])
      names(components_members)<-c("component_size","occurrences")
      components_members$name<-repertoire[[i]][[9]]

      if(log_plot==TRUE){
        components_members$component_size<-log(components_members$component_size)
        components_members$occurrences<-log(components_members$occurrences)

        x_var_name_plot<-paste("log(",x_var_name,")",sep="")
        y_var_name_plot<-paste("log(",y_var_name,")",sep="")

      }

      list_arch_df[[i]]<-components_members

    }
    arch_plot_df<-data.table::setDF(data.table::rbindlist(list_arch_df))

    arch_plot <- ggplot2::ggplot(arch_plot_df, ggplot2::aes(x=component_size,y=occurrences,colour=name))+
      ggplot2::geom_line(alpha=0.4,show.legend=FALSE)+
      ggplot2::theme_bw()+
      ggplot2::labs(x = x_var_name_plot, y = y_var_name_plot) +
      #ggplot2::scale_fill_manual(values = c("black"))+
      .theme.akbar()
  }

  return(arch_plot)

}
