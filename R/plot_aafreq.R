#RETURNS Amino Acid frequency plot either for one repertoire
# or for a list of repertoires (need to plot using wrapping)

#plot_aafreq(repertoires_analyzed[[1]],mode="single",length_seq='14')
#plot_aafreq(repertoires_analyzed[1:5],mode="list",length_seq='14')

.plot_aafreq<-function(repertoire,mode="single",length_seq='14'){

  colorset<-c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','black')
  color_list <- list()
  for(i in 1:length(colorset)){
    color_list[[i]] <- .darken_color(colorset[i],factor=1.)
  }
  darkened_colors<- unlist(color_list)


  curr_length<-length_seq

  if(mode=="single"){

    #calculate AA freq for AA length of interest
    #select the frequencies of that length
    AA_freqs_input<-repertoire[[2]][[length_seq]]

    #melt into dataframe and make the values into %
    AA_freq_plot <- reshape2::melt(AA_freqs_input)
    AA_freq_plot$value <- AA_freq_plot$value*100

    #name repertoire and plot
    AA_freq_plot$name_repertoire <- repertoire[[9]]
    names(AA_freq_plot) <- c('amino_acid','freqs','position','name_repertoire')
    number_of_positions <- length(AA_freq_plot$amino_acid)/21
    #set position for labels
    AA_freq_plot <- base::transform(AA_freq_plot, mid_y = stats::ave(AA_freq_plot$freqs, AA_freq_plot$position, FUN = function(val) base::cumsum(val) - (0.5 * val)))
    AA_freq_plot$mid_y <- 100-AA_freq_plot$mid_y

    #initialize variables
    position <- amino_acid <- mid_y <- freqs <-  NULL

    #plot
    aa_freq_curr_length <- ggplot2::ggplot(AA_freq_plot, ggplot2::aes(x=position,y=freqs,fill=amino_acid,label=amino_acid))+
      ggplot2::geom_bar(stat="identity",show.legend = FALSE)+
      ggplot2::geom_text(ggplot2::aes(y = mid_y),size=1,colour='black',fontface='bold',show.legend=FALSE)+
      ggplot2::geom_hline(yintercept=c(25,50,75),linetype="dotted")+
      ggplot2::theme_bw()+
      ggplot2::labs(x = "CDR3 sequence position", y = "Occurrence AA [%]", fill = "", colour = "") +
      #ggtitle(paste("AA composition per position (length=",unique_lengths[i],")"))+
      ggplot2::scale_fill_manual(values = colorset)+
      ggplot2::scale_colour_manual(values = darkened_colors) +
      ggplot2::scale_x_discrete(limit = factor(c(1:number_of_positions)))+
      .theme.akbar()

    list_aa_plots <- aa_freq_curr_length

  }else if(mode=="list"){

    list_aa_plots<-list()
    for(i in 1:length(repertoire)){

      AA_freqs_input<-repertoire[[i]][[2]][[length_seq]]
      #melt into dataframe and make the values into %
      AA_freq_plot <- reshape2::melt(AA_freqs_input)
      AA_freq_plot$value <- AA_freq_plot$value*100

      #name repertoire and plot
      AA_freq_plot$name_repertoire <- repertoire[[i]][[9]]
      names(AA_freq_plot) <- c('amino_acid','freqs','position','name_repertoire')
      number_of_positions <- length(AA_freq_plot$amino_acid)/21
      #set position for labels
      AA_freq_plot <- base::transform(AA_freq_plot, mid_y = stats::ave(AA_freq_plot$freqs, AA_freq_plot$position, FUN = function(val) base::cumsum(val) - (0.5 * val)))
      AA_freq_plot$mid_y <- 100-AA_freq_plot$mid_y

      #initialize variables
      position <- amino_acid <- mid_y <- freqs <-  NULL

      if(i==1){
        #plot
        aa_freq_curr_length <- ggplot2::ggplot(AA_freq_plot, ggplot2::aes(x=position,y=freqs,fill=amino_acid,label=amino_acid))+
          ggplot2::geom_bar(stat="identity",show.legend = FALSE)+
          #ggplot2::geom_text(ggplot2::aes(y = mid_y),size=1,colour='black',fontface='bold',show.legend=FALSE)+
          ggplot2::geom_hline(yintercept=c(25,50,75),linetype="dotted")+
          ggplot2::theme_bw()+
          ggplot2::labs(x = "CDR3 sequence position", y = "Occurrence AA [%]", fill = "", colour = "") +
          ggplot2::scale_fill_manual(values = colorset)+
          ggplot2::scale_colour_manual(values = darkened_colors) +
          ggplot2::scale_x_discrete(limit = factor(c(1:number_of_positions)))+
          ggplot2::theme(axis.title.x = ggplot2::element_text(size=20),axis.title.y = ggplot2::element_text(size=20),axis.text.x = ggplot2::element_text(size=15, angle = 0),axis.text.y = ggplot2::element_text(size=20, angle = 0))
          #.theme.akbar()
      }else if(i==length(repertoire)){
        aa_freq_curr_length <- ggplot2::ggplot(AA_freq_plot, ggplot2::aes(x=position,y=freqs,fill=amino_acid,label=amino_acid))+
          ggplot2::geom_bar(stat="identity")+
          #ggplot2::geom_text(ggplot2::aes(y = mid_y),size=1,colour='black',fontface='bold',show.legend=FALSE)+
          ggplot2::geom_hline(yintercept=c(25,50,75),linetype="dotted")+
          ggplot2::theme_bw()+
          ggplot2::labs(x = "", y = "", fill = "", colour = "") +
          ggplot2::scale_fill_manual(values = colorset,  
                                     guide = ggplot2::guide_legend(
                                      direction = "vertical",
                                      label.position = "right",
                                      ncol=1,
                                      override.aes = list(size=7),
                                      label.theme = ggplot2::element_text(angle = 0)))+
          ggplot2::scale_colour_manual(values = darkened_colors) +
          ggplot2::scale_x_discrete(limit = factor(c(1:number_of_positions)))+
          ggplot2::theme(axis.title.x = ggplot2::element_text(size=20),axis.title.y = ggplot2::element_text(size=20),axis.text.x = ggplot2::element_text(size=15, angle = 0),axis.text.y = ggplot2::element_text(size=20, angle = 0))
          #.theme.akbar()+
          ggplot2::theme(legend.position="right")
      }else{
        aa_freq_curr_length <- ggplot2::ggplot(AA_freq_plot, ggplot2::aes(x=position,y=freqs,fill=amino_acid,label=amino_acid))+
          ggplot2::geom_bar(stat="identity",show.legend = FALSE)+
          #ggplot2::geom_text(ggplot2::aes(y = mid_y),size=1,colour='black',fontface='bold',show.legend=FALSE)+
          ggplot2::geom_hline(yintercept=c(25,50,75),linetype="dotted")+
          ggplot2::theme_bw()+
          ggplot2::labs(x = "", y = "", fill = "", colour = "") +
          ggplot2::scale_fill_manual(values = colorset)+
          ggplot2::scale_colour_manual(values = darkened_colors) +
          ggplot2::scale_x_discrete(limit = factor(c(1:number_of_positions)))+
          ggplot2::theme(axis.title.x = ggplot2::element_text(size=20),axis.title.y = ggplot2::element_text(size=20),axis.text.x = ggplot2::element_text(size=15, angle = 0),axis.text.y = ggplot2::element_text(size=20, angle = 0))
          #.theme.akbar()
      }

      list_aa_plots[[i]]<-aa_freq_curr_length

    }
  }
  return(list_aa_plots)
}







