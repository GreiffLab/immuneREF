#RETURNS Amino Acid frequency plot either for one repertoire
# or for a list of repertoires (need to plot using wrapping)

#plot_aafreq(repertoires_analyzed[[1]],mode="single",length_seq='14')
#plot_aafreq(repertoires_analyzed[1:5],mode="list",length_seq='14')

.plot_maxmin<-function(repertoires,mode="list"){

  ev_plot <- .plot_evenness(repertoires,mode=mode,var="evenness")
  aa_plot <-.plot_aafreq(repertoires,mode=mode,length_seq='14')
  arch_plot <- .plot_arch(repertoires,mode=mode,plot_var="components",log_plot=TRUE)
  vdj_plot <- .plot_vdj_diversity(repertoires,mode)

  plot_out <- list(evenness=ev_plot,aa=aa_plot,arch=arch_plot,vdj=vdj_plot)
  return(plot_out)
}
