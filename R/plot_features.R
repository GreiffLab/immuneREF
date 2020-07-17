#plot all plottable layers

.plot_features <- function(repertoire,mode="single",selected_features=""){

  ev_plot <- .plot_evenness(repertoire,mode="single",var="evenness")
  aa_plot <- .plot_aafreq(repertoire,mode="single",length_seq='14')
  arch_plot <- .plot_arch(repertoire,mode="single",plot_var="components",log_plot=TRUE)
  vdj_plot <- .plot_vdj_diversity(repertoire,"single")

}
