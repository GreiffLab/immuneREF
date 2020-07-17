#' Outputs repertoire comparison plot
#'
#' @param list_repertoires List of AIRR-compliant repertoires.
#' @param mode Input format "list" or "single" repertoire
#' @param aa_freq_length Length for which AA-frequency plot should be output
#' @param chain Chain of receptor (Default="H", alternatives not implemented)
#' @param path_figure Path where figures should be saved
#' @param name_plots Naming for plot
#' @return TRUE
#' @examples
#' \dontrun{
#' print_repertoire_comparison(list_repertoires)
#' }

print_repertoire_comparison <- function(list_repertoires,mode="list",aa_freq_length=14,chain="H",path_figure="",name_plots=""){
	ev_plot <- .plot_evenness(repertoire=list_repertoires,mode=mode,var="evenness")
	aa_plot <-.plot_aafreq(repertoire=list_repertoires,mode=mode,length_seq=as.character(aa_freq_length))
	vdj_plots <- .plot_vdj_diversity(repertoire=list_repertoires,chain=chain,mode=mode)

	grDevices::pdf(file.path(path_figure,paste("all_",name_plots,".pdf",sep="")), width = 15,  height = 20,onefile=FALSE)
	grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,26), height=1,width=1)) # 3 rows, 1 columns
	vplayout<-function(x,y) grid::viewport(layout.pos.row=x,layout.pos.col=y)
	base::print(aa_plot[[1]], vp=vplayout(1,1:12))
	base::print(aa_plot[[2]], vp=vplayout(1,13:26))
	base::print(vdj_plots[["V"]], vp=vplayout(2,1:26))
	base::print(vdj_plots[["J"]], vp=vplayout(3,1:13))
	base::print(ev_plot, vp=vplayout(3,14:26))
	grDevices::dev.off()


	return(TRUE)
}
