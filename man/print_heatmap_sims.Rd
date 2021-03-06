% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_heatmap_sims.R
\name{print_heatmap_sims}
\alias{print_heatmap_sims}
\title{Prints heatmap of all immuneREF layers in input list}
\usage{
print_heatmap_sims(
  list_similarity_matrices,
  annotation_list,
  path_figure = "figures"
)
}
\arguments{
\item{list_similarity_matrices}{list of immuneREF layers}

\item{annotation_list}{List with entries 'categories' and 'colors' for heatmap annotation}

\item{path_figure}{Path where plots should be stored}
}
\value{
TRUE
}
\description{
Prints heatmap of all immuneREF layers in input list
}
\examples{
\dontrun{
print_heatmap_sims(list_similarity_matrices,annotation_list,path_figure="figures")
}
}
