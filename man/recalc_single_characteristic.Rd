% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recalc_single_characteristic.R
\name{recalc_single_characteristic}
\alias{recalc_single_characteristic}
\title{Function to recalculate single immuneREF layer}
\usage{
recalc_single_characteristic(
  analyzed_list,
  characteristic,
  repertoire_df,
  models,
  species,
  receptor,
  chain,
  identifier_rep,
  vdj_list,
  dictionary,
  architecture_depth = 1,
  aa_range = c(8:20)
)
}
\arguments{
\item{analyzed_list}{prepared and named repertoires_analyzed list}

\item{characteristic}{layer of interest}

\item{repertoire_df}{repertoire that needs to recalculated}

\item{models}{immunosignature models}

\item{species}{species of repertoire}

\item{receptor}{receptor of repertoire}

\item{chain}{chain of repertoire}

\item{identifier_rep}{name of repertoire}

\item{vdj_list}{reference germline gene list}

\item{dictionary}{dictionary for kmer calculation}

\item{architecture_depth}{maximal LD for network calculation (a network is calculated for each distance 1:architecture_depth)}

\item{aa_range}{range of sequence lengths for which amino acid frequencies should be evaluated}
}
\value{
List of repertoire layers
}
\description{
Function to recalculate single immuneREF layer
}
\examples{
\dontrun{
recalc_single_characteristic(analyzed_list,characteristic,
repertoire_df, svm, "mm", "igh", "repertoire_1",
list_germline_genes_allele_01, kmer_dict)
}
}
