#' List of gapped kmer patterns.
#'
#' A list that contains templates for the counting of gapped-kmers
#' with parameters k = 1, m <= 1 for amino acid
#' and parameters k = 3, m <= 3 for nucleotide sequences.
#'
#' @format A list of with 4 entries :
#' \describe{
#'   \item{nt}{List of pairs of nucleotide 3-mers}
#'   \item{aa}{List of pairs of amino acid  1-mers}
#'   \item{gap_nt}{Numeric vector containing allowed gap sizes for nucleotide sequences}
#'   \item{gap_aa}{Numeric vector containing allowed gap sizes for amino acid sequences}
#' }
"dictionary_counts"


#' Collection of germline genes and frequencies
#'
#' A list containing sublists for species ("hs","mm") which in turn
#' contain sublists for receptors ("ig","tr") which are subset in
#' chains ("h", "k", "l" and "b", "a", respectively). Each entry
#' contains a list of three dataframes ("V","D" and "J") with the major IMGT
#' annotated germline genes including name, sequence based on IMGT and frequencies based on
#' experimental data from DeWitt(2017), Emerson (2017), Greiff (2017) and Madi (2017).
#'
#' @format A list of lists containing dataframes with up to 126 entries:
#' \describe{
#'   \item{gene}{name of germline gene}
#'   \item{allele}{allele number (presently restricted to allele 01)}
#'   \item{sequence}{nucleotide sequence of germline gene}
#'   \item{species}{name of species}
#'   \item{frequency}{Frequencies of germline genes based on experimental data}
#' }
#' @source {
#'    \url{http://www.imgt.org/vquest/refseqh.html}
#'    \url{https://doi.org/10.1371/journal.pone.0160853}
#'    \url{https://doi.org/10.1038/ng.3822}
#'    \url{https://doi.org/10.1016/j.celrep.2017.04.054}
#'    \url{https://doi.org/10.7554/eLife.22057}
#'    }
"list_germline_genes"




#' A pretrained Support Vector Machine model for the prediction of public clones
#'
#' An SVM model for the prediction of public clones trained on murine B-cell data.
#' As published by Greiff et al. (2017)
#'
#' @format A kebabs SVM model.
#' \describe{
#'   \item{svm_model_public}{Trained on 109'568 nucleotide sequences with gapped-kmers k = 3,m = 1}
#' }
#' @source \url{https://doi.org/10.4049/jimmunol.1700594}
"svm_model_public"


#' Tutorial repertoires for code examples
#'
#' A list containing five simulated repertoires using the immuneSIM package.
#'
#' @format A list containing 4 dataframes containing 10'000 sequences each.
#' \describe{
#'   \item{Repertoire 1}{sim_mm_ig_h_2_0__A}
#'   \item{Repertoire 2}{sim_mm_ig_h_4_0__A}
#'   \item{Repertoire 3}{sim_hs_ig_h_2_0__A}
#'   \item{Repertoire 4}{sim_hs_ig_h_4_0__A}
#' }
"tutorial_repertoires"
