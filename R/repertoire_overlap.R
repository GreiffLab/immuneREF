#' CALCULATES PAIRWISE OVERLAP BETWEEN REPERTOIRES
#'
#' @param list_simulated_repertoires list of repertoires to be analyzed
#' @param basis Basis on which overlap is calculated (default="CDR3_aa)
#' @return overlap_stat
#' @examples
#' \dontrun{
#' repertoire_overlap(list_simulated_repertoires=tutorial_repertoires,basis="CDR3_aa")
#' }

repertoire_overlap <- function(list_simulated_repertoires, basis="CDR3_aa"){
  #calculates the overlap between every pair of repertoires based on either CDR3 aa,nt or VDJ aa,nt. Default (CDR3 a.a.)
  #output: matrix of every pairwise comparison
  intersect_percentage_function <- function(x,y){
    c(length(intersect(x,y)))/mean(c(length(x), length(y)))
  }


  #save repnames and make combos
  rep_names <- names(list_simulated_repertoires)
  name_combo <- combn(rep_names, 2)

  #calculate overlaps based on chosen sequence
  #overlap is calculated as intersect(x,y)/min(length(x),length(y))
  if(basis=="CDR3_aa"){
    overlap_tmp <- sapply(split(name_combo, col(name_combo)),
                          function(x) intersect_percentage_function(unique(as.character(list_simulated_repertoires[[x[1]]]$junction_aa)),unique(as.character(list_simulated_repertoires[[x[2]]]$junction_aa))))
  }else if(basis == "CDR3_nt"){
    overlap_tmp <- sapply(split(name_combo, col(name_combo)),
                          function(x) intersect_percentage_function(unique(as.character(list_simulated_repertoires[[x[1]]]$junction)),unique(as.character(list_simulated_repertoires[[x[2]]]$junction))))
  }else if(basis == "VDJ_aa"){
    overlap_tmp <- sapply(split(name_combo, col(name_combo)),
                          function(x) intersect_percentage_function(unique(as.character(list_simulated_repertoires[[x[1]]]$sequence_aa)),unique(as.character(list_simulated_repertoires[[x[2]]]$sequence_aa))))
  }else if(basis == "VDJ_nt"){
    overlap_tmp <- sapply(split(name_combo, col(name_combo)),
                          function(x) intersect_percentage_function(unique(as.character(list_simulated_repertoires[[x[1]]]$sequence)),unique(as.character(list_simulated_repertoires[[x[2]]]$sequence))))
  }else if(basis == "V_CDR3_J_nt"){
    overlap_tmp <- sapply(split(name_combo, col(name_combo)),
                          function(x) intersect_percentage_function(unique(paste(list_repertoires[[x[1]]]$v_call,as.character(list_repertoires[[x[1]]]$junction),list_repertoires[[x[1]]]$j_call,sep="_")),unique(paste(list_repertoires[[x[2]]]$v_call,as.character(list_repertoires[[x[2]]]$junction),list_repertoires[[x[2]]]$j_call,sep="_"))))
  }else if(basis == "V_CDR3_J_aa"){
    overlap_tmp <- sapply(split(name_combo, col(name_combo)),
                          function(x) intersect_percentage_function(unique(paste(list_repertoires[[x[1]]]$v_call,as.character(list_repertoires[[x[1]]]$junction_aa),list_repertoires[[x[1]]]$j_call,sep="_")),unique(paste(list_repertoires[[x[2]]]$v_call,as.character(list_repertoires[[x[2]]]$junction_aa),list_repertoires[[x[2]]]$j_call,sep="_"))))
  }

  ## initialize matrix (initialize to 1 as diagonal is 100%overlap)
  overlap_stat <- matrix(1, ncol = length(list_simulated_repertoires), nrow = length(list_simulated_repertoires))
  ## fill in lower triangle
  overlap_stat[lower.tri(overlap_stat)] <- overlap_tmp
  ## transpose and fill in lower triangle of transposed
  toverlap_stat <- t(overlap_stat)
  toverlap_stat[lower.tri(toverlap_stat)] <- overlap_tmp
  ## transpose back to original order
  overlap_stat <- t(toverlap_stat)
  ## add col and rownames
  colnames(overlap_stat) <- rownames(overlap_stat) <- rep_names

  return(overlap_stat)
}
