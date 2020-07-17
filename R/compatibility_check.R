#' Tests immuneREF-compatibility of repertoire
#'
#' @param repertoire Input immune repertoire
#' @param species Species of input repertoire
#' @param receptor Receptor chain of immune repertoire ("igh" or "trb")
#' @return Outputs ERRORs when encountered and returns number of errors.
#' @examples
#' \dontrun{
#' compatibility_check(repertoire, "mm", "igh")
#' }

#tests whether a repertoire is compatible with immuneREF or is likely gonna throw an error

compatibility_check <- function(repertoire,species,receptor){

  error_sum<-0

  #create tables of V, D, and J gene and evaluate whether the genes are in the list of germline genes used.
  repertoire_v_genes <- as.character(repertoire$v_call)
  is_in_table_V <- table(repertoire_v_genes %in% list_germline_genes[[species]][[substr(receptor,1,2)]][[substr(receptor,3,3)]][["V"]][["gene"]])

  repertoire_d_genes <- as.character(repertoire$d_call)
  is_in_table_D <- table(repertoire_d_genes %in% list_germline_genes[[species]][[substr(receptor,1,2)]][[substr(receptor,3,3)]][["D"]][["gene"]])

  repertoire_j_genes <- as.character(repertoire$j_call)
  is_in_table_J <- table(repertoire_j_genes %in% list_germline_genes[[species]][[substr(receptor,1,2)]][[substr(receptor,3,3)]][["J"]][["gene"]])

  #if some genes are missing report to user and add error
  if("FALSE" %in% names(is_in_table_V)){
    cat("ERROR: Check V names against germline_gene list\n")
    error_sum<-error_sum+1
  }

  if("FALSE" %in% names(is_in_table_D)){
    cat("ERROR: Check D names against germline_gene list\n")
    error_sum<-error_sum+1
  }

  if("FALSE" %in% names(is_in_table_J)){
    cat("ERROR: Check J names against germline_gene list\n")
    error_sum<-error_sum+1
  }


  #check all needed columns are present and with correct name
  if(!(table(c("sequence", "junction", "junction", "junction_aa","freqs","v_call","d_call","j_call") %in% names(repertoire))[[TRUE]] == 8)){

    cat("ERROR: Need to rename/add columns\n")
    error_sum<-error_sum+1

  }

  #check clone count is a freq that sums up to one.
  if(sum(repertoire$freqs) != 1){

    if(round(sum(repertoire$freqs),3) == 1){
      cat("NOTE: Clonal frequency rounding inaccuracy\n")
    }else{
      cat("ERROR: Check clonal frequency (doesnt sum up to 1)\n")
      error_sum<-error_sum+1

    }

  }

  return(error_sum)

}

