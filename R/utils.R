library(flipscores)

#' @title Utility: Extract Info from List
#' @description Funzione ausiliaria per estrarre oggetti specifici (es. p.values) dalle liste di risultati.
#' @param my.list Lista di risultati.
#' @param object Nome oggetto da estrarre.
#' @param data Dataframe riferimento.
#' @export
extract_info <- function(my.list, object, data) {
  
  extract <- sapply(my.list, `[[`, object)
  
  if(is.list(extract)){
    
    out <- as.data.frame(do.call(cbind, lapply(extract, `length<-`, max(lengths(extract)))))
    
  }
  
  else {
    
    out <- extract
    
  }
  
  colnames(out) <- names(my.list)
  
  return(out)
  
}

#' @title Utility: Extract Info 2
#' @description Variante semplificata per estrazione diretta.
#' @param my.list Lista risultati.
#' @param object Oggetto da estrarre.
#' @param data Dataframe.
#' @export
extract_info_2 <- function(my.list, object, data) {
  
  out <- sapply(my.list, `[[`, object)
  
  colnames(out) <- names(my.list)
  
  return(out)
  
}

#' @title Max-T Correction
#' @description Calcola i p-value corretti per molteplicità usando la statistica Max-T basata su permutazioni.
#' @param permTabs Matrice delle statistiche test permutate.
#' @export
maxT <- function(permTabs){
  mxt_flips=apply(permTabs[-1,],1,max)
  sapply(permTabs[1,],function(tob)flipscores:::.t2p(c(tob,mxt_flips)))
}
