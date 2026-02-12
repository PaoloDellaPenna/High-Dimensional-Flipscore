library(hdi)
library(stats)

#' @title SVD-Lasso Selection and Single Flipscore
#' @description Combina decomposizione SVD e selezione Lasso sulle componenti principali, seguito dal test flipscores sulla variabile target.
#' @param response_var_name Nome variabile risposta.
#' @param reference_var Nome variabile di test.
#' @param x Matrice predittori.
#' @param y Vettore risposta.
#' @param data Dataframe.
#' @param n_comp Numero componenti.
#' @param family Famiglia modello.
#' @param score_type Tipo di score.
#' @param seed Seed.
#' @param n_flips Numero di flips.
#' @param flips Matrice flips.
#' @param alternative Tipo di test.
#' @param ... Altri parametri.
#' @return Risultati del test con riduzione dimensionale SVD.
#' @export
Lasso_SVD_flipscores <- function(response_var_name,
                                 reference_var,
                                 x, 
                                 y, 
                                 data,
                                 n_comp = round(nrow(data)/2), 
                                 family = "gaussian",
                                 score_type = "standardized",
                                 seed = NULL,
                                 n_flips = 5000,
                                 flips = NULL,
                                 alternative = "two.sided",
                                 ...) {
  
  U <- svd(x)$u
  
  pos_sel_comp <- hdi::lasso.firstq(x = U, y = y, q = n_comp, family = family)
  
  if(!is.character(reference_var)){
    reference_var <- names(data)[reference_var]
  }
  
  model <- formula(paste(response_var_name, "~", paste0("U[,",pos_sel_comp ,"]", collapse = " + "), "+", reference_var))
  
  if(!is.null(flips)){
    n_flips <- nrow(flips)
  }
  
  res <- flipscores_single(formula = model, 
                           data = data,
                           to_be_tested = reference_var,
                           family = family, 
                           score_type = score_type,
                           seed = seed,
                           n_flips = n_flips, 
                           flips = flips,
                           alternative = alternative)
  results <- NULL
  results$scores <- res$scores
  results$Tspace <- res$Tspace
  results$coefficients <- res$coefficients
  results$residuals <- res$residuals
  results$fitted.values <- res$fitted.values
  results$p.values <- res$p.values
  results$formula <- res$formula
  results$offset <- res$offset
  
  return(results)
}
