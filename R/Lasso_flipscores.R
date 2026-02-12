library(hdi)
library(stats)

#' @title Lasso Selection and Single Flipscore
#' @description Esegue una selezione delle variabili di disturbo tramite Lasso (pacchetto hdi) e successivamente testa una singola variabile target usando i flipscores.
#' @param response_var_name Nome della variabile risposta (Y).
#' @param reference_var Nome o indice della variabile di test.
#' @param x Matrice dei predittori.
#' @param y Vettore della risposta.
#' @param data Dataframe completo.
#' @param n_comp Numero di componenti da selezionare col Lasso (default n/2).
#' @param family Famiglia del modello (es. "gaussian").
#' @param score_type Tipo di score per flipscores.
#' @param seed Seed per riproducibilità.
#' @param n_flips Numero di permutazioni.
#' @param flips Matrice flips pre-calcolata.
#' @param alternative Alternativa del test.
#' @param ... Altri argomenti.
#' @return Risultati del test inclusa la selezione Lasso e p-value.
#' @export
Lasso_flipscores <- function(response_var_name, 
                             reference_var, # nel dataframe
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
  
    pos_sel_comp <- hdi::lasso.firstq(x = x, y = y, q = n_comp, family = family)
    
    selected_variables <- names(x)[pos_sel_comp]
    
    if(!is.character(reference_var)){
      reference_var <- names(data)[reference_var]
    }
    
    model <- formula(paste(response_var_name, "~", 
                           paste(colnames(x[,pos_sel_comp]), 
                                 collapse = " + "), "+", reference_var))
    
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
    results$selection <- selected_variables
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
