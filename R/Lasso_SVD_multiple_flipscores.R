library(future)
library(future.apply)
library(hdi)

#' @title Multiple SVD-Lasso Selection and Flipscores
#' @description Applica `Lasso_SVD_flipscores` su più variabili in parallelo o sequenziale.
#' @param response_var_name Nome variabile risposta.
#' @param data Dataset.
#' @param vars_to_test Variabili da testare.
#' @param n_comp Componenti SVD-Lasso.
#' @param family Famiglia modello.
#' @param score_type Tipo score.
#' @param seed Seed.
#' @param n_flips Numero flips.
#' @param flips Matrice flips.
#' @param alternative Alternativa.
#' @param parallel Calcolo parallelo.
#' @param n_cores Core da usare.
#' @param max_cores Max core.
#' @param ... Altri parametri.
#' @return Lista risultati multipli con correzione MaxT.
#' @export
Lasso_SVD_multiple_flipscores <- function(response_var_name, 
                                          data,
                                          vars_to_test = NULL,
                                          n_comp = round(nrow(data)/2), # minimo 2 componenti
                                          family,
                                          score_type = "standardized", 
                                          seed = NULL, 
                                          n_flips = 5000,
                                          flips = NULL,
                                          alternative = "two.sided",
                                          parallel = TRUE,
                                          n_cores,
                                          max_cores,
                                          ...
){
  response_var_pos <- which(names(data) == response_var_name)
  
  if(is.null(vars_to_test)) {
    vars_to_test_pos <- setdiff(1:length(names(data)), response_var_pos)
    vars_to_test <- names(data)[vars_to_test_pos]
  }
  
  else {
    if (is.character(vars_to_test)){
      vars_to_test_pos <- which(names(data) %in% vars_to_test)
    }
    else  {
      vars_to_test_pos <- vars_to_test
      vars_to_test <- names(data)[vars_to_test_pos]
    }
  }
  
  if(!isTRUE(parallel)) {
    
    if(length(vars_to_test) == 1) {
      
      results <- Lasso_SVD_flipscores(response_var_name, reference_var = vars_to_test, 
                                  x = data[,-c(response_var_pos,vars_to_test_pos)],
                                  y = data[,response_var_pos], data = data, 
                                  n_comp = n_comp, family = family, score_type = score_type, 
                                  seed = seed, n_flips = n_flips, flips = flips, 
                                  alternative = alternative)
      results <- list(results)
    }
    
    else{
      
      results <- lapply(vars_to_test_pos, function(i) Lasso_SVD_flipscores(response_var_name, reference_var = names(data)[i],
                                                                       x = data[,-c(response_var_pos,i)],
                                                                       y = data[,response_var_pos], data = data,
                                                                       n_comp = n_comp, family = family, 
                                                                       score_type = score_type, seed = seed,
                                                                       n_flips = n_flips, flips = flips, alternative = alternative))
      names(results) <- vars_to_test
      
    }
    
  }
  
  else {
    
    if(is.null(n_cores)) {
      n_cores <- min(future::availableCores() - 1, 8)  # Max 8 core per evitare overhead
    }
    
    future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(future::sequential))
    
    results <- future.apply::future_lapply(vars_to_test_pos, function(i) {
      Lasso_SVD_flipscores(response_var_name, reference_var = names(data)[i],
                       x = data[,-c(response_var_pos,i)],
                       y = data[,response_var_pos], data = data,
                       n_comp = n_comp, family = family,
                       score_type = score_type, seed = seed,
                       n_flips = n_flips, flips = flips,
                       alternative = alternative, ...)
    }, future.seed = TRUE)
    
    names(results) <- vars_to_test 
    
  }
  
  scores <- extract_info(results, "scores", data = data)
  Tspace <- extract_info(results, "Tspace", data = data)
  residuals <- extract_info_2(results, "residuals", data = data)
  fitted.values <- extract_info_2(results, "fitted.values", data = data)
  p.values_raw <- data.frame("p_value" = sapply(results, `[[`, "p.values"))
  row.names(p.values_raw) <- vars_to_test
  p.values_maxT <- data.frame("p_value" = maxT(as.matrix(abs(Tspace))))
  row.names(p.values_maxT) <- vars_to_test
  
  model <- list(
    scores = scores,
    Tspace = Tspace,
    residuals = residuals,
    fitted.values = fitted.values,
    p.values = p.values_raw,
    p.values_maxT = p.values_maxT
  )
  
  return(model)
}
