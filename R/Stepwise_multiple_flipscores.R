library(furrr)
library(future)

#' @title Multiple Stepwise Selection and Flipscores
#' @description Applica `Stepwise_flipscores` su più variabili.
#' @param response_var_name Risposta.
#' @param data Dataset.
#' @param vars_to_test Variabili da testare.
#' @param n_comp Numero componenti (qui usato come max_variables).
#' @param direc Direzione stepwise.
#' @param family Famiglia.
#' @param score_type Score.
#' @param seed Seed.
#' @param n_flips Numero flips.
#' @param flips Flips precalcolati.
#' @param alternative Alternativa.
#' @param parallel Parallelo.
#' @param n_cores Cores.
#' @param max_cores Max cores.
#' @param ... Altri argomenti.
#' @return Risultati multipli con MaxT.
#' @export
Stepwise_multiple_flipscores <- function(response_var_name, 
                                         data,
                                         vars_to_test = NULL,
                                         n_comp = round(nrow(data) / 2), 
                                         direc = "both", 
                                         family,
                                         score_type = "standardized", 
                                         seed = NULL, 
                                         n_flips = 5000,
                                         flips = NULL,
                                         alternative = "two.sided",
                                         parallel = TRUE,
                                         n_cores = NULL,
                                         max_cores = 8,
                                         ...) {
  
  all_vars <- names(data)
  response_var_pos <- which(all_vars == response_var_name)
  
  if (is.null(vars_to_test)) {
    vars_to_test <- all_vars[-response_var_pos]
  } else if (is.numeric(vars_to_test)) {
    vars_to_test <- all_vars[vars_to_test]
  }
  
  y <- data[[response_var_name]]
  
  # Precompute flip matrix once
  if (is.null(flips)) {
    flips <- flipscores:::.make_flips(n = nrow(data), B = n_flips)
  }
  
  # Precompute excluded variables for each var
  excluded_list <- lapply(vars_to_test, function(v) setdiff(all_vars, c(response_var_name, v)))
  
  run_flips <- function(i) {
    reference_var <- vars_to_test[i]
    excluded_vars <- excluded_list[[i]]
    
    Stepwise_flipscores(
      response_var_name = response_var_name,
      reference_var = reference_var,
      x = data[excluded_vars],
      y = y,
      data = data,
      max_variables = n_comp,
      direc = direc,
      family = family,
      score_type = score_type,
      seed = seed,
      n_flips = n_flips,
      flips = flips,
      alternative = alternative,
      ...
    )
  }
  
  if (!parallel) {
    results <- lapply(seq_along(vars_to_test), run_flips)
  } else {
    if (is.null(n_cores)) {
      n_cores <- min(future::availableCores() - 1, max_cores)
    }
    
    future::plan(multisession, workers = n_cores)
    on.exit(future::plan(sequential))
    
    results <- furrr::future_map(seq_along(vars_to_test), run_flips, .options = furrr_options(seed = TRUE))
  }
  
  names(results) <- vars_to_test
  
  # Extraction
  selected_variables <- extract_info(results, "selection")
  scores <- extract_info(results, "scores", data = data)
  Tspace <- extract_info(results, "Tspace", data = data)
  residuals <- extract_info_2(results, "residuals", data = data)
  fitted.values <- extract_info_2(results, "fitted.values", data = data)
  p.values_raw <- data.frame("p_value" = sapply(results, `[[`, "p.values"))
  row.names(p.values_raw) <- vars_to_test
  p.values_maxT <- data.frame("p_value" = maxT(as.matrix(abs(Tspace))))
  row.names(p.values_maxT) <- vars_to_test
  
  model <- list(
    selection = selected_variables,
    scores = scores,
    Tspace = Tspace,
    residuals = residuals,
    fitted.values = fitted.values,
    p.values = p.values_raw,
    p.values_maxT = p.values_maxT
  )
  
  return(model)
}
