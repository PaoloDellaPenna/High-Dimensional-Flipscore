library(stats)

#' @title Stepwise Selection and Single Flipscore
#' @description Esegue una selezione Stepwise (forward/backward) delle variabili di disturbo e successivamente testa la variabile target.
#' @param response_var_name Variabile risposta.
#' @param reference_var Variabile di test.
#' @param x Predittori.
#' @param y Risposta.
#' @param data Dataset.
#' @param direc Direzione stepwise ("both", "forward", "backward").
#' @param max_variables Numero massimo variabili selezionabili.
#' @param family Famiglia modello.
#' @param score_type Tipo score.
#' @param seed Seed.
#' @param n_flips Numero flips.
#' @param flips Matrice flips.
#' @param alternative Alternativa.
#' @param ... Altri parametri.
#' @return Risultati selezione stepwise e test.
#' @export
Stepwise_flipscores <- function(response_var_name, 
                                reference_var, 
                                x, 
                                y, 
                                data, 
                                direc,
                                max_variables = round(nrow(data)/2), 
                                family = "gaussian",
                                score_type = "standardized", 
                                seed = NULL, 
                                n_flips = 5000,
                                flips = NULL,
                                alternative = "two.sided",
                                ...) {
  
  # Combine response and predictors for stepwise
  full_data <- data
  full_data[[response_var_name]] <- y
  
  # Construct initial and full models
  predictors <- names(x)
  reference_var <- if (!is.character(reference_var)) names(data)[reference_var] else reference_var
  predictors <- setdiff(predictors, reference_var)
  
  null_model <- as.formula(paste(response_var_name, "~", 1))
  full_model <- as.formula(paste(response_var_name, "~", 
                                 paste(predictors, collapse = " + ")))
  
  # Fit initial model
  base_fit <- glm(null_model, data = full_data, family = family)
  
  # Stepwise selection (forward)
  step_fit <- step(base_fit, 
                   scope = list(lower = null_model, upper = full_model), 
                   direction = direc, 
                   trace = 0)
  
  # Limit number of selected variables (excluding the reference_var)
  selected_vars <- setdiff(attr(terms(step_fit), "term.labels"), reference_var)
  if (length(selected_vars) > max_variables) {
    selected_vars <- selected_vars[1:max_variables]
  }
  
  final_model_formula <- as.formula(paste(response_var_name, "~", 
                                          paste(c(selected_vars, reference_var), collapse = " + ")))
  
  # Adjust number of flips if flips are provided
  if (!is.null(flips)) {
    n_flips <- nrow(flips)
  }
  
  # Run flip test
  res <- flipscores_single(formula = final_model_formula, 
                           data = data,
                           to_be_tested = reference_var,
                           family = family, 
                           score_type = score_type, 
                           seed = seed,
                           n_flips = n_flips, 
                           flips = flips,
                           alternative = alternative,
                           ...)
  
  # Collect results
    results <- list(
    selection = selected_vars,
    scores = res$scores,
    Tspace = res$Tspace,
    coefficients = res$coefficients,
    residuals = res$residuals,
    fitted.values = res$fitted.values,
    p.values = res$p.values,
    formula = res$formula,
    offset = res$offset
  )
  
  return(results)
}
