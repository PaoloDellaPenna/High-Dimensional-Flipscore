#' @title High-Dimensional Flipscore Function
#' @description Funzione principale (wrapper) per eseguire test High-Dimensional Flipscores con diversi metodi di selezione.
#' @param method Metodo di selezione: "lasso", "svd-lasso", o "stepwise".
#' @param response_var_name Name of the response variable.
#' @param data A data frame containing all variables.
#' @param vars_to_test A vector of variables to test. If NULL, all variables except the response are tested.
#' @param n_comp Number of components to select (default: n/2).
#' @param family Model family (e.g., "gaussian", "binomial").
#' @param score_type Type of score to compute: "standardized", "orthogonalized", etc.
#' @param seed Seed for reproducibility.
#' @param n_flips Number of permutations (flips) to perform.
#' @param flips Optional matrix of precomputed flips.
#' @param alternative Type of hypothesis test: "two.sided", "greater", or "less".
#' @param parallel Logical. Whether to run computations in parallel.
#' @param n_cores Number of cores to use if parallel is TRUE.
#' @param direc Direction of selection for stepwise: "forward", "backward", or "both".
#' @param ... Additional arguments passed to internal functions.
#' @return A list containing the results of the flipscore test, including selected variables (if applicable), scores, T-space, residuals, fitted values, and p-values.
#' @export
hd_flipscore <- function(method = c("lasso", "svd-lasso", "stepwise"),
                         response_var_name,
                         data,
                         vars_to_test = NULL,
                         n_comp = round(nrow(data) / 2),
                         family = "gaussian",
                         score_type = "standardized",
                         seed = NULL,
                         n_flips = 5000,
                         flips = NULL,
                         alternative = "two.sided",
                         parallel = TRUE,
                         n_cores = NULL,
                         direc = "both",
                         ...) {
  
  method <- match.arg(method)
  
  if (!response_var_name %in% names(data)) {
    stop("response_var_name non presente nel dataset.")
  }
  
  if (method == "lasso") {
    result <- Lasso_multiple_flipscores(response_var_name = response_var_name,
                                        data = data,
                                        vars_to_test = vars_to_test,
                                        n_comp = n_comp,
                                        family = family,
                                        score_type = score_type,
                                        seed = seed,
                                        n_flips = n_flips,
                                        flips = flips,
                                        alternative = alternative,
                                        parallel = parallel,
                                        n_cores = n_cores,
                                        ...)
    
  } else if (method == "svd-lasso") {
    result <- Lasso_SVD_multiple_flipscores(response_var_name = response_var_name,
                                            data = data,
                                            vars_to_test = vars_to_test,
                                            n_comp = n_comp,
                                            family = family,
                                            score_type = score_type,
                                            seed = seed,
                                            n_flips = n_flips,
                                            flips = flips,
                                            alternative = alternative,
                                            parallel = parallel,
                                            n_cores = n_cores,
                                            ...)
    
  } else if (method == "stepwise") {
    result <- Stepwise_multiple_flipscores(response_var_name = response_var_name,
                                           data = data,
                                           vars_to_test = vars_to_test,
                                           n_comp = n_comp,
                                           direc = direc,
                                           family = family,
                                           score_type = score_type,
                                           seed = seed,
                                           n_flips = n_flips,
                                           flips = flips,
                                           alternative = alternative,
                                           parallel = parallel,
                                           n_cores = n_cores,
                                           ...)
  } else {
    stop("Method not recognised. Use “lasso”, “svd-lasso” or “stepwise”.")
  }
  
  return(result)
}
