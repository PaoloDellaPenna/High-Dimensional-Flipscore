library(flipscores)
library(stats)

#' @title Flipscores Single Variable Test
#' @description Funzione wrapper interna per testare una singola variabile utilizzando la metodologia dei flipscores.
#' @param formula Formula del modello da testare.
#' @param data Dataframe contenente i dati.
#' @param to_be_tested Nome o indice della variabile sotto test H0.
#' @param family Famiglia esponenziale per il GLM (es. "gaussian", "binomial").
#' @param score_type Tipo di score ("standardized", "orthogonalized", ecc.).
#' @param seed Seed opzionale.
#' @param n_flips Numero di permutazioni.
#' @param flips Matrice di flips pre-calcolata (opzionale).
#' @param precompute_flips Logico. Se pre-calcolare i flips internamente.
#' @param id Vettore per osservazioni clusterizzate.
#' @param alternative Tipo di alternativa ("two.sided", "greater", "less").
#' @param ... Altri argomenti passati a metodi interni.
#' @return Un oggetto di classe "flipscores".
#' @export
flipscores_single <- function(formula, 
                              data, 
                              to_be_tested = NULL,
                              family, 
                              score_type = "standardized", 
                              seed = NULL, 
                              n_flips = 5000, 
                              flips = NULL, 
                              precompute_flips = TRUE,
                              id = NULL, 
                              alternative = "two.sided",
                              ...){
  fs_call <- mf <- match.call()
  score_type = match.arg(score_type, c("orthogonalized", "standardized", 
                                       "effective", "basic", "my_lab"))
  if (missing(score_type)) 
    stop("test type is not specified or recognized")
  m <- match(c("score_type", "n_flips", "alternative", "id", 
               "seed", "flips", "precompute_flips"), names(mf), 0L)
  m <- m[m > 0]
  flip_param_call = mf[c(1L, m)]
  flip_param_call[[1L]] = flipscores:::.flip_test
  flip_param_call$id = eval(flip_param_call$id, parent.frame())
  flip_param_call$alternative = eval(flip_param_call$alternative, 
                                     parent.frame())
  flip_param_call$flips <- eval(flip_param_call$flips, parent.frame())
  if (!is.null(flip_param_call$flips)) 
    flip_param_call$n_flips = nrow(flip_param_call$flips)
  else {
    flip_param_call$n_flips <- eval(flip_param_call$n_flips, 
                                    parent.frame())
    if (is.null(flip_param_call$n_flips)) 
      flip_param_call$n_flips = 5000
  }
  if (is.null(flip_param_call$precompute_flips)) 
    flip_param_call$precompute_flips = TRUE
  flip_param_call$family <- eval(flip_param_call$family, parent.frame())
  family <- flip_param_call$family
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())()
  }
  flip_param_call$score_type <- eval(flip_param_call$score_type, 
                                     parent.frame())
  if (is.null(flip_param_call$score_type)) 
    flip_param_call$score_type = "standardized"
  flip_param_call$seed <- eval(flip_param_call$seed, parent.frame())
  flip_param_call$nobservations = eval(mf$nobservations, parent.frame())
  mf$nobservations = NULL
  if (!is.null(list(...)$parms_DV)) {
    flip_param_call$parms_DV = list(...)$parms_DV
    mf$parms_DV = NULL
  }
  m2 <- match(c("to_be_tested"), names(mf), 0L)
  if (m2 == 0) 
    to_be_tested = NULL
  else {
    m <- c(m, m2)
    to_be_tested = mf[[m2]]
  }
  if (!is.null(flip_param_call$id) && (score_type == "orthogonalized")) {
    print(warning("WARNING: Use of id is not possible with score_type=='orthogonalized', yet.\n Nothing done."))
    return(NULL)
  }
  if (length(m) > 0) 
    mf <- mf[-m]
  mf$offset = eval(mf$offset, parent.frame())
  model = eval(mf$formula, parent.frame())
  if ("formula" %in% is(model)) {
    mf[[1L]] <- if (!is.null(family) && identical(family$family, "negbinom")) quote(glm.nb) else quote(glm)
    
    mf$family <- family  # <-- QUESTA RIGA È FONDAMENTALE
    param_x_ORIGINAL <- mf$x
    mf$x <- TRUE
    model <- eval(mf, parent.frame())
  }
  else {
    param_x_ORIGINAL <- TRUE
    model <- update(model, x = TRUE)
  }
  if (is.null(model$y)) 
    model$y = model$model[, 1]
  to_be_tested = eval(to_be_tested, parent.frame())
  if (is.null(to_be_tested)) 
    to_be_tested = colnames(model[["x"]])
  else {
    if (is.numeric(to_be_tested)) 
      to_be_tested = colnames(model[["x"]])[to_be_tested]
    to_be_tested = eval(to_be_tested, parent.frame())
  }
  if (!is.null(flip_param_call$flips)) {
    flip_param_call$precompute_flips = FALSE
    flip_param_call$n_flips = nrow(eval(flip_param_call$flips, 
                                        parent.frame()))
  }
  else if (flip_param_call$precompute_flips) {
    set.seed(seed)
    flip_param_call$flips = flipscores:::.make_flips(max(nrow(model$model), 
                                                         flip_param_call$nobservations), flip_param_call$n_flips, 
                                                     flip_param_call$id)
  }
  results = flipscores:::socket_compute_scores_and_flip(to_be_tested, model,flip_param_call = flip_param_call)
  model$scores = data.frame(results[[1]]$scores)
  # model$scores = data.frame(lapply(results, function(x) x[[1]]$scores))
  nrm = attributes(results[[1]]$scores)$scale_objects$nrm
  # nrm = sapply(results, function(x) attributes(x[[1]]$scores)$scale_objects$nrm)
  std_dev = attributes(results[[1]]$scores)$sd
  # std_dev = sapply(results, function(x) attributes(x[[1]]$scores)$sd)
  model$Tspace = data.frame(results[[1]]$Tspace)
  # model$Tspace = data.frame(lapply(results, function(x) x[[1]]$Tspace))
  model$p.values = results[[1]]$p.values
  # model$p.values = sapply(results, function(x) x[[1]]$p.values)
  flip_param_call$flips = NULL
  attr(model$scores, "nrm") = nrm
  attr(model$scores, "sd") = std_dev
  attr(model$scores, "resid_std") = data.frame(attr(results[[1]]$scores, "resid_std"))
  # attr(model$scores, "resid_std") = data.frame(lapply(results, 
  #                                                     function(x) attr(x[[1]]$scores, "resid_std")))
  names(attributes(model$scores)$resid_std) <- names(nrm) <- names(std_dev) <- names(model$scores) <- names(model$Tspace) <- names(model$p.values) <- to_be_tested
  model$call = fs_call
  model$flip_param_call = flip_param_call
  model$score_type = score_type
  if (is.null(param_x_ORIGINAL) || (!param_x_ORIGINAL)) 
    model$x = NULL
  class(model) <- c("flipscores", class(model))
  return(model)
}
