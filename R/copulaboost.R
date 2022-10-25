copulaboost <- function(y, x, cov_types, n_models=100, n_covs=5,
                        learning_rate=0.33, eps=5e-2, verbose=FALSE,
                        cont_method="Localmedian",
                        family_set=c("gaussian", "clayton", "gumbel"),
                        jitter_sel=TRUE, ml_update=FALSE, ml_sel=FALSE,
                        max_ml_scale=1, keep_sel_struct=TRUE, approx_order=2,
                        parametric_margs=TRUE, parallel = FALSE,
                        par_method_sel = "itau", update_intercept=TRUE,
                        model = NULL, xtreme = FALSE) {
  ##' copulaboost
  ##' @name copulaboost
  ##' @aliases copulaboost
  ##' @description This is the main function of the package, which
  ##' fits an additive model with a fixed number of components, each
  ##' involving a fixed number of covariates, where each component is a 
  ##' copula regression model.
  ##' @param y A vector of n observations of the (univariate) binary outcome
  ##' variable y
  ##' @param x A (n x p) matrix of n observations of p covariates
  ##' @param cov_types A vector of p characters that have to take the value
  ##' "c" or "d" to indicate whether each margin of the covariates is discrete
  ##' or continuous.
  ##' @param n_models The number of model components to fit.
  ##' @param n_covs The number of covariates included in each component.
  ##' @param learning_rate Factor to scale (down) the each component.
  ##' @param eps Control parameter for the approximation to the conditional
  ##' expectation (the prediction) for each copula model (component), which
  ##' splits the interval [-1, 1] into equal pieces of eps length. 
  ##' @param verbose Logical indicator of whether a progressbar should be shown
  ##' in the terminal.
  ##' @param cont_method Method to use for the approximation of each conditional
  ##' expectation, can either be "Localmedian" or "Trapezoidalsurv", for the 
  ##' former, see section 3.2 of 
  ##' https://arxiv.org/ftp/arxiv/papers/2208/2208.04669.pdf. The latter uses 
  ##' the so called "Darth vader rule" in conjuction with a simple translative
  ##' transformation to write the conditional expectation as an integral along
  ##' the conditional survival function, which is then approximated by the 
  ##' trapezoidal method.
  ##' @param family_set A vector of strings that specifies the set of
  ##' pair-copula families that the fitting algorithm chooses from. For an
  ##' overview of which values that can be specified, see the documentation for
  ##' \link[rvinecopulib]{bicop}.
  ##' @param jitter_sel Logical indicator of whether jittering should
  ##' be used for any discrete covariates when selecting the variables for each
  ##' component (improves computational speed).
  ##' @param ml_update Logical indicator of whether each new component should
  ##' be scaled by a number between 0 and max_ml_scale by maximising the 
  ##' log-likelihood of the scaling factor given the current model and the new
  ##' component.
  ##' @param ml_sel The same as ml_update, but for the variable selection
  ##' algorithm.
  ##' @param max_ml_scale The maximum scaling factor allowed for each component.
  ##' @param keep_sel_struct Logical indicator of whether the d-vine structures found
  ##' by the model selection algorithm should be kept when fitting the
  ##' components.
  ##' @param approx_order The order of the approximation used for evaluating the
  ##' conditional expectations when selecting covariates for each component. The
  ##' allowed values for approx_order are 1, 2, 3, 4, 5, and 6.
  ##' @param parametric_margs Logical indicator of whether parametric (gaussian
  ##' or bernoulli) models should be used for the marginal distributions of the
  ##' covariates. 
  ##' @param parallel (Experimental) Logical indicator of whether
  ##' parallelization should be used when selecting covariates.
  ##' @param par_method_sel Estimation method for copulas used when selecting
  ##' the model components, either "itau" or "mle", see the documentation for
  ##' \link[rvinecopulib]{bicop}.
  ##' @param update_intercept Logical indicator of whether the intercept
  ##' parameter should be updated (by univariate maximum likelihood) after
  ##' each component is added. 
  ##' @param model Initial copulaboost-model. If model is a copulaboost model
  ##' with k components, the resulting model will have k + n_models components.
  ##' @param xtreme (Experimental) Logical indicator of whether a second order 
  ##' expansion of the log-likelihood should be used in each gradient boosting 
  ##' step, similar to the xgboost algorithm.
  ##' @return A copulaboost object, which contains a nested list 'object$model'
  ##' which contains all of the model components. The first element of each
  ##' list contains a copulareg object, and the second element contains a vector
  ##' listing the indexes of the covariates that are a part of the component.
  ##' The object also contains a list of the updated intercepts 
  ##' 'object$f0_updated' at each stage of the fitting process, so that the
  ##' j-th intercept is the intercept for the model that is the weighted sum of 
  ##' the j first components. 'object$scaling' contains a vector of weights for
  ##' each components, equal to the learning rate, possibly multiplied by an
  ##' individual factor if ml_update = TRUE. In addition the object contains
  ##' the values of the arguments learning_rate, cov_types, and eps that where
  ##' used when calling copulaboost().
  ##' @examples 
  ##' # Compile some test data
  ##' data('ChickWeight')
  ##' set.seed(10)
  ##' tr <- sample(c(TRUE, FALSE), nrow(ChickWeight), TRUE, c(0.7, 0.3))
  ##' y_tr <- as.numeric(ChickWeight$weight[tr] > 100)
  ##' y_te <- as.numeric(ChickWeight$weight[!tr] > 100)
  ##' x_tr <- apply(ChickWeight[tr, -1], 2, as.numeric)
  ##' x_te <- apply(ChickWeight[!tr, -1], 2, as.numeric)
  ##' cov_types <- apply(x_tr, 2,
  ##'                    function(x) if(length(unique(x)) < 10) "d" else "c")
  ##'
  ##' # Fit model to training data
  ##' md <- copulaboost::copulaboost(y_tr, x_tr, cov_types, n_covs = 2, 
  ##'                                n_models = 5, verbose = TRUE)
  ##'
  ##' # Out of sample predictions for a new data matrix
  ##' preds <- predict(md, new_x = x_te, all_parts = TRUE)
  ##'
  ##' # Plot log-likelihood
  ##' plot(apply(preds, 2,
  ##'            function(eta) {
  ##'              sum(stats::dbinom(y_te, 1, stats::plogis(eta), log = TRUE))
  ##'              }),
  ##'      type = "s")
  ##'
  ##' @export

  if (ncol(x) != length(cov_types)) {
    stop(
        "length(cov_types) must be equal to ncol(x)!"
      )
  }

  if (parallel) {
    cl <- parallel::makeCluster(parallel::detectCores())
  } else {
    cl <- NULL
  }

  # Compute marginal distributions
  distr_x <- .compute_distrbs(x, cov_types, parametric = parametric_margs)

  if (is.null(model)) {
    # Make an empty version of the object that will be returned
    model <- lapply(1:(n_models), function(i) vector("list", 2))
    f_0 <- 0
    f_0_updated <- rep(0, n_models)

    # First, compute an intercept (ml)
    f_0 <- stats::nlm(
      f = function(f0) (
        log(1 - stats::plogis(f0)) + 
          log(stats::plogis(f0) / (1 - stats::plogis(f0))) * mean(y)
      ) ** 2, p = 0
    )$estimate

    if (!update_intercept) {
      f_0_updated <- rep(f_0, n_models)
    }
    current_prediction <- rep(f_0, length(y))
    scaling <- rep(learning_rate, n_models)
    start_ind <- 1
    stop_ind <- n_models

  } else {

    # FIX
    current_prediction <- predict.copulareg(model, new_x=x)
    scaling <- c(model$scaling, rep(learning_rate, n_models))
    model_old <- model$model
    f_0_updated <- c(model$f_0_updated, rep(0, n_models))
    f_0 <- model$f_0
    if (!update_intercept) {
      f_0_updated[seq(length(f_0_updated))] <- f_0
    }
    model <- lapply(seq(n_models + length(model_old)),
                    function(i) vector("list", 2))
    for (m in 1:length(model_old)) {
      model[[m]][[1]] <- model_old[[m]][[1]]
      model[[m]][[2]] <- model_old[[m]][[2]]
    }

    start_ind <- length(model_old) + 1
    stop_ind <- n_models + length(model_old)
    rm(model_old)
  }

  if (verbose) {
    pb <- utils::txtProgressBar(min = start_ind - 2, max = stop_ind, style = 3,
                                initial = start_ind - 1)
  }

  ll <- function(eta, y) {
    mean(y * eta - log(1 + exp(eta)))
  }

  if (jitter_sel) {
    x_jitt <- x
    x_jitt[, cov_types == "d"] <- apply(
      as.matrix(x_jitt[, cov_types == "d"]),
      2, function(x) x + stats::rnorm(length(x), 0, 1e-3)
    )
    distr_x_jitt <- .compute_distrbs(x_jitt, rep("c", ncol(x_jitt)),
                                     parametric_margs)
  }

  # Fit the rest of the models
  for (m in start_ind:stop_ind) {
    if (jitter_sel) {
      model[[m]][[2]] <- .select_n_cov(
        current_prediction, y, x_jitt, rep("c", length(cov_types)), n_covs,
        m = approx_order, ystar_cont = if(m == 1) FALSE else TRUE,
        ml_update = ml_sel, max_ml_scale=max_ml_scale,
        par_method = par_method_sel, dx = distr_x_jitt,
        dy = .compute_distrbs(
          if(!xtreme) {
            y - stats::plogis(current_prediction)
          } else {
            (y - stats::plogis(current_prediction)) / (
              stats::plogis(current_prediction) * 
                (1 - stats::plogis(current_prediction))
            )
          }, if(m == 1) "d" else "c",
          if (TRUE) FALSE else parametric_margs
        ), cl=cl, xtreme=xtreme
      )

    } else {
      model[[m]][[2]] <- .select_n_cov(
        current_prediction, y, x, cov_types, n_covs, m = approx_order,
        ystar_cont = if(m == 1) FALSE else TRUE, ml_update = ml_sel,
        max_ml_scale=max_ml_scale, par_method = par_method_sel, dx = distr_x,
        dy =.compute_distrbs(
          if(!xtreme) {
            y - stats::plogis(current_prediction)
          } else {
            (y - stats::plogis(current_prediction)) / (
              stats::plogis(current_prediction) * 
                (1 - stats::plogis(current_prediction))
            )
          }, if(m == 1) "d" else "c",
          if (TRUE) FALSE else parametric_margs
        ), cl=cl, xtreme=xtreme
      )
    }

    model[[m]][[1]] <- copulaboost::copulareg(
      if(!xtreme) {
        y - stats::plogis(current_prediction)
        } else {
          (y - stats::plogis(current_prediction)) / (
            stats::plogis(current_prediction) * 
              (1 - stats::plogis(current_prediction))
          )
        }, x = x[, model[[m]][[2]]],
      var_type_y = if(m == 1) "d" else "c",
      var_type_x = cov_types[model[[m]][[2]]], family_set = family_set,
      dvine = keep_sel_struct,
      distr_x = .extract_margs(distr_x, model[[m]][[2]]),
      distr_y = .compute_distrbs(
        if(!xtreme) {
          y - stats::plogis(current_prediction)
        } else {
          (y - stats::plogis(current_prediction)) / 
            (stats::plogis(current_prediction) * 
               (1 - stats::plogis(current_prediction))
          )
        }, if(m == 1) "d" else "c",
        if (TRUE) FALSE else parametric_margs
      )
    )

    curr_incr <- predict.copulareg(
      model[[m]][[1]], eps = eps, cont_method = cont_method
    )

    scaling[m] <- if (ml_update) {
      stats::optim(
        fn = function(theta) {
          - ll(current_prediction + theta * curr_incr, y = y)
        }, par = 1, method = "Brent", lower = 0, upper = max_ml_scale
      )$par*learning_rate
    } else {
      learning_rate
    }

    current_prediction <- (
      current_prediction + scaling[m] * curr_incr
    )

    if (update_intercept) {
      upd <- stats::nlm(
        function(upd_to_ic) - sum(
          stats::dbinom(y, 1, stats::plogis(upd_to_ic + current_prediction), T)
        ),
        p=0
      )$estimate
      f_0_updated[m] <- if(m == 1) f_0 + upd else f_0_updated[m-1] + upd
      current_prediction <- current_prediction + upd
    }

    if (verbose) {
      utils::setTxtProgressBar(pb, value = m)
    }

  }

  res <- list(
    model = model, learning_rate = learning_rate, cov_types = cov_types,
    eps = eps, scaling = scaling, f_0_updated = f_0_updated,
    f_0 = f_0)

  if (parallel) {
    parallel::stopCluster(cl)
  }

  class(res) <- "copulaboost"
  res
}

predict.copulaboost <- function(object, new_x=NULL, verbose=FALSE,
                                all_parts=FALSE, impute_na=TRUE, ...) {
  ##' predict.copulaboost
  ##' @name predict.copulaboost
  ##' @aliases predict.copulaboost
  ##' @description Function that computes predictions for a copulaboost model.
  ##' @param object A copulaboost model object returned from the copulaboost 
  ##' function.
  ##' @param new_x A matrix of covariate values to compute predictions for. If
  ##' new_x is not provided, the function will return predictions for the
  ##' data used for fitting the model.
  ##' @param verbose A logical indicator of whether a progressbar should be
  ##' shown in the terminal.
  ##' @param all_parts A logical indicator of whether predictions for the k
  ##' first components for k = 1, ..., n_models should be returned as a matrix.
  ##' If all_parts = FALSE, only the prediction with all of the model components
  ##' is returned.
  ##' @param impute_na A logical indicator of whether any potential NA values 
  ##' from any of the predictions from each component should be replaced by the
  ##' median prediction for that component.
  ##' @param ... further arguments passed to or from other methods.
  ##' @return If all_parts=FALSE this function will return a list of predictions
  ##' for each row of new_x (if specified, otherwise predictions for the 
  ##' training data will be returned). If all_parts=TRUE a matrix will be
  ##' returned, where the j-th column contains predictions using the j first
  ##' model components.

  n_models <- length(object$model)

  predictions <- matrix(
    0, nrow = if (!is.null(new_x)) {
      nrow(new_x)
    } else {
      length(predict.copulareg(object$model[[2]][[1]]))
    }, ncol = n_models
  )

  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = n_models, style = 3, initial = 1)
  }

  for (j in seq(n_models)) {
    predictions[, j] <- (
      object$scaling[j] *
        predict.copulareg(object$model[[j]][[1]],
                          new_x[, object$model[[j]][[2]]],
                          eps = object$eps)
    )
    if (any(is.na(predictions[, j])) ) {
      predictions[is.na(predictions[, j]), j] <- stats::median(
        predictions[!is.na(predictions[, j]), j]
      )
    }
    if (verbose) {
      utils::setTxtProgressBar(pb, value = j)
    }
  }

  if (all_parts) {
    all_preds <- t(apply(predictions, 1, cumsum))
    for (j in 1:n_models){
      all_preds[, j] <- all_preds[, j] + object$f_0_updated[j]
    }
    all_preds
  } else {
    apply(predictions, 1, sum) + object$f_0_updated[n_models]
  }
}
