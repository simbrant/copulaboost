copulaboost <- function(y, x, cov_types, n_models=100, n_covs=5,
                        learning_rate=0.33, eps=5e-2, verbose=F,
                        cont_method="Localmedian",
                        family_set=c("gaussian", "clayton", "gumbel"),
                        jitter_sel=T, ml_update=F, ml_sel=F, max_ml_scale=1,
                        keep_sel_struct=T, approx_order=2, parametric_margs=F,
                        parallel = F, par_method_sel = "itau",
                        update_intercept=F, model = NULL) {

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
    model <- lapply(1:(n_models + 1), function(i) vector("list", 2))

    # First, compute an intercept (ml)
    intercept <- stats::nlm(
      f = function(f0) (
        log(1 - plogis(f0)) + log(plogis(f0) / (1 - plogis(f0))) * mean(y)
      ) ** 2, p = 0
    )$estimate

    model[[1]][[1]] <- intercept
    model[[1]][[2]] <- NULL

    current_prediction <- rep(intercept, length(y))
    scaling <- rep(learning_rate, n_models + 1)
    scaling[1] <- 1
    start_ind <- 2
    stop_ind <- n_models + 1

  } else {

    current_prediction <- predict(model, new_x=x)
    scaling <- c(model$scaling, rep(learning_rate, n_models))
    model_old <- model$model
    model <- lapply(seq(n_models + length(model_old)),
                    function(i) vector("list", 2))
    model[[1]][[1]] <- model_old[[1]][[1]]

    for (m in 2:length(model_old)) {
      model[[m]][[1]] <- model_old[[m]][[1]]
      model[[m]][[2]] <- model_old[[m]][[2]]
    }

    start_ind <- length(model_old) + 1
    stop_ind <- n_models + length(model_old)
    rm(model_old)

  }

  if (verbose) {
    pb <- txtProgressBar(min = start_ind - 2, max = stop_ind, style = 3,
                         initial = start_ind - 1)
  }

  ll <- function(eta, y) {
    mean(y * eta - log(1 + exp(eta)))
  }

  if (jitter_sel) {
    x_jitt <- x
    x_jitt[, cov_types == "d"] <- apply(
      as.matrix(x_jitt[, cov_types == "d"]),
      2, function(x) x + rnorm(length(x), 0, 1e-3)
    )
    distr_x_jitt <- .compute_distrbs(x_jitt, rep("c", ncol(x_jitt)),
                                     parametric_margs)
  }

  # Fit the rest of the models
  for (m in start_ind:stop_ind) {
    if (jitter_sel) {
      model[[m]][[2]] <- .select_n_cov(
        current_prediction, y, x_jitt, rep("c", length(cov_types)), n_covs,
        m = approx_order, ystar_cont = if(m == 2) FALSE else TRUE,
        ml_update = ml_sel, max_ml_scale=max_ml_scale,
        par_method = par_method_sel, dx = distr_x_jitt,
        dy = .compute_distrbs(y - plogis(current_prediction),
                              "c", parametric_margs), cl=cl
      )

    } else {
      model[[m]][[2]] <- .select_n_cov(
        current_prediction, y, x, cov_types, n_covs, m = approx_order,
        ystar_cont = if(m == 2) FALSE else TRUE, ml_update = ml_sel,
        max_ml_scale=max_ml_scale, par_method = par_method_sel, dx = distr_x,
        dy = .compute_distrbs(y - plogis(current_prediction),
                              "c", parametric_margs), cl=cl
      )
    }

    model[[m]][[1]] <- copulaboost::copulareg(
      y - plogis(current_prediction), x = x[, model[[m]][[2]]],
      var_type_y = if(m == 2) "d" else "c",
      var_type_x = cov_types[model[[m]][[2]]], family_set = family_set,
      dvine = keep_sel_struct,
      distr_x = .extract_margs(distr_x, model[[m]][[2]]),
      distr_y = .compute_distrbs(y - plogis(current_prediction),
                                 if(m == 2) "d" else "c", parametric_margs)
    )

    curr_incr <- predict(
      model[[m]][[1]], eps = eps, cont_method = cont_method
    )

    scaling[m] <- if (ml_update) {
      optim(
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
      upd <- nlm(
        function(upd_to_ic) - sum(
          dbinom(y, 1, plogis(upd_to_ic + current_prediction), T)
        ),
        p=0
      )$estimate
      model[[1]][[1]] <- model[[1]][[1]] + upd
      current_prediction <- current_prediction + upd
    }

    if (verbose) {
      setTxtProgressBar(pb, value = m)
    }

  }

  res <- list(
    model = model, learning_rate = learning_rate, cov_types = cov_types,
    eps = eps, scaling = scaling
  )

  if (parallel) {
    parallel::stopCluster(cl)
  }
  class(res) <- "copulaboost"
  res
}

predict.copulaboost <- function(model, new_x=NULL, verbose=F, all_parts=F, impute_na=T) {

  n_models <- length(model$model) - 1

  predictions <- matrix(
    0, nrow = if (!is.null(new_x)) {
      nrow(new_x)
    } else {
      length(predict(model$model[[2]][[1]]))
    }, ncol = n_models + 1
  )

  predictions[, 1] <- model$model[[1]][[1]]

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_models, style = 3, initial = 1)
  }

  for (j in seq(n_models)) {
    predictions[, j + 1] <- (
      model$scaling[j + 1] *
        predict(model$model[[j + 1]][[1]], new_x[, model$model[[j + 1]][[2]]],
                eps = model$eps)
    )
    if (any(is.na(predictions[, j + 1])) ) {
      predictions[is.na(predictions[, j + 1]), j + 1] <- median(
        predictions[!is.na(predictions[, j + 1]), j + 1]
      )
    }
    if (verbose) {
      setTxtProgressBar(pb, value = j)
    }
  }

  if (all_parts) {
    t(apply(predictions, 1, cumsum))
  } else {
    apply(predictions, 1, sum)
  }

}
