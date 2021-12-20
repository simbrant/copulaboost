copulaboost <- function(y, x, cov_types, n_models=100, n_covs=5,
                        learning_rate=0.33, eps=5e-2, verbose=F,
                        cont_method="Localmedian",
                        family_set=c("gaussian", "clayton", "gumbel"),
                        jitter_sel=T, ml_update=F, ml_sel=F, max_ml_scale=1,
                        keep_sel_struct=T, approx_order=2, parametric_margs=T,
                        parallel = F, par_method_sel = "itau",
                        update_intercept=T, model = NULL,
                        xtreme = F) {

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
        log(1 - plogis(f0)) + log(plogis(f0) / (1 - plogis(f0))) * mean(y)
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
    current_prediction <- predict(model, new_x=x)
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
        m = approx_order, ystar_cont = if(m == 1) FALSE else TRUE,
        ml_update = ml_sel, max_ml_scale=max_ml_scale,
        par_method = par_method_sel, dx = distr_x_jitt,
        dy = .compute_distrbs(
          if(!xtreme) {
            y - plogis(current_prediction)
          } else {
            (y - plogis(current_prediction)) / (
              plogis(current_prediction)*(1 - plogis(current_prediction))
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
            y - plogis(current_prediction)
          } else {
            (y - plogis(current_prediction)) / (
              plogis(current_prediction)*(1 - plogis(current_prediction))
            )
          }, if(m == 1) "d" else "c",
          if (TRUE) FALSE else parametric_margs
        ), cl=cl, xtreme=xtreme
      )
    }

    model[[m]][[1]] <- copulaboost::copulareg(
      if(!xtreme) {
        y - plogis(current_prediction)
        } else {
          (y - plogis(current_prediction)) / (
            plogis(current_prediction) * (1 - plogis(current_prediction))
          )
        }, x = x[, model[[m]][[2]]],
      var_type_y = if(m == 1) "d" else "c",
      var_type_x = cov_types[model[[m]][[2]]], family_set = family_set,
      dvine = keep_sel_struct,
      distr_x = .extract_margs(distr_x, model[[m]][[2]]),
      distr_y = .compute_distrbs(
        if(!xtreme) {
          y - plogis(current_prediction)
        } else {
          (y - plogis(current_prediction)) / (
            plogis(current_prediction)*(1 - plogis(current_prediction))
          )
        }, if(m == 1) "d" else "c",
        if (TRUE) FALSE else parametric_margs
      )
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
      f_0_updated[m] <- if(m == 1) f_0 + upd else f_0_updated[m-1] + upd
      current_prediction <- current_prediction + upd
    }

    if (verbose) {
      setTxtProgressBar(pb, value = m)
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

predict.copulaboost <- function(model, new_x=NULL, verbose=F, all_parts=F, impute_na=T) {

  n_models <- length(model$model)

  predictions <- matrix(
    0, nrow = if (!is.null(new_x)) {
      nrow(new_x)
    } else {
      length(predict(model$model[[2]][[1]]))
    }, ncol = n_models
  )

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_models, style = 3, initial = 1)
  }

  for (j in seq(n_models)) {
    predictions[, j] <- (
      model$scaling[j] *
        predict(model$model[[j]][[1]], new_x[, model$model[[j]][[2]]],
                eps = model$eps)
    )
    if (any(is.na(predictions[, j])) ) {
      predictions[is.na(predictions[, j]), j] <- median(
        predictions[!is.na(predictions[, j]), j]
      )
    }
    if (verbose) {
      setTxtProgressBar(pb, value = j)
    }
  }

  if (all_parts) {
    all_preds <- t(apply(predictions, 1, cumsum))
    for (j in 1:n_models){
      all_preds[, j] <- all_preds[, j] + model$f_0_updated[j]
    }
    all_preds
  } else {
    apply(predictions, 1, sum) + model$f_0_updated[n_models]
  }
}
