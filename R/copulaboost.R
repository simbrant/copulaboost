

copulaboost <- function(y, x, cov_types, n_models = 2, n_covs=5,
                        learning_rate = .1, eps=5e-2, extra_x=NULL,
                        verbose=F, cont_method="Ingrid_mod", subsample_rows = F,
                        prop_sample=.1, n_sample = NULL) {


  # Make an empty version of the object that will be returned
  model <- lapply(1:(n_models + 1), function(i) vector("list", 2))

  # First, compute an intercept (ml)
  intercept <- stats::nlm(f=function(f0) (log(1- plogis(f0)) +
                                     log(plogis(f0)/(1 -plogis(f0))) *
                                     mean(y))**2,
                   p=0.5)$estimate

  model[[1]][[1]] <- intercept
  model[[1]][[2]] <- NULL

  current_prediction <- rep(intercept, length(y))

  # This piece of code compiles the p-values of the regression coefficients
  # of a univariate linear fit, for each continuous and binary variable,
  # and the minimum p-value of the regression coefficients of categorical
  # variables.
  compute_p_value_weights <- function(y, x, x_type){
    x <- as.data.frame(x)
    vals <- rep(1, ncol(x))
    vals[x_type == "d"] <- apply(x[, x_type=="d"], 2,
                                 function(x) length(unique(x)) - 1)
    colind <- cumsum(vals)
    x[, which(x_type=="d")] <- apply(x[, which(x_type == "d")], 2, as.factor)
    x[, which(x_type=="c")] <- apply(x[, which(x_type == "c")], 2, as.numeric)

    x_ <- model.matrix(~., data=as.data.frame(x))
    lmod <- lm(y~x_-1)
    t_vals <- summary(lmod)$coef[-1, 4]
    c(min(t_vals[1:colind[1]]),
      sapply(2:length(colind),
             function(j) min(t_vals[(colind[j-1] + 1):(colind[j])])))
  }

  # Compute the weights, and draw the covariates to include.
  p_val <- compute_p_value_weights(y - plogis(current_prediction), x, cov_types)
  model[[2]][[2]] <- sample(dim(x)[2], n_covs, prob = 1 - p_val)

  # Fit first non-constant model
  model[[2]][[1]] <- copulareg::copulareg(
    y - plogis(current_prediction),
    x = x[, model[[2]][[2]]],
    var_type_y = "d",
    var_type_x = cov_types[model[[2]][[2]]],
    extra_x = extra_x[, model[[2]][[2]]]
    )

  time <- proc.time()
  current_prediction <- (current_prediction +
                           copulareg::predict.copulareg(model[[2]][[1]]))
  if (verbose){
    cat(paste0("First integration took ", round(proc.time()[3] - time[3], 4),
               " seconds.\n"))
  }
  # Fit the rest of the models
  if (n_models >= 2){

    for (m in 3:(n_models + 1)){

      # Choose covariates to be included in the model
      # random sampling as a placeholder for now
      #model[[m]][[2]] <- sample(dim(x)[2], n_covs)
      p_val <- compute_p_value_weights(y-plogis(current_prediction),
                                       x, cov_types)
      model[[m]][[2]] <- sample(dim(x)[2], n_covs, prob = 1 - p_val)

      # Fit model
      if(subsample_rows){

        if (!is.null(n_sample)){
        sample_rows <- sample(length(y), n_sample)
        } else {
          sample_rows <- sample(length(y), ceiling(length(y)*prop_sample))
        }

        model[[m]][[1]] <- copulareg::copulareg((y - plogis(current_prediction))[sample_rows],
                                                x=x[sample_rows, model[[m]][[2]]], var_type_y = "c",
                                                var_type_x = cov_types[model[[m]][[2]]],
                                                extra_x = rbind(extra_x[, model[[m]][[2]]],
                                                                x[-sample_rows, model[[m]][[2]]]),
                                                extra_y = (y - plogis(current_prediction))[-sample_rows])

      } else {
        model[[m]][[1]] <- copulareg::copulareg(y - plogis(current_prediction),
                                                x=x[, model[[m]][[2]]], var_type_y = "c",
                                                var_type_x = cov_types[model[[m]][[2]]],
                                                extra_x = extra_x[, model[[m]][[2]]])
      }

      time <- proc.time()
      if (subsample_rows){
        current_prediction <- (
          current_prediction +
            copulareg::predict.copulareg(model[[m]][[1]],
                                         newdata = x[, model[[m]][[2]]],
                                         eps=eps,
                                         cont_method=cont_method)
        )
      } else{
        current_prediction <- (
          current_prediction +
            copulareg::predict.copulareg(model[[m]][[1]], eps=eps,
                                         cont_method=cont_method)
        )
      }
      if (verbose){
        cat(paste0("Integration no. ", m-1, " took ", round(proc.time()[3] - time[3], 4), " seconds.\n"))
      }
    }
  }

  res <- list(model=model, learning_rate=learning_rate, cov_types=cov_types, eps=eps)
  class(res) <- "copulaboost"
  res
}


predict.copulaboost <- function(model, new_x=NULL, eps=NULL) {

  if (is.null(eps)){
    eps <- model$eps
  }
  if (!is.null(new_x)){
    prediction <- matrix(0, nrow=nrow(new_x), ncol=length(model$model))
    prediction[, 1] <- model$model[[1]][[1]]
    prediction[, 2:length(model$model)] <- vapply(
      2:length(model$model),
      function(m) model$learning_rate*copulareg::predict.copulareg(
        model$model[[m]][[1]],
        newdata = new_x[, model$model[[m]][[2]]],
        eps=eps
      ),
      FUN.VALUE = prediction[, 1]
    )
  } else {
    prediction <- matrix(0, nrow=model$model[[2]][[1]]$model$nobs,
                         ncol=length(model$model))
    prediction[, 1] <- model$model[[1]][[1]]
    prediction[, 2:length(model$model)] <- vapply(
      2:length(model$model),
      function(m) model$learning_rate*copulareg::predict.copulareg(
        model$model[[m]][[1]],
        eps=eps
      ),
      FUN.VALUE = prediction[, 1]
    )
  }

  for (i in seq(nrow(prediction))){
    prediction[i, is.na(prediction[i, ])] <- 0
  }

  list(individual=prediction, total = t(apply(prediction, 1, cumsum)))

}
