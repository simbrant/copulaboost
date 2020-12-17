

copulaboost <- function(y, x, cov_types, n_models = 100, n_covs=5,
                        learning_rate = 0.33, eps=5e-2, extra_x=NULL,
                        verbose=F, cont_method="Ingrid", subsample_rows = F,
                        prop_sample=.1, n_sample = NULL, ml_update=T,
                        sel_method="p-value", lambda_ridgesel = 1) {

  if (sel_method == "p-value-ridge"){
    invmat <- .comp_invmat_ridge(x, cov_types, lambda = lambda_ridgesel)
  } else {
    invmat <- NULL
  }

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
  scaling <- rep(1, n_models + 1)

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_models + 1, style = 3, initial = 1)
  }

  ll <- function(eta, y){
    mean(y*eta - log(1 + exp(eta)))
  }

  # Compute the weights, and draw the covariates to include.
  model[[2]][[2]] <- .sel_covs(y - plogis(current_prediction), x, cov_types,
                               n_covs, method=sel_method, invmat=invmat,
                               lambda=lambda_ridgesel)

  # Fit first non-constant model
  model[[2]][[1]] <- copulareg::copulareg(
    y - plogis(current_prediction),
    x = x[, model[[2]][[2]]],
    var_type_y = "d",
    var_type_x = cov_types[model[[2]][[2]]],
    extra_x = extra_x[, model[[2]][[2]]]
    )


  curr_incr <- copulareg::predict.copulareg(model[[2]][[1]])
  if (ml_update){
    scaling[2] <- optim(par=1,
                        fn=function(theta) -ll(current_prediction +
                                                 theta*curr_incr, y=y),
                        method = "Brent", lower = 0, upper = 100
                        )$par
  }

  current_prediction <- current_prediction + learning_rate*scaling[2]*curr_incr

  if (verbose) {
    setTxtProgressBar(pb, value = 2)
  }

  # Fit the rest of the models
  if (n_models >= 2){

    for (m in 3:(n_models + 1)){

      model[[m]][[2]] <- .sel_covs(y - plogis(current_prediction), x, cov_types,
                                   n_covs, method=sel_method, invmat=invmat,
                                   lambda=lambda_ridgesel)
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
        model[[m]][[1]] <- copulareg::copulareg(
          y - plogis(current_prediction), x=x[, model[[m]][[2]]],
          var_type_y = "c", var_type_x = cov_types[model[[m]][[2]]],
          extra_x = extra_x[, model[[m]][[2]]]
          )
      }

      if (subsample_rows){
        curr_incr <- copulareg::predict.copulareg(
          model[[m]][[1]],
          newdata = x[, model[[m]][[2]]],
          eps=eps,
          cont_method=cont_method
        )
      } else{
        curr_incr <- copulareg::predict.copulareg(
          model[[m]][[1]], eps=eps, cont_method=cont_method
          )
      }
      if (ml_update) {
        scaling[m] <- optim(par=1,
                            fn=function(theta) -ll(current_prediction +
                                                     theta*curr_incr, y=y),
                            method = "Brent", lower = 0, upper = 100
        )$par

      }
      current_prediction <- (current_prediction +
                               learning_rate*scaling[m]*curr_incr)
      if (verbose) {
        setTxtProgressBar(pb, value = m)
      }
    }
  }

  res <- list(model=model, learning_rate=learning_rate, cov_types=cov_types,
              eps=eps, scaling=scaling)
  class(res) <- "copulaboost"
  res
}


predict.copulaboost <- function(model, new_x=NULL, eps=NULL,
                                cont_method="Ingrid", verbose = F) {

  if (is.null(eps)){
    eps <- model$eps
  }

  if (verbose){
    pb <- txtProgressBar(min=0, max = length(model$model), style = 3)
  } else{
    pb <- NULL
  }

  pred_m <- function(m, learning_rate, pb = NULL, model, new_x=NULL, eps){
    # Function that implements the code that is looped
    # using vapply
    if (!is.null(pb)){
      setTxtProgressBar(pb, value=m)
    }
    if (!is.null(new_x)){
      learning_rate[m]*copulareg::predict.copulareg(
        model[[m]][[1]], newdata = new_x[, model[[m]][[2]]], eps=eps,
        cont_method=cont_method
      )
    } else {
      learning_rate[m]*copulareg::predict.copulareg(model[[m]][[1]], eps=eps,
                                                 cont_method=cont_method)
    }
  }

  nrowoutput <- switch(1*(!is.null(new_x)) + 1,
                       model$model[[2]][[1]]$model$nobs,
                       nrow(new_x))

  prediction <- matrix(0, nrow=nrowoutput, ncol=length(model$model))
  prediction[, 1] <- model$model[[1]][[1]]

  if (!is.null(new_x)){
    prediction[, 2:length(model$model)] <- vapply(
      2:length(model$model), pred_m,
      learning_rate=model$learning_rate*model$scaling, pb=pb,
      model=model$model, new_x=new_x, eps=eps, FUN.VALUE = prediction[, 1]
    )
  } else {
    prediction[, 2:length(model$model)] <- vapply(
      2:length(model$model), pred_m,
      learning_rate=model$learning_rate*model$scaling, pb=pb,
      model=model$model, eps=eps, FUN.VALUE = prediction[, 1]
    )
  }

  for (i in seq(nrow(prediction))){
    prediction[i, is.na(prediction[i, ])] <- 0
  }

  list(individual=prediction, total = t(apply(prediction, 1, cumsum)))

}

# Covariate selection methods
.sel_covs <- function(y, x, x_type, n_covs, method="p-value", invmat=NULL,
                      lambda = 1){

  if (method == "p-value"){
    .p_val_sel(y, x, x_type, n_covs)
  } else if (method == "p-value-ridge") {
    .ridge_p_val_sel(y, x, x_type, n_covs, invmat = invmat, lambda = lambda)
  } else if (method == "coef-ridge") {
    .ridge_coef_sel(y, x, x_type, n_covs, lambda)
  } else if (method == "cor-greedy") {
    .greedy_cor_sel(y, x, n_covs)
  } else if (method == "random") {
    sample(ncol(x), n_covs)
  } else {
    stop(paste0("Method '", method, "' not implemented!"))
  }

}

.p_val_sel <- function(y, x, x_type, n_covs) {

  # This piece of code compiles the p-values of the regression coefficients
  # of a univariate linear fit, for each continuous and binary variable,
  # and the minimum p-value of the regression coefficients of categorical
  # variables.

  compute_p_value_weights <- function(y, x, x_type){

    x <- as.data.frame(x)
    vals <- rep(1, ncol(x))
    vals[x_type == "d"] <- apply(as.matrix(x[, x_type == "d"]), 2,
                                 function(x) length(unique(x)) - 1)
    colind <- cumsum(vals)
    x[, which(x_type=="d")] <- apply(as.matrix(x[, which(x_type == "d")]),
                                     2, as.factor)
    x[, which(x_type=="c")] <- apply(as.matrix(x[, which(x_type == "c")]),
                                     2, as.numeric)

    x_ <- model.matrix(~., data=as.data.frame(x))
    lmod <- lm(y~x_-1)
    t_vals <- summary(lmod)$coef[-1, 4]

    c(min(t_vals[1:colind[1]]),
      sapply(2:length(colind),
             function(j) min(t_vals[(colind[j-1] + 1):(colind[j])])))

  }

  p_val <- compute_p_value_weights(y, x, x_type)

  sample(dim(x)[2], n_covs, prob = 1 - p_val)

}

.comp_invmat_ridge <- function(x, x_type, lambda=1){

  x <- as.data.frame(x)
  vals <- rep(1, ncol(x))
  vals[x_type == "d"] <- apply(as.matrix(x[, x_type=="d"]), 2,
                               function(x) length(unique(x)) - 1)
  colind <- cumsum(vals)
  x[, which(x_type=="d")] <- apply(as.matrix(x[, which(x_type == "d")]), 2, as.factor)
  x[, which(x_type=="c")] <- apply(as.matrix(x[, which(x_type == "c")]), 2, as.numeric)

  X <- model.matrix(~.-1, data=as.data.frame(x))

  # First compute inverse of (X^tX + lambdaI), which should be provided, really.
  invmat <- solve(t(X)%*%X + diag(rep(lambda, ncol(X))))

}

.ridge_p_val_sel <- function(y, x, x_type, n_covs, lambda = 1, invmat = NULL){

  x <- as.data.frame(x)
  vals <- rep(1, ncol(x))
  vals[x_type == "d"] <- apply(as.matrix(x[, x_type=="d"]), 2,
                               function(x) length(unique(x)) - 1)
  colind <- cumsum(vals)
  x[, which(x_type=="d")] <- apply(as.matrix(x[, which(x_type == "d")]), 2, as.factor)
  x[, which(x_type=="c")] <- apply(as.matrix(x[, which(x_type == "c")]), 2, as.numeric)

  X <- model.matrix(~.-1, data=as.data.frame(x))

  # First compute inverse of (X^tX + lambdaI), which should be provided, really.
  if(is.null(invmat)){
    invmat <- .comp_invmat_ridge(x, x_type, lambda=lambda)
  }
  hat <- X %*% invmat %*% t(X)
  nu <- nrow(X) - sum(diag(2*hat - hat %*% t(hat)))
  betahat <- invmat %*% t(X) %*% (y - mean(y))
  sigmahat <- as.numeric((t((y -mean(y)) - X %*% betahat) %*%
                            ((y - mean(y)) - X %*% betahat))/nu)
  var_betahat <- diag(sigmahat * (invmat %*% (t(X) %*% X) %*% invmat))

  vals <- 2*pt(-abs(betahat / sqrt(var_betahat)), df=nu)

  p_fin <- c(min(vals[1:colind[1]]),
             sapply(2:length(colind),
                    function(j) min(vals[(colind[j-1] + 1):(colind[j])])))


  sample(dim(x)[2], n_covs, prob = 1 - p_fin)
}

.ridge_coef_sel <- function(y, x, x_type, n_covs, lambda = 1) {

  x <- as.data.frame(x)
  vals <- rep(1, ncol(x))
  vals[x_type == "d"] <- apply(as.matrix(x[, x_type=="d"]), 2,
                               function(x) length(unique(x)) - 1)
  colind <- cumsum(vals)
  x[, which(x_type=="d")] <- apply(as.matrix(x[, which(x_type == "d")]), 2, as.factor)
  x[, which(x_type=="c")] <- apply(as.matrix(x[, which(x_type == "c")]), 2, as.numeric)

  X <- model.matrix(~.-1, data=as.data.frame(x))
  betas <- as.numeric(glmnet::glmnet(X, y, alpha=0, lambda=lambda)$beta)[-1]

  b_fin <- c(max(abs(betas[1:colind[1]])),
             sapply(2:length(colind),
                    function(j) max(abs(betas[(colind[j-1] + 1):(colind[j])]))))
  sample(dim(x)[2], n_covs, prob = b_fin)
}

.greedy_cor_sel <- function(y, x, n_covs, start="random") {

  # Initital covariate
  if (start=="random"){
    k <- sample(ncol(x), 1)
  } else if (start == "maxcor"){
    k <- which.max(abs(cor(x, y, method="kendall")))
  } else {
    stop(paste0("Start: ", start, " not implemented."))
  }

  for (j in 2:n_covs){
    add_2_p1 <- sapply((1:ncol(x))[-k], function(i) abs(cor(x[, i], y,
                                                            method = "kendall")))
    add_2_p2 <- sapply((1:ncol(x))[-k],
                       function(i) sum(abs(cor(x[, k], x[, i],
                                               method = "kendall"))))
    cov2add <- which.max((add_2_p1 + obj_p1) / (obj_p2 + add_2_p2))
    k <- c(k, (1:ncol(x))[-k][cov2add])
    obj_p1 <- obj_p1 + add_2_p1[k]
    obj_p2 <- obj_p2 + add_2_p2[k]
  }

  return(k)

}
