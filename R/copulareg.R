copulareg <- function(y, x, var_type_y, var_type_x, distr_x=NULL, distr_y=NULL,
                      dvine=FALSE,
                      family_set=c("gaussian", "clayton", "gumbel")) {
  ##' copulareg
  ##' @aliases copulareg
  ##' @description This function fits joint distributions with an R-vine
  ##' pair copula structure, that is constructed in a specific way so that
  ##' the conditional density and distribution of the variable y can be
  ##' computed explicitly.
  ##' @param y A vector of n observations of the (univariate)
  ##' outcome variable y
  ##' @param x A (n x p) matrix of n observations of p covariates
  ##' @param var_type_y A character that has to be specified as "d" or "c"
  ##' to indicate whether y is discrete or continuous, respectively.
  ##' @param var_type_x A vector of p characters that have to take the value
  ##' "c" or "d" to indicate whether each margin of the covariates is discrete
  ##' or continuous.
  ##' @param distr_x Internally created object that contains a nested list.
  ##' The first element of the 'outer' list is a list where each element
  ##' is a distribution object similar to that created by ecdf, the second
  ##' element is a function 'transform' that takes a matrix of values of x, 
  ##' and returns the corresponding cumulative distributions F(x).
  ##' @param distr_y Similar to distr_x, but for the outcome y.
  ##' @param dvine Logical variable indicating whether the function should fit
  ##' a canonical d-vine to (y, x) with the ordering (x_1, ..., x_p, y).
  ##' @param family_set A vector of strings that specifies the set of
  ##' pair-copula families that the fitting algorithm chooses from. For an
  ##' overview of which values that can be specified, see the documentation
  ##' for \link[rvinecopulib]{bicop}.
  ##' @return A copulareg object. Consists of 'model', an rvinecopulib
  ##' 'vinecop' object for the copula-model for (x_1, ..., x_p, y),
  ##' a hash table containing all of the conditionals for the model (for the
  ##' training set), objects distr_x and distr_y that contain the marginal
  ##' distributions of the covariates and the outcome, and y, the y-values for
  ##' the training data.

  fit <- .fit_model_xy(y, x, var_type_y, var_type_x, family_set, dvine, distr_x,
                       distr_y)
  class(fit) <- "copulareg"

  fit
}

.fit_model_xy <- function(y, x, var_type_y, var_types_x, family_set, dvine,
                          distr_x, distr_y) {

  # Fits a pair copula model, by first fitting a model to the p-dimensional X,
  # and then augmenting this with the variable y, unless dvine=TRUE or
  # ncol(x) = 1. If dvine = TRUE, the function will fit a canonical d-vine to
  # (y, x) with the ordering (x_1, ..., x_p, y). If ncol(x) = 1 a bivariate
  # copula will be fitted to (x, y).

  if (is.null(distr_x)) {
    distr_x <- .compute_distrbs(x, var_types_x)
  }

  if (is.null(distr_y)) {
    distr_y <- .compute_distrbs(y, var_type_y)
  }

  if (ncol(as.matrix(x)) > 1) {

    if (dvine) {

      u_x <- distr_x$transform(x, var_types_x)
      u_y <- distr_y$transform(y, var_type_y)

      u_all <- cbind(u_x[, seq_len(ncol(x))], u_y[, 1])

      if (any(var_types_x == "d")) {
        u_all <- cbind(u_all, u_x[, seq(from = ncol(x) + 1, to = ncol(u_x))])
      }

      if (var_type_y == "d") {
        u_all <- cbind(u_all, u_y[, 2])
      }

      model_xy <- rvinecopulib::vinecop(
        u_all, var_types = c(var_types_x, var_type_y), family_set = family_set,
        structure = rvinecopulib::dvine_structure(
          c(ncol(x) + 1, seq_len(ncol(x)))
        )
      )

      list(
        model = model_xy,
        conditionals = .compute_conditionals(u_all, model_xy),
        distr_x = distr_x, distr_y = distr_y, y = y
      )

    } else {
      model_x <- .fit_model_x(x, distr_x, var_types_x, family_set = family_set)
      .fit_model_y(y, distr_y, var_type_y, family_set, x, distr_x, model_x)
    }

  } else {

    u <- .weave_transformed(distr_x$transform(x, var_types_x),
                           distr_y$transform(y, var_type_y))

    model_xy <- rvinecopulib::vinecop(
      u, c(var_types_x, var_type_y), family_set = family_set,
      structure = rvinecopulib::dvine_structure(c(2, 1))
    )

    list(
      model = model_xy,
      conditionals = .compute_conditionals(u, model_xy),
      distr_x = distr_x, distr_y = distr_y, y = y
    )

  }
}


.fit_model_x <- function(x, distr_x, var_types_x, family_set) {
  # This function fits a pair copula model with an rvine structure
  # to the data x.
  u_x <- distr_x$transform(x, var_types_x)
  rvinecopulib::vinecop(data = u_x, var_types = var_types_x,
                        family_set = family_set)
}

.fit_model_y <- function(y, distr_y, var_type_y, family_set, x, distr_x,
                         model_x) {

  # Extract number of covariates, p
  p <- model_x$structure$d

  # Augment the matrix with a new zero row and column so that its dimension is
  # increased by 1 in both directions. The new variable (no. 4 in this case)
  # must be positioned in the bottom left corner, so this zero is replaced by
  # the new variable.
  #
  # Example:
  #    Original model matrix:
  #                          2 2 2
  #                          3 3 0
  #                          1 0 0
  #    New model matrix :
  #                          0 2 2 2
  #                          0 3 3 0
  #                          0 1 0 0
  #                          4 0 0 0

  new_mat <- cbind(c(rep(0, p), p + 1),
                   rbind(rvinecopulib::as_rvine_matrix(model_x$structure),
                         rep(0, p)))

  # Compute transformed distributions
  u_x <- distr_x$transform(x, model_x$var_types)

  # make a vinecopula object to store the (x, y) model in,
  # initially a copy of the model of x.
  model_xy <- model_x

  # Compute all the conditionals from the x - model
  conditionals <- .compute_conditionals(u_x, model_x)

  # Insert the pseudo-observations for y
  assign(.make_hfunc_key(p + 1, c()), distr_y$transform(y, var_type_y),
         envir = conditionals)

  for (t in 1:p) {

    # Compute which potential edges satisfy the proximity condition
    valid <- .get_valid_edges(t, new_mat)

    # Find the conditioning set, empty for t=1
    cond_set <- switch(1 + (t > 1), c(), new_mat[1:(t - 1), 1])

    if (length(valid) == 1) {
      # Only one possible edge
      new_mat[t, 1] <- valid
    } else {

      # More than one possibility, have to compute conditional Kendall's
      # coefficients, and select the variable that has the largest absolute
      # conditional Kendall's tau with the response.
      u_x <- sapply(valid, function(j) .get_hfunc(j, cond_set, conditionals))
      u_y <- .get_hfunc(p + 1, cond_set, conditionals)
      new_mat[t, 1] <- valid[which.max(
        abs(stats::cor(u_y, u_x, method = "kendall"))
        )]
    }

    # Combine the two variables that the new pair copula will be fit to.
    u_i_t <- .weave_transformed(get(.make_hfunc_key(p + 1, cond_set),
                                    envir = conditionals),
                                get(.make_hfunc_key(new_mat[t, 1],
                                                    cond_set),
                                    envir = conditionals))

    # Fit the new pair-copula
    if (t > length(model_xy$pair_copulas)) {
      model_xy$pair_copulas <- append(
        model_xy$pair_copulas,
        list(list(rvinecopulib::bicop(
          data = u_i_t,
          var_types = c(var_type_y,
                        model_x$var_types[new_mat[t, 1]]),
          family_set = family_set,
          par_method = "mle",
          selcrit = "aic"))), after = t)

    } else {
      model_xy$pair_copulas[[t]] <- append(
        model_xy$pair_copulas[[t]],
        list(rvinecopulib::bicop(
          data = u_i_t,
          var_types = c(var_type_y,
                        model_x$var_types[new_mat[t, 1]]),
          family_set = family_set,
          par_method = "mle",
          selcrit = "aic")), after = 0)
    }

    # update the log-likelihood
    model_xy$loglik <- model_xy$loglik + model_xy$pair_copulas[[t]][[1]]$loglik

    # Compute the new conditionals
    assign(.make_hfunc_key(p + 1, new_mat[1:t, 1]),
           .bicop_2_hbicop(u_i_t, bicop_obj = model_xy$pair_copulas[[t]][[1]],
                           cond_var = 2, return_u_minus = var_type_y == "d"),
           envir = conditionals)

  }

  # Complete the model object by setting the structure, the variable types,
  # and the total number of parameters
  model_xy$structure <- rvinecopulib::as_rvine_structure(new_mat)
  model_xy$var_types <- c(model_xy$var_types, var_type_y)
  model_xy$npars <- sum(
    sapply(1:(model_xy$structure$d - 1),
           function(t) sum(
             sapply(1:(model_xy$structure$d - t),
                    function(j) model_xy$pair_copulas[[t]][[j]]$npars))))

  list(
    model = model_xy, conditionals = conditionals, distr_x = distr_x,
    distr_y = distr_y, y = y
  )

}

predict.copulareg <- function(object, new_x=NULL, eps=1E-2,
                              cont_method="Localmedian", ...) {
  ##' predict.copulareg
  ##' @aliases predict.copulareg
  ##' @description Computes predictions based on a fitted copulareg model.
  ##' @param object Model fit as returned by copulareg
  ##' @param new_x optional matrix of covariate values to compute the predicted
  ##' values of the outcome for. If not specified, the predicted values for the
  ##' training sample is returned.
  ##' @param eps Interval between each interpolation point when integrating to
  ##' evaluate the predicted value of y, in the case where y is continuous. If
  ##' y is discrete this parameter is ignored.
  ##' @param cont_method Specifies the method used to compute the expected
  ##' values. Can be specified as 'Localmedian' or 'Trapezoidalsurv'. The
  ##' first method divides the range of the observed values of y into
  ##' subintervals according to the argument 'eps', where the sub-integral
  ##' is approximated as the measure of the interval weighted by the local
  ##' median on the interval. The second method computes the integral by
  ##' integrating the survival function using the trapezoidal rule, by
  ##' transforming the outcome into a positive variable by adding a constant.
  ##' @param ... further arguments passed to or from other methods.
  ##' @return A list of predicted values for each row of new_x, if specified,
  ##' otherwise, predictions for each row of the training data is returned.

  eval_at_u_y <- function(u_y, model, n_test, cond_x, y_type) {

    vinemat <- rvinecopulib::as_rvine_matrix(model$structure)
    h <- rep(0, n_test)
    u_ <- .weave_transformed(
      u_y, get(.make_hfunc_key(vinemat[1, 1], c()), envir = cond_x)
    )
    h <- .bicop_2_hbicop(u_, model$pair_copulas[[1]][[1]],
                         return_u_minus = (y_type == "d"))

    if (nrow(vinemat) > 2) {
      for (t in 2:(nrow(vinemat) - 1)) {
        u_ <- .weave_transformed(
          h, get(.make_hfunc_key(vinemat[t, 1], vinemat[seq_len(t - 1), 1]),
                 envir = cond_x)
        )
        h <- .bicop_2_hbicop(u_, model$pair_copulas[[t]][[1]],
                             return_u_minus = (y_type == "d"))
      }
    }
    h
  }

  # Extract x and y types from model
  y_type <- object$model$var_types[object$model$structure$d]
  x_type <- object$model$var_types[1:(object$model$structure$d - 1)]

  # Transform covariates
  if (!is.null(new_x)) {

    n_obs <- nrow(as.matrix(new_x))
    u_x <- object$distr_x$transform(new_x, x_type)

    if (object$model$structure$d > 2) {
      # Extract covariate model from model
      sub_mod <- object$model
      sub_mod$structure <- rvinecopulib::as_rvine_structure(
        rvinecopulib::as_rvine_matrix(
          sub_mod$structure
        )[seq_len(sub_mod$structure$d - 1),
          2:sub_mod$structure$d]
      )
      sub_mod$pair_copulas <- lapply(
        seq_len(sub_mod$structure$d - 1),
        function(l) {
          sub_mod$pair_copulas[[l]][2:(sub_mod$structure$d - (l - 1))]
        }
      )
      sub_mod$var_types <- x_type

      # Compute all the conditionals of the covariate model for the
      # data points where we want to compute the expectation
      cond_x <- .compute_conditionals(u_x, sub_mod)

    } else {
      cond_x <- new.env(hash = TRUE)
      assign(.make_hfunc_key(1, c()),
             u_x, envir = cond_x)
    }

  } else {
    cond_x <- object$conditionals
    n_obs <- object$model$nobs
  }

  if (y_type == "c") {

    ###################
    # Continuous case #
    ###################

    if (!(cont_method %in% c("Localmedian", "Trapezoidalsurv"))) {
      stop(paste0(c("Argument cont_method must be either 'Localmedian' or",
                    " 'Trapezoidalsurv'!"),
                  collapse = ""))
    }

    if (cont_method == "Localmedian") {
      u_y <- seq(eps, 1, eps)

      # The local median within each subinterval
      y_u <- stats::quantile(
        object$distr_y$margins[[1]], probs = u_y - eps / 2
        )

      # Compute the integral, compute the conditional CDF at the right endpoints
      # of the intervals corresponding to each evaluation point first (the last
      # has to be 1)
      uu <- cbind(
        sapply(
          u_y[-length(u_y)],
          function(u) eval_at_u_y(
            rep(u, n_obs), object$model, n_obs, cond_x, y_type
          )
        ),
        rep(1, n_obs)
      )

      # Compute the
      udiff <- t(apply(uu, 1, diff))
      uu[, 2:length(u_y)] <- udiff

      uu %*% y_u

    } else {
      u_y <- seq(0, 1, eps)

      # Inverse CDF of Y at evaluation points
      y_u <- stats::quantile(object$distr_y$margins[[1]], probs = u_y)

      # Compute the conditional survival function at each evaluation point
      surv <- sapply(u_y, function(u) 1 - eval_at_u_y(rep(u, n_obs),
                                                      object$model, n_obs,
                                                      cond_x, y_type))

      (stats::quantile(y_u, 0) +
          sapply(1:(length(u_y) - 1),
                 function(j) (surv[, j] + surv[, j + 1]) / 2)
        %*% (y_u[-1] - y_u[-length(y_u)]))

    }
  } else {

    #################
    # Discrete case #
    #################

    if (is.null(object$y)) {
      stop(
        "Must supply y-values used for fitting the model, when y is discrete"
      )
    }

    yval <- sort(unique(object$y))
    u_y <- object$distr_y$transform(yval, "d")
    uu <- vapply(
      1:(nrow(u_y)),
      function(i) eval_at_u_y(
        matrix(rep(u_y[i, ], n_obs), ncol = 2, byrow = TRUE), object$model,
        n_obs, cond_x, y_type = "d"
      ), FUN.VALUE = matrix(0, nrow = n_obs, ncol = 2)
    )

    p_y <- sapply(
      seq_len(nrow(u_y)),
      function(k) apply(uu[, , k], 1, function(x) diff(rev(x)))
    )

    p_y %*% yval

  }
}



.compute_conditionals <- function(u, model) {

  # Helper function that computes all of the transformed variables at every
  # level of the pair-copula model

  pick_u_cols <- function(var_inds, u, var_types) {
    # Small function that just selects the relevant columns of U for a specific
    # variable pair, used when computing the transformed variables for the
    # first tree, to avoid that the code in the loop below becomes too messy.
    if (any(var_types[var_inds] == "d")) {
      u[, c(var_inds, length(var_types) + cumsum(var_types == "d")[var_inds])]
    } else {
      u[, var_inds]
    }
  }

  # Create the hash table that we will store the transformed variables in,
  # and extract the R-vine matrix into a separate object to make the code in
  # the loop below easier to read.
  transformed_variables <- new.env(hash = TRUE)
  vinemat <- rvinecopulib::as_rvine_matrix(model$structure)
  d <- model$structure$d

  # First, add all the original variables
  for (j in 1:d) {
    assign(.make_hfunc_key(j, c()),
           pick_u_cols(c(j), u, model$var),
           envir = transformed_variables)
  }

  for (t in 1:(d - 1)) {
    for (e in 1:(d - t)) {
      # Extract the U-pair, and make the keys where we will insert the
      # new transformed variables we compute.
      if (t > 1) {
        # From the second tree, there are conditional variables D, stored in
        # the variable cond_set, and the variables U are new located in the
        # hash-table, so we'll have to extract them from there.
        cond_set <- vinemat[1:(t - 1), e]
        u_i_t <- .weave_transformed(
          get(.make_hfunc_key(vinemat[d + 1 - e, e], cond_set),
              envir = transformed_variables),
          get(.make_hfunc_key(vinemat[t, e], cond_set),
              envir = transformed_variables)
        )

        key1 <- .make_hfunc_key(vinemat[d + 1 - e, e],
                                c(cond_set, vinemat[t, e]))
        key2 <- .make_hfunc_key(vinemat[t, e],
                                c(cond_set, vinemat[d + 1 - e, e]))

      } else {
        # In the second tree, there are no conditional variables D,
        # and the variables U are the empirical margins of the original
        # variables.
        u_i_t <- pick_u_cols(vinemat[c(d + 1 - e, t), e], u, model$var)
        key1 <- .make_hfunc_key(vinemat[d + 1 - e, e], vinemat[t, e])
        key2 <- .make_hfunc_key(vinemat[t, e], vinemat[d + 1 - e, e])
      }

      # Compute and insert the new transformed variables
      assign(key1,
             .bicop_2_hbicop(u_i_t,
                             model$pair_copulas[[t]][[e]], cond_var = 2,
                             return_u_minus = model$var[vinemat[d + 1 - e,
                                                                e]] == "d"),
             envir = transformed_variables)

      assign(key2,
             .bicop_2_hbicop(u_i_t,
                             model$pair_copulas[[t]][[e]], cond_var = 1,
                             return_u_minus = model$var[vinemat[t, e]] == "d"),
             envir = transformed_variables)
    }
  }
  transformed_variables
}
