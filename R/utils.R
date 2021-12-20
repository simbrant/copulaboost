.make_hfunc_key <- function(i, cond_set) {

  #
  # This function creates a string key to be used when organising the evaluated
  # h-functions of every level of an R-vine pair copula model in a hash table,
  # using R's built-in enviroments.
  #

  cond_set <- sort(cond_set)
  paste0(c("(", i, "| ",
           paste0(c(sapply(cond_set[-length(cond_set)],
                           function(cond_var) paste0(c(cond_var, ", "),
                                                     collapse = "")),
                    cond_set[length(cond_set)]),
                  collapse = ""),
           ")"), collapse = "")
}

.weave_transformed <- function(u1, u2) {

  # This function weaves two transformed variables (so that the resulting
  # variable has the form [u+, u-], if any of the two margins are discrete,
  # u- only contains columns from discrete variables),
  # where u_j will be a 2xn matrix if the conditioned variable
  # in the j-th transformed variable is discrete, and a numeric vector
  # of length n if continuous.

  u <- cbind(u1, u2)

  if (is.null(ncol(u1)) || ncol(u1) == 1) {
    u
  } else if (is.null(ncol(u2)) || ncol(u2) == 1) {
    u[, c(1, 3, 2)]
  } else {
    u[, c(1, 3, 2, 4)]
  }
}

.bicop_2_hbicop <- function(u, bicop_obj, cond_var=2, return_u_minus=F) {

  # Function that lets you compute h-functions, without having to
  # specify the u_1^- when computing C(u_1 | u_2), when u_1 is
  # discrete. This is because u_1^- is redundant in that case, but rvinecopulibs
  # hbicop() function demands that it is provided. In addition, the specifying
  # return_u_minus = T, will make the function output the n x 2 matrix
  # [C(u_2 | u_1), C(u_2^- | u_1)] if cond_var = 1, or
  # [C(u_1 | u_2), C(u_1^- | u_2)] if cond_var = 2.

  if (!bicop_obj$var_types[c(2, 1)[cond_var]] == "d") {
    # In this case (that the conditioned variable is continuous), the
    # h-function in rvinecopulib can be called directly, as it then behaves as
    # expected

    if (return_u_minus) {
      cbind(rvinecopulib::hbicop(u, cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types),
            rvinecopulib::hbicop(u, cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types))
    } else {
      rvinecopulib::hbicop(u, cond_var = cond_var, family = bicop_obj$family,
                           rotation = bicop_obj$rotation,
                           parameters = bicop_obj$parameters,
                           var_types = bicop_obj$var_types)
    }
  } else {

    # This is more complicated. There are four_cases.

    u_columns_1 <- switch(1 + 2 * (cond_var - 1) + 1 * (ncol(u) == 4),
                          c(1, 2, 1, 3),
                          c(1, 2, 3, 4),
                          c(1, 2, 3, 2),
                          c(1, 2, 3, 4))

    if (return_u_minus) {

      u_columns_2 <- switch(1 + 2 * (cond_var - 1) + 1 * (ncol(u) == 4),
                            c(1, 3, 1, 3),
                            c(1, 4, 3, 4),
                            c(3, 2, 3, 2),
                            c(3, 2, 3, 4))

      cbind(rvinecopulib::hbicop(u[, u_columns_1], cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types),
            rvinecopulib::hbicop(u[, u_columns_2], cond_var = cond_var,
                                 family = bicop_obj$family,
                                 rotation = bicop_obj$rotation,
                                 parameters = bicop_obj$parameters,
                                 var_types = bicop_obj$var_types))

    } else {
      rvinecopulib::hbicop(u[, u_columns_1], cond_var = cond_var,
                           family = bicop_obj$family,
                           rotation = bicop_obj$rotation,
                           parameters = bicop_obj$parameters,
                           var_types = bicop_obj$var_types)
    }
  }
}

.approx_cond_exp <- function(y, u, v, m, r, truncate=T) {

  first_six_moments <- function(mu, sigma) {
    cbind(mu,
          mu**2 + sigma**2,
          mu**3 + 3 * mu * sigma**2,
          mu**4 + 6 * mu**2 * sigma**2 + 3 * sigma**4,
          mu**5 + 10 * mu**3 * sigma**2 + 15 * mu * sigma**4,
          mu**6 + 15 * mu**4 * sigma**2 + 45 * mu**2 + sigma**4 +
            15 * sigma**6)
  }

  if (! m %in% 1:6) {
    stop("m must be one of 1, 2, 3, 4, 5, or 6")
  }

  n <- length(y)

  gammas <- lm(
    y~. - 1,
    data = as.data.frame(
      sapply(seq_len(m + 1), function(j) qnorm(u) ** (j - 1))
    )
  )$coef

  if (any(is.na(gammas))){
    print(gammas)
    print(y)
    ptint(u)
  }
  means <- c(r) * qnorm(v)
  sd <- rep(sqrt(1 - r**2), n)

  if (m == 1) {
    est <- (
      gammas[1] +
        as.numeric(first_six_moments(means, sd)[, 1] * gammas[-1])
    )
  } else {
    est <- (
      gammas[1] +
        as.numeric(first_six_moments(means, sd)[, 1:m] %*%  gammas[-1])
    )
  }
  if (truncate) {
    est[est > 1] <- 1
    est[est < -1] <- -1
    est
  } else {
    est
  }
}

.get_hfunc <- function(j, cond_set, conditionals) {

  u <- get(.make_hfunc_key(j, cond_set),
           envir = conditionals)
  if (is.null(dim(u))) {
    u
  } else {
    u[, 1]
  }
}

.get_valid_edges <- function(t, rvine_mat, flipped_mat=F) {

  # This function finds all valid values that can
  # be entered into m_{d + 1 -t, 1} for t from 2 to d-1,
  # given an R-vine matrix

  if (t == 1) {
    seq(ncol(rvine_mat) - 1)
  } else {

    # Todo: This will always be false, so I should perhaps
    # rewrite the code below.
    if (!flipped_mat) {
      rvine_mat <- rvine_mat[seq(nrow(rvine_mat), 1), ]
    }

    k <- nrow(rvine_mat) + 2 - t
    match_vec <- rvine_mat[k:nrow(rvine_mat), 1]

    valid_values <- c()

    for (l in 2:(nrow(rvine_mat) + 1 - t)) {
      set_l_k <- rvine_mat[k:nrow(rvine_mat), l]
      if (setequal(set_l_k, match_vec)) {
        valid_values <- c(valid_values, rvine_mat[l, l])
      }
      if (k == nrow(rvine_mat)) {
        tilde_set_l_k <- rvine_mat[l, l]
      } else {
        tilde_set_l_k <- c(rvine_mat[l, l],
                           rvine_mat[(k + 1):nrow(rvine_mat), l])
      }
      if (setequal(tilde_set_l_k, match_vec)) {
        valid_values <- c(valid_values, rvine_mat[k, l])
      }
    }

    valid_values
  }
}

.compute_distrbs <- function(x, x_type, parametric=FALSE) {

  #
  # This function takes a matrix x, and returns a list containing empirical
  # cdf functions for each column, as well as a function that can transform
  # a new matrix that has the same number of columns, by these empirical cdf
  # functions. If parametric = TRUE, the function will fit and return marginal
  # cdfs that are
  #

  # Check, in case x is 1d
  if (is.null(ncol(x))) {
    x <- as.matrix(x)
  }

  ecdf_np1 <- function(x, x_type) {
    #
    # copy of stats::ecdf, but with a different
    # dividing constant: (1/n) is replaced by (1/(n+1)), if x is continuous
    #
    x <- sort(x)
    n <- length(x)

    if (n < 1) {
      stop("'x' must have 1 or more non-missing values")
    }

    weight <- switch(x_type, c = 1 / (n + 1), d = 1 / n)
    vals <- unique(x)
    rval <- approxfun(
      vals, cumsum(tabulate(match(x, vals))) * weight, method = "constant",
      yleft = 0, yright = 1, f = 0, ties = "ordered"
    )
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    rval
  }

  margins <- sapply(
    seq(ncol(x)), function(j) {
      if (x_type[j] == "d") {
        if (all(sort(unique(x[, j])) == c(0, 1)) & parametric) {
          bernoulli_marg(x[, j])
        } else {
          ecdf_np1(x[, j], x_type[j])
        }
      } else {
        if (parametric) {
          normal_marg(x[, j])
        } else {
          ecdf_np1(x[, j], x_type[j])
        }
      }
    }
  )

  transform <- function(x, var_types_x) {
    # Check, in case x is 1d
    if (is.null(ncol(x))) {
      x <- as.matrix(x)
    }
    distr_plus <- sapply(seq_len(ncol(x)), function(i) margins[[i]](x[, i]))
    if (any(var_types_x == "d")) {
      distr_min <- sapply(which(var_types_x == "d"),
                          function(i) margins[[i]](x[, i] - 1))
      return(cbind(distr_plus, distr_min))
    } else{
      return(distr_plus)
    }
  }

  list(margins = margins, transform = transform)

}

.extract_margs <- function(f, mrgs) {
  margins <- sapply(mrgs, function(j) f$margins[[j]])

  transform <- function(x, var_types_x) {
    # Check, in case x is 1d
    if (is.null(ncol(x))) {
      x <- as.matrix(x)
    }
    distr_plus <- sapply(seq_len(ncol(x)), function(i) margins[[i]](x[, i]))
    if (any(var_types_x == "d")) {
      distr_min <- sapply(which(var_types_x == "d"),
                          function(i) margins[[i]](x[, i] - 1))
      return(cbind(distr_plus, distr_min))
    } else{
      return(distr_plus)
    }
  }

  list(margins = margins, transform = transform)
}

.pick_u_cols <- function(var_inds, u, var_types) {
  # Small function that just selects the relevant columns of U for a specific
  # variable pair, used when computing the transformed variables for the
  # first tree, to avoid that the code in the loop below becomes too messy.
  if (any(var_types[var_inds] == "d")) {
    u[, c(var_inds, length(var_types) + cumsum(var_types == "d")[var_inds])]
  } else {
    u[, var_inds]
  }
}

normal_marg <- function(x) {
  f <- function(z) {
    mu <- attr(sys.function(), "mu")
    sigma <- attr(sys.function(), "sigma")
    pnorm(z, mu, sigma)
  }
  attr(f, "mu") <- mean(x)
  attr(f, "sigma") <- sd(x)
  class(f) <- c("normal_marg")
  attr(f, "call") <- sys.call()
  f
}

quantile.normal_marg <- function(f, probs) {
  qnorm(probs, attr(f, "mu"), attr(f, "sigma"))
}

bernoulli_marg <- function(x) {
  f <- function(z) {
    prob <- attr(sys.function(), "prob")
    pbinom(z, 1, prob)
  }
  attr(f, "prob") <- mean(x)
  class(f) <- c("bernoulli_marg")
  attr(f, "call") <- sys.call()
  f
}

quantile.bernoulli_marg <- function(f, probs) {
  qbinom(probs, 1, attr(f, "prob"))
}
