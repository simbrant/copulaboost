
.select_n_cov <- function(eta, y, x, x_type, n_cov, m, ystar_cont, ml_update,
                          max_ml_scale, par_method, dx, dy, cl) {
  if (is.null(cl)) {
    .select_n_cov_serial(eta, y, x, x_type, n_cov, m, ystar_cont, ml_update,
                         max_ml_scale, par_method, dx, dy)
  } else {
    .select_n_cov_parallel(eta, y, x, x_type, n_cov, m, ystar_cont, ml_update,
                           max_ml_scale, par_method, dx, dy, cl)
  }
}

.select_n_cov_parallel <- function(eta, y, x, x_type, n_cov, m, ystar_cont,
                                   ml_update, max_ml_scale, par_method, dx, dy,
                                   cl) {

  model <- list(
    ystar = y - plogis(eta),
    vtyp_sel = rep(NULL, n_cov + 1),
    selected_covs = rep(0, n_cov),
    dy = if (is.null(dy)) {
      .compute_distrbs(y - plogis(eta), if (ystar_cont) "c" else "d")
    } else {
      dy
    },
    dx = if (is.null(dx)) {
      .compute_distrbs(x, x_type)
    } else {
      dx
    },
    v_mat = rvinecopulib::as_rvine_matrix(
      rvinecopulib::dvine_structure(c(n_cov + 1, seq(n_cov)))
    ),
    p_cops = lapply(
      1:n_cov, function(l) vector("list", length = n_cov + 1 - l)
    ),
    transformed_variables = new.env(hash = TRUE)
  )

  # Add the margin of y to the object that holds the margins
  assign(.make_hfunc_key(model$v_mat[n_cov + 1, 1], c()),
         if (ystar_cont) {
           model$dy$margins[[1]](model$ystar)
         } else {
           cbind(model$dy$margins[[1]](model$ystar),
                 model$dy$margins[[1]](model$ystar - 1e-3))
         },
         envir = model$transformed_variables)
  model$vtyp_sel[n_cov + 1] <- (if (ystar_cont) "c" else "d")

  # Greedily update the model
  for (nc in seq(n_cov)) {

    rs_nc <- rep(-Inf, ncol(x))

    loop_inds <- (
      if (nc == 1) {
        seq_len(ncol(x))
      } else{
        seq(ncol(x))[-model$selected_covs[seq_len(nc - 1)]]
      }
    )

    loop_update <- function(j){
      model <- .update_model(j, nc, model, x, x_type, n_cov, par_method)

      # Compute approximation to E(ystar | x_s_curr)
      if (ystar_cont == T) {

        utop <- cbind(
          .get_hfunc(
            j = model$v_mat[n_cov + 1, 1],
            cond_set = if (nc > 1) model$v_mat[seq_len(nc - 1), 1] else c(),
            conditionals = model$transformed_variables
          ),
          .get_hfunc(
            j = model$v_mat[nc, 1],
            cond_set = if (nc > 1) model$v_mat[seq_len(nc - 1), 1] else c(),
            conditionals = model$transformed_variables
          )
        )

        h_k <- .approx_cond_exp(
          model$ystar, utop[, 1], utop[, 2], m = m,
          if (sum(dim(model$p_cops[[nc]][[1]]$parameters)) > 0) {
            as.numeric(model$p_cops[[nc]][[1]]$parameters)
          } else {
            0
          }
        )

      } else {
        yval <- sort(unique(model$ystar))
        u_y <- cbind(model$dy$margins[[1]](yval),
                     model$dy$margins[[1]](yval - 1))

        eval_at_u_y <- function(u_y) {
          u_ <- .weave_transformed(
            u_y,
            get(.make_hfunc_key(model$v_mat[1, 1], c()),
                envir = model$transformed_variables)
          )
          h <- .bicop_2_hbicop(u_, model$p_cops[[1]][[1]],
                               return_u_minus = T)
          if (nc > 1) {
            for (t in seq(2, nc)) {
              u_ <- .weave_transformed(
                h, get(.make_hfunc_key(model$v_mat[t, 1],
                                       model$v_mat[seq_len(t - 1), 1]),
                       envir = model$transformed_variables)
              )
              h <- .bicop_2_hbicop(u_, model$p_cops[[t]][[1]],
                                   return_u_minus = T)
            }
          }

          h
        }

        uu <- vapply(1:(nrow(u_y)),
                     function(i) eval_at_u_y(matrix(rep(u_y[i, ], length(y)),
                                                    ncol = 2, byrow = T)),
                     FUN.VALUE = matrix(0, nrow = length(y), ncol = 2))

        p_y <- sapply(seq_len(nrow(u_y)),
                      function(k) apply(uu[, , k], 1,
                                        function(x) diff(rev(x))))

        h_k <- p_y %*% yval
      }

      scale <- if (ml_update) {
        optim(
          fn = function(theta) -sum(dbinom(y, 1, plogis(eta + theta * h_k))),
          par = 1, method = "Brent", lower = 0, upper = max_ml_scale)$par
      } else {
        1
      }
      sum(dbinom(y, 1, plogis(eta + scale * h_k), log = T))
    }

    rs_nc[loop_inds] <- parallel::parSapply(cl, X = loop_inds,
                                            FUN = loop_update)

    # Update the model
    model <- .update_model(
      which.max(rs_nc), nc, model, x, x_type, n_cov, par_method
    )

  }

  model$selected_covs

}

.select_n_cov_serial<- function(eta, y, x, x_type, n_cov, m, ystar_cont,
                                ml_update, max_ml_scale, par_method, dx, dy) {

  model <- list(
    ystar = y - plogis(eta),
    vtyp_sel = rep(NULL, n_cov + 1),
    selected_covs = rep(0, n_cov),
    dy = if (is.null(dy)) {
      .compute_distrbs(y - plogis(eta), if (ystar_cont) "c" else "d")
    } else {
      dy
    },
    dx = if (is.null(dx)) {
      .compute_distrbs(x, x_type)
    } else {
      dx
    },
    v_mat = rvinecopulib::as_rvine_matrix(
      rvinecopulib::dvine_structure(c(n_cov + 1, seq(n_cov)))
    ),
    p_cops = lapply(
      1:n_cov, function(l) vector("list", length = n_cov + 1 - l)
    ),
    transformed_variables = new.env(hash = TRUE)
  )

  # Add the margin of y to the object that holds the margins
  assign(.make_hfunc_key(model$v_mat[n_cov + 1, 1], c()),
         if (ystar_cont) {
           model$dy$margins[[1]](model$ystar)
         } else {
           cbind(model$dy$margins[[1]](model$ystar),
                 model$dy$margins[[1]](model$ystar - 1e-3))
         },
         envir = model$transformed_variables)
  model$vtyp_sel[n_cov + 1] <- (if (ystar_cont) "c" else "d")

  # Greedily update the model
  for (nc in seq(n_cov)) {

    rs_nc <- rep(-Inf, ncol(x))

    loop_inds <- (
      if (nc == 1) {
        seq_len(ncol(x))
      } else{
        seq(ncol(x))[-model$selected_covs[seq_len(nc - 1)]]
      }
    )

    for (j in loop_inds) {
      model <- .update_model(j, nc, model, x, x_type, n_cov, par_method)

      # Compute approximation to E(ystar | x_s_curr)
      if (ystar_cont == T) {

        utop <- cbind(
          .get_hfunc(
            j = model$v_mat[n_cov + 1, 1],
            cond_set = if (nc > 1) model$v_mat[seq_len(nc - 1), 1] else c(),
            conditionals = model$transformed_variables
          ),
          .get_hfunc(
            j = model$v_mat[nc, 1],
            cond_set = if (nc > 1) model$v_mat[seq_len(nc - 1), 1] else c(),
            conditionals = model$transformed_variables
          )
        )

        h_k <- .approx_cond_exp(
          model$ystar, utop[, 1], utop[, 2], m = m,
          if (sum(dim(model$p_cops[[nc]][[1]]$parameters)) > 0) {
            as.numeric(model$p_cops[[nc]][[1]]$parameters)
          } else {
            0
          }
        )

      } else {
        yval <- sort(unique(model$ystar))
        u_y <- cbind(model$dy$margins[[1]](yval),
                     model$dy$margins[[1]](yval - 1))

        eval_at_u_y <- function(u_y) {
          u_ <- .weave_transformed(
            u_y,
            get(.make_hfunc_key(model$v_mat[1, 1], c()),
                envir = model$transformed_variables)
          )
          h <- .bicop_2_hbicop(u_, model$p_cops[[1]][[1]],
                               return_u_minus = T)
          if (nc > 1) {
            for (t in seq(2, nc)) {
              u_ <- .weave_transformed(
                h, get(.make_hfunc_key(model$v_mat[t, 1],
                                       model$v_mat[seq_len(t - 1), 1]),
                       envir = model$transformed_variables)
              )
              h <- .bicop_2_hbicop(u_, model$p_cops[[t]][[1]],
                                   return_u_minus = T)
            }
          }

          h
        }

        uu <- vapply(1:(nrow(u_y)),
                     function(i) eval_at_u_y(matrix(rep(u_y[i, ], length(y)),
                                                    ncol = 2, byrow = T)),
                     FUN.VALUE = matrix(0, nrow = length(y), ncol = 2))

        p_y <- sapply(seq_len(nrow(u_y)),
                      function(k) apply(uu[, , k], 1,
                                        function(x) diff(rev(x))))

        h_k <- p_y %*% yval
      }

      scale <- if (ml_update) {
        optim(
          fn = function(theta) -sum(dbinom(y, 1, plogis(eta + theta * h_k))),
          par = 1, method = "Brent", lower = 0, upper = max_ml_scale)$par
      } else {
        1
      }
      rs_nc[j] <- sum(dbinom(y, 1, plogis(eta + scale * h_k), log = T))
    }

    # Update the model
    model <- .update_model(
      which.max(rs_nc), nc, model, x, x_type, n_cov, par_method
    )

  }
  model$selected_covs
}

.update_model <- function(j_new, nc, model, x, x_type, n_cov, par_method) {
  # This function contains all the code that adds one new variable
  # to the d-vine model
  assign(.make_hfunc_key(nc, c()),
         if (x_type[j_new] == "c") {
           model$dx$margins[[j_new]](x[, j_new])
         } else {
           cbind(model$dx$margins[[j_new]](x[, j_new]),
                 model$dx$margins[[j_new]](x[, j_new] - 1e-2))},
         envir = model$transformed_variables)

  model$selected_covs[nc] <- j_new
  model$vtyp_sel[nc] <- x_type[j_new]

  for (k in seq(nc)) {

    u_curr <- .weave_transformed(
      get(.make_hfunc_key(model$v_mat[(n_cov + 1) - (nc - k),
                                      nc - (k - 1)],
                          model$v_mat[if (k == 1) c() else seq_len(k - 1),
                                      nc - (k - 1)]),
          envir = model$transformed_variables),
      get(.make_hfunc_key(model$v_mat[k, nc - (k - 1)],
                          model$v_mat[if (k == 1) c() else seq_len(k - 1),
                                      nc - (k - 1)]),
          envir = model$transformed_variables)
    )

    model$p_cops[[k]][[nc - (k - 1)]] <- rvinecopulib::bicop(
      u_curr,
      var_types = c(model$vtyp_sel[model$v_mat[(n_cov + 1) - (nc - k),
                                               nc - (k - 1)]],
                    model$vtyp_sel[model$v_mat[k, nc - (k - 1)]]),
      family_set = "gaussian", par_method = par_method
    )

    # Add transformed variables
    assign(.make_hfunc_key(model$v_mat[(n_cov + 1) - (nc - k),
                                       nc - (k - 1)],
                           c(model$v_mat[if (k == 1) c() else seq_len(k - 1),
                                         nc - (k - 1)],
                             model$v_mat[k, nc - (k - 1)])),
           .bicop_2_hbicop(u_curr,
                           bicop_obj = model$p_cops[[k]][[nc - (k - 1)]],
                           cond_var = 2,
                           return_u_minus = (
                             model$vtyp_sel[
                               model$v_mat[(n_cov + 1) - (nc - k),
                                           nc - (k - 1)]
                             ] == "d"
                           )
           ), envir = model$transformed_variables
    )
    assign(.make_hfunc_key(model$v_mat[k, nc - (k - 1)],
                           c(model$v_mat[if (k == 1) c() else seq_len(k - 1),
                                         nc - (k - 1)],
                             model$v_mat[(n_cov + 1) - (nc - k),
                                         nc - (k - 1)])),
           .bicop_2_hbicop(
             u_curr,
             bicop_obj = model$p_cops[[k]][[nc - (k - 1)]],
             cond_var = 1,
             return_u_minus = (
               model$vtyp_sel[model$v_mat[k, nc - (k - 1)]] == "d"
             )), envir = model$transformed_variables
    )
  }
  model
}
