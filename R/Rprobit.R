#' Rprobit Object
#'
#' @description
#' An object of the `Rprobit` R6 class consists of a [`mod`] object (defines the
#' model structure) and a [`choice_data`] object (interprets choice data).
#'
#' - Estimation can be performed in place via the method `Rprobit$fit()` or
#'   without modifying the original object via [`fit()`].
#'
#' - Methods can be used to analyze the model fit, see [`Rprobit-methods`].
#'
#' @param mod \[`mod`\]\cr
#' Model specification object.
#'
#' @param choice_data \[`choice_data`\]\cr
#' Choice data object.
#'
#' @section Fitting workflows:
#' Calling `object$fit()` uses the mutable R6 method and stores the fitted
#' quantities in `object` itself. Calling `fit(object)` first creates a deep
#' clone of `object`, then calls `$fit()` on that clone and returns it. Use
#' `fit(object)` when the original `Rprobit` object should remain unchanged.
#' For more details on R6-specific modify-in-place logic, see
#' <https://r-statistics.co/R6-Classes-in-R.html>.
#'
#' @export

Rprobit <- R6::R6Class(
  classname = "Rprobit",
  public = list(

    #' @description
    #' Create a new `Rprobit` object.
    #'
    #' @return
    #' A new `Rprobit` object.
    initialize = function(mod, choice_data) {
      checkmate::assert_r6(mod, classes = "mod")
      checkmate::assert_r6(choice_data, classes = "choice_data")
      private$.mod <- mod
      private$.choice_data <- choice_data
      invisible(self)
    },

    #' @description
    #' Fit the model in place.
    #'
    #' @param theta_init \[`numeric()` | `NULL`\]\cr
    #' Initial parameter vector. If `NULL`, an initialization procedure is used.
    #'
    #' @param cml_type \[`character(1)`\]\cr
    #' Defines the type of composite marginal likelihood used. Possible values:
    #' - `"probit"` for the full probit likelihood,
    #' - `"FP"` for full pairwise composite marginal likelihood,
    #' - `"AP"` for adjacent pairwise composite marginal likelihood,
    #' - `"AC"` for cyclical adjacent pairwise composite marginal likelihood.
    #'
    #' @param ... Additional arguments passed to [stats::nlm()].
    #'
    #' @param calc_hess \[`logical(1)`\]\cr
    #' If `TRUE`, calculate the analytical Hessian at the estimate. If `FALSE`,
    #' use the outer-product approximation returned by [ll_cml()]. The
    #' analytical Hessian can be considerably more expensive to compute.
    #'
    #' @param verbose \[`logical(1)`\]\cr
    #' If `TRUE`, progress messages are shown.
    #'
    #' @return
    #' The `Rprobit` object, invisibly.
    #'
    #' The following fields are updated:
    #' - `theta_init` with the initial values,
    #' - `theta_fit` with the estimates,
    #' - `ll_fit` with the log-likelihood value at the estimate,
    #' - `optimization_details` with the output of [stats::nlm()], the
    #'   log-likelihood value at the initial values, and initialization and
    #'   optimization runtimes in seconds,
    #' - `hessian_fit` with either the analytical or approximated Hessian and
    #'   `Vhat_fit` as returned by [ll_cml()] at the estimate,
    #' - `theta_vcov` with the covariance matrix of `theta_fit`.
    fit = function(
      theta_init = NULL, cml_type = "probit", ..., calc_hess = FALSE,
      verbose = TRUE
    ) {

      ### input checks
      if (!is.null(theta_init)) {
        checkmate::assert_numeric(
          theta_init, len = self$mod$length_theta, any.missing = FALSE,
          finite = TRUE
        )
      }
      checkmate::assert_choice(
        cml_type, choices = c("probit", "FP", "AP", "AC")
      )
      checkmate::assert_flag(calc_hess)
      checkmate::assert_flag(verbose)

      ### initialize parameters if needed
      sec_init <- 0
      init_type <- if (is.null(theta_init)) "auto" else "custom"
      if (is.null(theta_init)) {
        if (isTRUE(verbose)) {
          message("Looking for initial values...")
        }
        initialization_start <- proc.time()[["elapsed"]]
        tryCatch(
          private$.initialize_theta(),
          error = function(e) private$.theta_init <- NULL
        )
        sec_init <- unname(proc.time()[["elapsed"]] - initialization_start)
      } else {
        if (isTRUE(verbose)) {
          message("Using custom initial values...")
        }
        private$.theta_init <- theta_init
      }

      ### build differenced design matrices
      design_diff <- self$choice_data$build_design(
        mod = self$mod,
        difference = TRUE
      )

      ### break panel structure if no random effects are present
      if (self$mod$b_count[2L] == 0L) {
        design_diff <- unlist(lapply(design_diff, function(decider) {
          lapply(seq_along(decider$X), function(t) {
            list(X = list(decider$X[[t]]), y = decider$y[t])
          })
        }), recursive = FALSE)
      }

      ### check that initial values are valid; if not, try to use unit values
      theta_init <- self$theta_init
      init_eval <- NULL
      for (attempt in seq_len(2L)) {
        init_eval <- tryCatch({
          checkmate::assert_numeric(
            theta_init, any.missing = FALSE, finite = TRUE,
            len = self$mod$length_theta
          )
          theta_init <- self$mod$make_theta_unique(theta_init)
          parameter <- self$mod$build_parameter(theta_init)
          checkmate::assert_true(
            all(is.finite(unlist(parameter, use.names = FALSE)))
          )
          ll_cml(
            theta = theta_init,
            design_diff = design_diff,
            mod_list = self$mod$mod_list,
            cml_type = cml_type,
            calc_hess = FALSE
          )
        }, error = function(e) NULL)
        if (checkmate::test_number(init_eval, na.ok = FALSE, finite = TRUE)) {
          break
        }
        if (attempt == 2L) {
          stop(
            "Optimization failed to start; please check model defintion.",
            call. = FALSE
          )
        }
        warning(
          "Initialization failed; using unit fallback.",
          call. = FALSE
        )
        ### using unit fallback
        init_type <- "fallback"
        theta_init <- tryCatch(
          suppressWarnings(self$mod$build_theta(
            b = numeric(sum(self$mod$b_count)),
            Omega = diag(self$mod$b_count[2L]),
            Sigma = diag(self$mod$alt)
          )),
          error = function(e) NULL
        )
      }
      private$.theta_init <- theta_init

      ### if the model has no free parameters, just evaluate the likelihood
      if (self$mod$length_theta == 0L) {
        private$.theta_fit <- numeric()
        private$.ll_fit <- -as.numeric(init_eval)
        private$.hessian_fit <- matrix(numeric(), nrow = 0L, ncol = 0L)
        private$.Vhat_fit <- matrix(numeric(), nrow = 0L, ncol = 0L)
        private$.theta_vcov <- matrix(numeric(), nrow = 0L, ncol = 0L)
        if (isTRUE(verbose)) {
          message("Done.")
        }
        return(invisible(self))
      }

      ### optimize negative log-likelihood
      if (isTRUE(verbose)) {
        message("Starting numerical optimization...")
      }
      nlm_args <- utils::modifyList(
        ### check.analyticals = FALSE to prevent nlm fails for tiny deviations
        list(check.analyticals = FALSE), list(...)
      )
      optimization_start <- proc.time()[["elapsed"]]
      opt <- do.call(stats::nlm, c(list(
        f = ll_cml, p = self$theta_init, design_diff = design_diff,
        mod_list = self$mod$mod_list, cml_type = cml_type
      ), nlm_args))
      sec_fit <- unname(proc.time()[["elapsed"]] - optimization_start)

      ### save results
      theta_fit <- self$mod$make_theta_unique(opt$estimate)
      private$.theta_fit <- theta_fit
      private$.ll_fit <- -opt$minimum
      private$.optimization_details <- c(
        opt,
        list(
          "ll_init" = -as.numeric(init_eval),
          "init_type" = init_type,
          "sec_fit" = sec_fit,
          "sec_init" = sec_init
        )
      )

      ### calculate Hessian matrices and covariance matrix at the estimate
      if (isTRUE(verbose)) {
        hessian_type <- if (isTRUE(calc_hess)) "analytical" else "approximated"
        message("Calculating ", hessian_type, " Hessian at estimate...")
      }
      fit_eval <- ll_cml(
        theta = self$theta_fit,
        design_diff = design_diff,
        mod_list = self$mod$mod_list,
        cml_type = cml_type,
        calc_hess = calc_hess
      )
      hessian_attr <- if (isTRUE(calc_hess)) "hessian" else "hessian_approx"
      hessian_fit <- attr(fit_eval, hessian_attr, exact = TRUE)
      Vhat_fit <- attr(fit_eval, "Vhat")
      theta_names <- paste0("theta_", seq_len(self$mod$length_theta))
      dimnames(hessian_fit) <- list(theta_names, theta_names)
      dimnames(Vhat_fit) <- list(theta_names, theta_names)
      adjusted_hessian <- adjust_covariance_matrix(
        hessian_fit, label = "Hessian"
      )
      adjusted_Vhat <- adjust_covariance_matrix(
        Vhat_fit, label = "Vhat matrix"
      )
      if (cml_type == "probit") {
        ### full probit case, covariance matrix from information equality
        theta_vcov <- adjusted_Vhat$inverse
      } else {
        ### CML case, covariance matrix from sandwich estimator
        theta_vcov <- adjusted_hessian$inverse %*%
          adjusted_Vhat$mat %*%
          t(adjusted_hessian$inverse)
        theta_vcov <- (theta_vcov + t(theta_vcov)) / 2
      }
      dimnames(theta_vcov) <- list(theta_names, theta_names)

      ### save results
      private$.hessian_fit <- hessian_fit
      private$.Vhat_fit <- Vhat_fit
      private$.theta_vcov <- theta_vcov
      if (isTRUE(verbose)) {
        message("Done.")
      }
      invisible(self)
    },

    #' @description
    #' Print a summary of the `Rprobit` object.
    print = function() {
      cat("<Rprobit object>\n\n")
      private$.choice_data$print(what = "dimensions")
      private$.mod$print(what = c("formula", "re", "utility_equations"))
      invisible(self)
    }
  ),

  private = list(
    .mod = NULL,
    .choice_data = NULL,
    .theta_init = NULL,
    .theta_fit = NULL,
    .ll_fit = NULL,
    .optimization_details = NULL,
    .hessian_fit = NULL,
    .Vhat_fit = NULL,
    .theta_vcov = NULL,

    # Link between mod object and initialization procedure.
    .initialize_theta = function() {

      ### trivial case
      if (self$mod$length_theta == 0L) {
        private$.theta_init <- numeric()
        return(invisible(NULL))
      }

      ### prepare required information
      mod <- self$mod
      P <- lengths(read_formula_rhs(mod$formula))
      P[2L] <- P[2L] + as.integer(read_formula_asc(mod$formula))
      design <- self$choice_data$build_design(mod = mod, difference = FALSE)
      first_X <- design[[1L]]$X[[1L]]
      J <- nrow(first_X)
      ref_alt_index <- match(mod$ref_alt, rownames(first_X))
      re_cols <- match(mod$var_names$random, colnames(first_X))

      ### initialize unrestricted model
      raw_init <- as.numeric(initialize(
        data = design, J = J, P = P, ref_alt_index = ref_alt_index,
        re_cols = re_cols
      ))
      b_count <- P[1L] + P[2L] * (J - 1L) + P[3L] * J
      O_count <- length(re_cols) * (length(re_cols) + 1L) / 2L
      L_count <- J * (J - 1L) / 2L
      b_init <- matrix(raw_init[seq_len(b_count)], ncol = 1)
      O_index <- b_count + seq_len(O_count)
      O_init <- matrix(raw_init[O_index], ncol = 1)
      L_index <- b_count + O_count + seq_len(L_count)
      L_diff_init <- vech_chol_to_chol(raw_init[L_index])
      Sigma_diff_init <- tcrossprod(L_diff_init)

      ### map error term covariance matrix
      H_L <- unname(mod$H_L)
      f_L <- as.numeric(mod$f_L)
      H_b <- unname(mod$H_b)
      f_b <- as.numeric(mod$f_b)
      difference <- cbind(-1, diag(J - 1L))
      L_init <- matrix(0, J, J)
      L_init[-1L, -1L] <- L_diff_init
      L_init <- L_init[lower.tri(L_init, diag = TRUE)]
      theta_L_init <- numeric()
      if (ncol(H_L) > 0L) {
        theta_L_init <- stats::lm.fit(H_L, L_init - f_L)$coefficients
        theta_L_init[!is.finite(theta_L_init)] <- 0
      }

      ### find common utility scale implied by restrictions
      inverse_scale <- 1
      theta_L_unscaled <- theta_L_init
      if (ncol(H_L) == 0L) {
        fixed_L_diff <- difference %*% vech_chol_to_chol(f_L)
        fixed_Sigma_diff <- tcrossprod(fixed_L_diff)
        fixed_Sigma_size <- sum(fixed_Sigma_diff^2)
        inverse_scale_sq <- if (fixed_Sigma_size > 0) {
          sum(fixed_Sigma_diff * Sigma_diff_init) / fixed_Sigma_size
        } else {
          NA_real_
        }
        if (is.finite(inverse_scale_sq) && inverse_scale_sq > 0) {
          inverse_scale <- sqrt(inverse_scale_sq)
        }
      } else {
        b_residuals <- cbind(fixed = f_b, init = as.numeric(b_init))
        if (nrow(H_b) > 0L && ncol(H_b) > 0L) {
          b_residuals <- stats::lm.fit(H_b, b_residuals)$residuals
        }
        fixed_b_resid <- b_residuals[, 1L]
        b_init_resid <- b_residuals[, 2L]
        Sigma_size <- max(sum(Sigma_diff_init^2), 1)
        b_size <- max(sum(b_residuals^2), 1)
        scale_solution <- c(log_inverse_scale = 0, theta_L_init)
        scale_loss <- function(candidate) {
          inverse_scale <- exp(candidate[1L])
          theta_L <- candidate[-1L]
          L <- H_L %*% theta_L + inverse_scale * f_L
          L_diff <- difference %*% vech_chol_to_chol(L)
          Sigma_error <- tcrossprod(L_diff) - Sigma_diff_init
          b_error <- b_init_resid - inverse_scale * fixed_b_resid
          Sigma_loss <- sum(Sigma_error^2) / Sigma_size
          b_loss <- sum(b_error^2) / b_size
          loss <- Sigma_loss + b_loss
          if (is.finite(loss)) loss else 1e100
        }
        scale_fit <- tryCatch(
          stats::nlm(
            scale_loss, p = scale_solution, gradtol = 1e-12,
            steptol = 1e-12, stepmax = 1, iterlim = 500L
          ),
          error = function(error) NULL
        )
        valid_scale <- !is.null(scale_fit) &&
          all(is.finite(c(scale_fit$estimate, scale_fit$minimum))) &&
          scale_fit$minimum <= scale_loss(scale_solution)
        if (valid_scale) {
          inverse_scale_fit <- exp(scale_fit$estimate[1L])
          b_sets_scale <- any(abs(fixed_b_resid) > sqrt(.Machine$double.eps))
          valid_scale <- is.finite(inverse_scale_fit) &&
            (b_sets_scale || inverse_scale_fit > .Machine$double.eps)
        }
        if (valid_scale) scale_solution <- scale_fit$estimate
        inverse_scale <- exp(scale_solution[1L])
        theta_L_unscaled <- scale_solution[-1L]
      }

      ### use common utility scale for all model parameters
      ### apply model restrictions
      suppress_warnings <- mod$suppress_warnings
      mod$suppress_warnings <- TRUE
      on.exit(mod$suppress_warnings <- suppress_warnings, add = TRUE)
      scale_candidates <- unique(list(
        c(inverse_scale, theta_L_unscaled), c(1, theta_L_init)
      ))
      theta_init <- NULL
      theta_L <- numeric()
      for (index in seq_along(scale_candidates)) {
        candidate <- scale_candidates[[index]]
        utility_scale <- 1 / candidate[1L]
        theta_L <- utility_scale * candidate[-1L]
        L <- H_L %*% theta_L + f_L
        Sigma <- tcrossprod(vech_chol_to_chol(L))
        Omega <- NULL
        if (length(O_init) > 0L) {
          Omega <- tcrossprod(vech_chol_to_chol(O_init))
          O_entries <- which(lower.tri(Omega, diag = TRUE), arr.ind = TRUE)
          O_offdiag <- which(O_entries[, 1L] != O_entries[, 2L])
          diagonal_O <- all(mod$H_O[O_offdiag, , drop = FALSE] == 0) &&
            all(mod$f_O[O_offdiag, , drop = FALSE] == 0)
          if (diagonal_O) {
            prec <- tryCatch(solve(Omega), error = function(error) NULL)
            cond_var <- if (is.null(prec)) diag(Omega) else 1 / diag(prec)
            invalid_var <- any(!is.finite(cond_var)) || any(cond_var <= 0)
            if (invalid_var) cond_var <- diag(Omega)
            Omega[row(Omega) != col(Omega)] <- 0
            diag(Omega) <- cond_var
          }
          Omega <- utility_scale^2 * Omega
        }
        b <- utility_scale * b_init

        ### apply model restrictions
        last_candidate <- index == length(scale_candidates)
        theta_init <- if (last_candidate) {
          mod$build_theta(b = b, Omega = Omega, Sigma = Sigma)
        } else {
          tryCatch(
            mod$build_theta(b = b, Omega = Omega, Sigma = Sigma),
            error = function(error) NULL
          )
        }
        if (!is.null(theta_init)) break
      }
      if (length(theta_L) > 0L) {
        theta_L_index <- mod$length_theta_b + mod$length_theta_O +
          seq_along(theta_L)
        theta_init[theta_L_index] <- theta_L
      }

      ### save initial parameter vector
      private$.theta_init <- theta_init
    }
  ),

  active = list(

    #' @field mod (mutable) \[`mod`\]\cr
    #' Model specification object.
    mod = function(value) {
      if (missing(value)) {
        private$.mod
      } else if (identical(value, private$.mod)) {
        # Nested replacement writes the same R6 object back to this binding.
        invisible(private$.mod)
      } else {
        stop("Field `mod` is read-only.", call. = FALSE)
      }
    },

    #' @field choice_data (mutable) \[`choice_data`\]\cr
    #' Choice data object.
    choice_data = function(value) {
      if (missing(value)) {
        private$.choice_data
      } else if (identical(value, private$.choice_data)) {
        # Nested replacement writes the same R6 object back to this binding.
        invisible(private$.choice_data)
      } else {
        stop("Field `choice_data` is read-only.", call. = FALSE)
      }
    },

    #' @field theta_init (read only) \[`numeric()` | `NULL`\]\cr
    #' Initial parameter vector used for fitting.
    theta_init = function(value) {
      if (missing(value)) {
        private$.theta_init
      } else {
        stop("Field `theta_init` is read-only.", call. = FALSE)
      }
    },

    #' @field theta_fit (read only) \[`numeric()` | `NULL`\]\cr
    #' Model estimate.
    theta_fit = function(value) {
      if (missing(value)) {
        private$.theta_fit
      } else {
        stop("Field `theta_fit` is read-only.", call. = FALSE)
      }
    },

    #' @field ll_fit (read only) \[`numeric(1)` | `NULL`\]\cr
    #' Log-likelihood value at the model estimate.
    ll_fit = function(value) {
      if (missing(value)) {
        private$.ll_fit
      } else {
        stop("Field `ll_fit` is read-only.", call. = FALSE)
      }
    },

    #' @field optimization_details (read only) \[`list` | `NULL`\]\cr
    #' Details of the numerical optimization. After fitting, this list contains
    #' all components returned by [stats::nlm()], `ll_init` (the log-likelihood
    #' at `theta_init`), `init_type` (`"custom"`, `"auto"`, or `"fallback"`),
    #' `sec_fit` (the [stats::nlm()] runtime in seconds), and `sec_init` (the
    #' `.initialize_theta()` runtime in seconds).
    optimization_details = function(value) {
      if (missing(value)) {
        private$.optimization_details
      } else {
        stop("Field `optimization_details` is read-only.", call. = FALSE)
      }
    },

    #' @field hessian_fit (read only) \[`matrix` | `NULL`\]\cr
    #' Hessian of the negative fitting criterion at `theta_fit`. This is the
    #' analytical Hessian when fitting with `calc_hess = TRUE` and an
    #' outer-product approximation otherwise.
    hessian_fit = function(value) {
      if (missing(value)) {
        private$.hessian_fit
      } else {
        stop("Field `hessian_fit` is read-only.", call. = FALSE)
      }
    },

    #' @field Vhat_fit (read only) \[`matrix` | `NULL`\]\cr
    #' Estimate of the variance of the score at `theta_fit`.
    Vhat_fit = function(value) {
      if (missing(value)) {
        private$.Vhat_fit
      } else {
        stop("Field `Vhat_fit` is read-only.", call. = FALSE)
      }
    },

    #' @field theta_vcov (read only) \[`matrix` | `NULL`\]\cr
    #' Covariance matrix of `theta_fit`.
    theta_vcov = function(value) {
      if (missing(value)) {
        private$.theta_vcov
      } else {
        stop("Field `theta_vcov` is read-only.", call. = FALSE)
      }
    }
  )
)

#' Fit an Rprobit Object Without Modifying It in Place
#'
#' @description
#' `fit()` is the functional counterpart to the mutable `object$fit()` method.
#' It deep-clones `object`, fits the clone, and returns the fitted clone,
#' see details.
#'
#' @param object \[`Rprobit`\]\cr
#' An `Rprobit` object.
#'
#' @param theta_init \[`numeric()` | `NULL`\]\cr
#' Initial parameter vector. If `NULL`, an initialization procedure is used.
#'
#' @param cml_type \[`character(1)`\]\cr
#' Defines the type of composite marginal likelihood used. Possible values:
#' - `"probit"` for the full probit likelihood,
#' - `"FP"` for full pairwise composite marginal likelihood,
#' - `"AP"` for adjacent pairwise composite marginal likelihood,
#' - `"AC"` for cyclical adjacent pairwise composite marginal likelihood.
#'
#' @param ... Additional arguments passed to [stats::nlm()].
#'
#' @param calc_hess \[`logical(1)`\]\cr
#' If `TRUE`, calculate the analytical Hessian at the estimate. If `FALSE`,
#' use the outer-product approximation returned by [ll_cml()]. The analytical
#' Hessian can be considerably more expensive to compute.
#'
#' @param verbose \[`logical(1)`\]\cr
#' If `TRUE`, progress messages are shown.
#'
#' @return
#' A fitted clone of `object`.
#'
#' @inheritSection Rprobit Fitting workflows
#'
#' @export

fit <- function(
    object, theta_init = NULL, cml_type = "probit", ..., calc_hess = FALSE,
    verbose = TRUE
) {
  checkmate::assert_r6(object, classes = "Rprobit")
  fitted <- object$clone(deep = TRUE)
  fitted$fit(
    theta_init = theta_init, cml_type = cml_type, calc_hess = calc_hess, ...,
    verbose = verbose
  )
  fitted
}
