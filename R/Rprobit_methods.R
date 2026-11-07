#' Rprobit Methods
#'
#' @description
#' Methods for fitted `Rprobit` objects.
#'
#' @param object \[`Rprobit`\]\cr
#' An `Rprobit` object.
#'
#' @param x \[`Rprobit` | `summary.Rprobit`\]\cr
#' An `Rprobit` object for `plot()` or a `summary.Rprobit` object for `print()`.
#'
#' @param ... \[`any`\]\cr
#' Further arguments passed to or from other methods.
#'
#' @param k \[`numeric(1)`\]\cr
#' Numeric penalty per parameter.
#'
#' @param parm \[`character()` | `numeric()` | `NULL`\]\cr
#' Coefficient names or positions selected by `confint()` or `plot()`.
#'
#' @param level \[`numeric(1)`\]\cr
#' Confidence level.
#'
#' @param cdf \[`logical(1)`\]\cr
#' If `TRUE`, plot the cumulative distribution function instead of the
#' probability density function.
#'
#' @param margin \[`numeric(1)`\]\cr
#' Non-negative factor by which the plotting range is extended.
#'
#' @param newdata \[`choice_data` | `NULL`\]\cr
#' Optional data for prediction. If `NULL`, the data stored in `object` are
#' used.
#'
#' @param type \[`character(1)`\]\cr
#' Type of prediction to return:
#' - `"probabilities"` for all choice probabilities,
#' - `"choice"` for the probability of the observed choice,
#' - `"frequencies"` for the frequencies of the predicted choices,
#' - `"confusion"` for a confusion matrix of observed and predicted choices.
#'
#' @param digits \[`integer(1)`\]\cr
#' Minimum number of significant digits used for printing.
#'
#' @param signif.stars \[`logical(1)`\]\cr
#' If `TRUE`, significance stars are printed for the coefficient tests.
#'
#' @details
#' `AIC()` computes \eqn{-2 \ell + k p}, where \eqn{\ell} is the fitted
#' log-likelihood, \eqn{p} is the number of free model parameters, and \eqn{k}
#' is the non-negative penalty applied to each free parameter. The default
#' \eqn{k = 2} gives the classical Akaike information criterion.
#'
#' `BIC()` computes \eqn{-2 \ell + \log(n) p}, where \eqn{n} is the total
#' number of choice occasions, that is, `sum(object$choice_data$Tp)`.
#'
#' `coef()` returns the regression coefficients \eqn{b}. Coefficient names are
#' derived from `object$mod$var_names`.
#'
#' `confint()` constructs Wald confidence intervals for \eqn{b}.
#'
#' `logLik()` returns the fitted log-likelihood value stored in `object$ll_fit`.
#' The `"df"` attribute is the number of free model parameters,
#' `object$mod$length_theta`, and the `"nobs"` attribute is the total number of
#' choice occasions, `sum(object$choice_data$Tp)`.
#'
#' `plot()` visualizes the fitted population distribution of one or two
#' random effect coefficients selected by `parm`.
#'
#' `predict()` computes marginal choice probabilities for each choice occasion.
#'
#' `residuals()` computes response residuals
#' \eqn{r_{ij} = 1(y_i = j) - \hat{p}_{ij}} for every choice occasion \eqn{i}
#' and alternative \eqn{j}, where \eqn{\hat{p}_{ij}} is the fitted choice
#' probability.
#'
#' `summary()` summarizes the model results.
#'
#' `vcov()` returns the covariance matrix of the regression coefficients
#' \eqn{b}, matching the parameters returned by `coef()`.
#'
#' @return
#' `AIC()` and `BIC()` return a numeric scalar.
#'
#' `coef()` returns a named numeric vector of regression coefficients.
#'
#' `confint()` returns a numeric matrix with one row per coefficient
#' and two columns containing the lower and upper confidence limits.
#'
#' `logLik()` returns the fitted log-likelihood as a numeric scalar of class
#' `"logLik"` with `"df"` and `"nobs"` attributes.
#'
#' `plot()` returns a [`ggplot2::ggplot`] object.
#'
#' For `predict()`, the return value depends on `type`:
#' - `"probabilities"` returns a `data.frame` with one row per choice occasion.
#'   Its columns contain the predicted probabilities for the alternatives in
#'   `object$mod$alt_names` order.
#' - `"choice"` returns a numeric vector containing the predicted probability
#'   of the observed choice for each choice occasion.
#' - `"frequencies"` returns a `table` of predicted choices.
#' - `"confusion"` returns a two-dimensional `table` with observed choices in
#'   rows and predicted choices in columns.
#'
#' `residuals()` returns a numeric matrix with one row per choice occasion and
#' one column per alternative in `object$mod$alt_names` order.
#'
#' `summary()` returns an object of class `"summary.Rprobit"`.
#'
#' `vcov()` returns a named numeric covariance matrix with one row and column
#' per regression coefficient returned by `coef()`.
#'
#' @name Rprobit-methods

NULL

#' @rdname Rprobit-methods
#' @exportS3Method

AIC.Rprobit <- function(object, ..., k = 2) {
  checkmate::assert_r6(object, classes = "Rprobit")
  checkmate::assert_number(k, lower = 0, finite = TRUE)
  ll <- logLik(object)
  -2 * as.numeric(ll) + k * attr(ll, "df")
}

#' @rdname Rprobit-methods
#' @exportS3Method

BIC.Rprobit <- function(object, ...) {
  checkmate::assert_r6(object, classes = "Rprobit")
  ll <- logLik(object)
  nobs <- attr(ll, "nobs")
  if (nobs <= 0L) {
    stop("BIC() requires at least one choice occasion.", call. = FALSE)
  }
  -2 * as.numeric(ll) + base::log(nobs) * attr(ll, "df")
}

#' @rdname Rprobit-methods
#' @exportS3Method

coef.Rprobit <- function(object, ...) {
  checkmate::assert_r6(object, classes = "Rprobit")
  if (is.null(object$theta_fit)) {
    stop(
      "Field `theta_fit` is unavailable; run `$fit()` first.",
      call. = FALSE
    )
  }
  b <- object$mod$build_parameter(object$theta_fit)$b
  stats::setNames(as.numeric(b), rownames(b))
}

#' @rdname Rprobit-methods
#' @exportS3Method

confint.Rprobit <- function(object, parm, level = 0.95, ...) {
  checkmate::assert_r6(object, classes = "Rprobit")
  checkmate::assert_number(level, lower = 0, upper = 1, finite = TRUE)
  if (level <= 0 || level >= 1) {
    stop(
      "Argument `level` must be strictly between 0 and 1.",
      call. = FALSE
    )
  }
  coefficients <- coef(object)
  coefficient_names <- names(coefficients)
  if (missing(parm)) {
    parm <- coefficient_names
  } else if (is.numeric(parm)) {
    parm <- coefficient_names[parm]
  }
  standard_errors <- sqrt(pmax(diag(vcov(object)), 0))
  names(standard_errors) <- coefficient_names
  fixed <- rowSums(abs(object$mod$H_b)) == 0
  names(fixed) <- coefficient_names
  standard_errors[fixed] <- NA_real_
  alpha <- (1 - level) / 2
  probabilities <- c(alpha, 1 - alpha)
  quantiles <- stats::qnorm(probabilities)
  percentages <- paste0(
    format(
      100 * probabilities,
      digits = 3,
      trim = TRUE,
      scientific = FALSE
    ),
    " %"
  )
  intervals <- matrix(
    NA_real_,
    nrow = length(parm),
    ncol = 2L,
    dimnames = list(parm, percentages)
  )
  intervals[] <- coefficients[parm] + standard_errors[parm] %o% quantiles
  intervals
}

#' @rdname Rprobit-methods
#' @exportS3Method

logLik.Rprobit <- function(object, ...) {
  checkmate::assert_r6(object, classes = "Rprobit")
  if (is.null(object$ll_fit)) {
    stop(
      "Field `ll_fit` is unavailable; run `$fit()` first.",
      call. = FALSE
    )
  }
  checkmate::assert_number(object$ll_fit, finite = TRUE)
  structure(
    object$ll_fit,
    df = object$mod$length_theta,
    nobs = sum(object$choice_data$Tp),
    class = "logLik"
  )
}

#' @rdname Rprobit-methods
#' @exportS3Method

plot.Rprobit <- function(x, parm = NULL, cdf = FALSE, margin = 0.1, ...) {
  checkmate::assert_r6(x, classes = "Rprobit")
  checkmate::assert_flag(cdf)
  checkmate::assert_number(margin, lower = 0, finite = TRUE)
  margin <- max(margin, 0.1)

  coefficients <- coef(x)
  coefficient_names <- names(coefficients)
  random_names <- x$mod$var_names$random
  if (length(random_names) == 0L) {
    stop(
      "Method `plot()` requires at least one random coefficient.",
      call. = FALSE
    )
  }
  if (is.null(parm)) {
    parm <- random_names[1L]
  }
  if (is.character(parm)) {
    checkmate::assert_character(
      parm,
      min.len = 1L,
      max.len = 2L,
      any.missing = FALSE,
      unique = TRUE
    )
    indices <- match(parm, coefficient_names)
    if (anyNA(indices)) {
      stop(
        "Argument `parm` contains unknown coefficient name(s): ",
        paste(
          sprintf("`%s`", parm[is.na(indices)]),
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }
  } else {
    checkmate::assert_integerish(
      parm,
      lower = 1L,
      upper = length(coefficients),
      min.len = 1L,
      max.len = 2L,
      any.missing = FALSE,
      unique = TRUE
    )
    indices <- as.integer(parm)
  }
  selected_names <- coefficient_names[indices]
  non_random <- setdiff(selected_names, random_names)
  if (length(non_random) > 0L) {
    stop(
      "Argument `parm` may only select random coefficients; non-random ",
      "coefficient(s): ", paste(sprintf("`%s`", non_random), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  parameters <- x$mod$build_parameter(x$theta_fit)
  means <- as.numeric(coefficients[indices])
  covariance <- parameters$Omega[
    selected_names,
    selected_names,
    drop = FALSE
  ]
  tolerance <- sqrt(.Machine$double.eps)
  eigenvalues <- eigen(
    covariance,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (any(eigenvalues <= tolerance)) {
    stop(
      paste(
        "The covariance matrix for the selected random coefficient(s)",
        "must be positive definite."
      ),
      call. = FALSE
    )
  }

  if (length(indices) == 1L) {
    standard_deviation <- sqrt(covariance[1L, 1L])
    half_width <- 4 * standard_deviation
    limits <- means + c(-1, 1) * half_width * (1 + margin)
    values <- seq(limits[1L], limits[2L], length.out = 201L)
    distribution <- if (isTRUE(cdf)) {
      stats::pnorm(values, mean = means, sd = standard_deviation)
    } else {
      stats::dnorm(values, mean = means, sd = standard_deviation)
    }
    plot_data <- data.frame(
      value = values,
      distribution = distribution
    )
    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = !!as.name("value"),
        y = !!as.name("distribution")
      )
    )
    plot <- plot + ggplot2::geom_line(linewidth = 0.8)
    plot <- plot + ggplot2::geom_vline(
      xintercept = means,
      linetype = "dashed"
    )
    if (isTRUE(cdf)) {
      plot <- plot + ggplot2::coord_cartesian(ylim = c(0, 1))
    }
    return(plot + ggplot2::labs(
      x = selected_names,
      y = if (isTRUE(cdf)) "Cumulative probability" else "Density",
      title = paste(
        if (isTRUE(cdf)) "CDF of" else "Density of",
        selected_names
      )
    ) + ggplot2::theme_minimal())
  }

  standard_deviations <- sqrt(diag(covariance))
  limits <- lapply(seq_len(2L), function(j) {
    means[j] + c(-1, 1) * 4 * standard_deviations[j] * (1 + margin)
  })
  grid <- expand.grid(
    value_1 = seq(limits[[1L]][1L], limits[[1L]][2L], length.out = 51L),
    value_2 = seq(limits[[2L]][1L], limits[[2L]][2L], length.out = 51L)
  )
  grid$distribution <- vapply(seq_len(nrow(grid)), function(row) {
    value <- as.numeric(grid[row, c("value_1", "value_2")])
    if (isTRUE(cdf)) {
      oeli::pmvnorm(x = value, mean = means, Sigma = covariance)
    } else {
      oeli::dmvnorm(x = value, mean = means, Sigma = covariance)
    }
  }, numeric(1))
  mean_data <- data.frame(value_1 = means[1L], value_2 = means[2L])
  ggplot2::ggplot(
    grid,
    ggplot2::aes(
      x = !!as.name("value_1"),
      y = !!as.name("value_2"),
      z = !!as.name("distribution")
    )
  ) +
    ggplot2::geom_contour() +
    ggplot2::geom_point(
      data = mean_data,
      mapping = ggplot2::aes(
        x = !!as.name("value_1"),
        y = !!as.name("value_2")
      ),
      inherit.aes = FALSE,
      size = 2
    ) +
    ggplot2::labs(
      x = selected_names[1L],
      y = selected_names[2L],
      title = if (isTRUE(cdf)) {
        "Joint coefficient CDF"
      } else {
        "Joint coefficient density"
      }
    ) +
    ggplot2::theme_minimal()
}

#' @rdname Rprobit-methods
#' @exportS3Method

predict.Rprobit <- function(
    object,
    newdata = NULL,
    type = c("probabilities", "choice", "frequencies", "confusion"),
    ...
) {
  checkmate::assert_r6(object, classes = "Rprobit")
  type <- match.arg(type)
  if (is.null(object$theta_fit)) {
    stop(
      "Field `theta_fit` is unavailable; run `$fit()` first.",
      call. = FALSE
    )
  }
  checkmate::assert_numeric(
    object$theta_fit,
    len = object$mod$length_theta,
    any.missing = FALSE,
    finite = TRUE
  )
  prediction_data <- if (is.null(newdata)) {
    object$choice_data
  } else if (inherits(newdata, "choice_data")) {
    newdata
  } else {
    stop(
      "Argument `newdata` must be `NULL` or a `choice_data` object.",
      call. = FALSE
    )
  }
  design <- prediction_data$build_design(
    mod = object$mod,
    difference = FALSE,
    check_choice = type %in% c("choice", "confusion")
  )
  probabilities <- if (length(design) == 0L) {
    matrix(
      numeric(),
      nrow = 0L,
      ncol = object$mod$alt,
      dimnames = list(NULL, object$mod$alt_names)
    )
  } else {
    do.call(rbind, lapply(design, function(design_single) {
      pred_probit_approx(
        theta = object$theta_fit,
        design_single = design_single,
        mod_list = object$mod$mod_list
      )
    }))
  }
  colnames(probabilities) <- object$mod$alt_names
  choices <- unlist(
    lapply(design, function(design_single) design_single$y),
    use.names = FALSE
  )
  out <- as.data.frame(probabilities)
  names(out) <- object$mod$alt_names
  choice_prob <- if (length(choices) == 0L) {
    numeric()
  } else {
    probabilities[cbind(seq_along(choices), choices)]
  }
  if (type == "probabilities") {
    return(out)
  }
  if (type == "choice") {
    return(choice_prob)
  }
  predicted_index <- max.col(probabilities, ties.method = "first")
  predicted_choices <- factor(
    object$mod$alt_names[predicted_index],
    levels = object$mod$alt_names
  )
  if (type == "frequencies") {
    return(table(predicted = predicted_choices))
  }
  observed_choices <- factor(
    object$mod$alt_names[choices],
    levels = object$mod$alt_names
  )
  table(observed = observed_choices, predicted = predicted_choices)
}

#' @rdname Rprobit-methods
#' @exportS3Method

residuals.Rprobit <- function(object, ...) {
  checkmate::assert_r6(object, classes = "Rprobit")
  probabilities <- as.matrix(
    predict(object, type = "probabilities")
  )
  design <- object$choice_data$build_design(
    mod = object$mod,
    difference = FALSE
  )
  choices <- unlist(
    lapply(design, function(design_single) design_single$y),
    use.names = FALSE
  )
  residuals <- -probabilities
  if (length(choices) > 0L) {
    indices <- cbind(seq_along(choices), choices)
    residuals[indices] <- residuals[indices] + 1
  }
  residuals
}

#' @rdname Rprobit-methods
#' @exportS3Method

summary.Rprobit <- function(object, ...) {
  checkmate::assert_r6(object, classes = "Rprobit")
  estimates <- coef(object)
  covariance <- vcov(object)
  standard_errors <- sqrt(pmax(diag(covariance), 0))
  fixed <- rowSums(abs(object$mod$H_b)) == 0
  standard_errors[fixed] <- NA_real_
  statistics <- estimates / standard_errors
  p_values <- 2 * stats::pnorm(abs(statistics), lower.tail = FALSE)
  coefficients <- cbind(
    "Estimate" = estimates,
    "Std. Error" = standard_errors,
    "z value" = statistics,
    "Pr(>|z|)" = p_values
  )
  parameters <- object$mod$build_parameter(object$theta_fit)
  covariance_components <- list(
    "Omega" = list(
      chol = parameters$Omega_chol,
      H = object$mod$H_O,
      theta_index = object$mod$length_theta_b +
        seq_len(object$mod$length_theta_O)
    ),
    "Sigma" = list(
      chol = parameters$Sigma_chol,
      H = object$mod$H_L,
      theta_index = object$mod$length_theta_b +
        object$mod$length_theta_O +
        seq_len(object$mod$length_theta_L)
    )
  )
  covariance_summaries <- vector(
    mode = "list",
    length = length(covariance_components)
  )
  names(covariance_summaries) <- names(covariance_components)
  for (label in names(covariance_components)) {
    component <- covariance_components[[label]]
    chol <- component$chol
    lower <- lower.tri(chol, diag = TRUE)
    if (!any(lower)) {
      covariance_summaries[[label]] <- matrix(
        numeric(),
        nrow = 0L,
        ncol = ncol(coefficients),
        dimnames = list(character(), colnames(coefficients))
      )
      next
    }

    covariance_estimates <- (chol %*% t(chol))[lower]
    positions <- which(lower, arr.ind = TRUE)
    parameter_names <- paste0(
      label,
      "[",
      rownames(chol)[positions[, 1L]],
      ",",
      rownames(chol)[positions[, 2L]],
      "]"
    )
    jacobian <- matrix(
      0,
      nrow = length(covariance_estimates),
      ncol = ncol(component$H)
    )
    if (ncol(jacobian) > 0L) {
      for (j in seq_len(ncol(jacobian))) {
        chol_derivative <- matrix(0, nrow = nrow(chol), ncol = ncol(chol))
        chol_derivative[lower] <- component$H[, j]
        covariance_derivative <-
          chol_derivative %*% t(chol) + chol %*% t(chol_derivative)
        jacobian[, j] <- covariance_derivative[lower]
      }
    }
    covariance_standard_errors <- rep(
      NA_real_,
      length(covariance_estimates)
    )
    if (ncol(jacobian) > 0L) {
      covariance_vcov <- object$theta_vcov[
        component$theta_index,
        component$theta_index,
        drop = FALSE
      ]
      covariance_variances <- rowSums(
        (jacobian %*% covariance_vcov) * jacobian
      )
      covariance_standard_errors <- sqrt(pmax(covariance_variances, 0))
      covariance_standard_errors[rowSums(abs(jacobian)) == 0] <- NA_real_
    }
    covariance_statistics <-
      covariance_estimates / covariance_standard_errors
    covariance_p_values <- 2 * stats::pnorm(
      abs(covariance_statistics),
      lower.tail = FALSE
    )
    covariance_summaries[[label]] <- cbind(
      "Estimate" = covariance_estimates,
      "Std. Error" = covariance_standard_errors,
      "z value" = covariance_statistics,
      "Pr(>|z|)" = covariance_p_values
    )
    rownames(covariance_summaries[[label]]) <- parameter_names
  }
  coefficients <- rbind(
    coefficients,
    covariance_summaries$Omega,
    covariance_summaries$Sigma
  )

  predictions <- predict(object, type = "probabilities")
  observed_frequencies <- stats::setNames(
    numeric(object$mod$alt),
    object$mod$alt_names
  )
  choice_frequencies <- object$choice_data$choice_frequencies
  observed_frequencies[names(choice_frequencies)] <-
    as.numeric(choice_frequencies)
  prediction_summary <- rbind(
    "Observed choice frequency" = observed_frequencies,
    "Mean predicted probability" = colMeans(predictions)
  )
  likelihood <- logLik(object)

  structure(
    list(
      formula = object$mod$formula,
      coefficients = coefficients,
      covariance = covariance,
      parameters = parameters,
      residuals = residuals(object),
      prediction_summary = prediction_summary,
      predictions = predictions,
      logLik = likelihood,
      AIC = AIC(object),
      BIC = BIC(object),
      nobs = attr(likelihood, "nobs"),
      N = object$choice_data$N,
      Tp = object$choice_data$Tp,
      alternatives = object$mod$alt_names
    ),
    class = c("summary.Rprobit", "list")
  )
}

#' @rdname Rprobit-methods
#' @exportS3Method

print.summary.Rprobit <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    signif.stars = getOption("show.signif.stars"),
    ...
) {
  checkmate::assert_class(x, "summary.Rprobit")
  checkmate::assert_number(digits, lower = 1, finite = TRUE)
  checkmate::assert_flag(signif.stars)
  digits <- as.integer(digits)

  cat("\nRprobit model summary\n")
  cat("\nFormula:\n")
  cat(paste(deparse(x$formula), collapse = "\n"), "\n")

  cat("\nData:\n")
  cat("  Deciders:", x$N, "\n")
  cat("  Choice occasions:", x$nobs, "\n")
  if (length(x$Tp) > 0L) {
    if (length(unique(x$Tp)) == 1L) {
      cat("  Choice occasions per decider:", unique(x$Tp), "\n")
    } else {
      cat(
        "  Choice occasions per decider:",
        min(x$Tp), "to", max(x$Tp), "\n"
      )
    }
  }
  cat(sprintf(
    "  Alternatives: %d (%s)\n",
    length(x$alternatives),
    paste(x$alternatives, collapse = ", ")
  ))

  cat("\nResponse residuals:\n")
  if (nrow(x$residuals) > 5L) {
    residual_summary <- apply(
      x$residuals,
      2L,
      stats::quantile,
      probs = c(0, 0.25, 0.5, 0.75, 1),
      names = FALSE
    )
    rownames(residual_summary) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(residual_summary, digits = digits)
  } else if (nrow(x$residuals) > 0L) {
    print(x$residuals, digits = digits)
  } else {
    cat("  No residuals.\n")
  }

  cat("\nCoefficients:\n")
  stats::printCoefmat(
    x$coefficients,
    digits = digits,
    signif.stars = signif.stars,
    na.print = "NA",
    ...
  )

  if (length(x$parameters$Omega) > 0L) {
    cat("\nRandom-effect covariance matrix (Omega):\n")
    print(x$parameters$Omega, digits = digits)
  }
  cat("\nError covariance matrix (Sigma):\n")
  print(x$parameters$Sigma, digits = digits)

  cat("\nPrediction summary:\n")
  print(x$prediction_summary, digits = digits)

  cat("\nFit statistics:\n")
  cat(
    "  Log-likelihood:",
    formatC(as.numeric(x$logLik), digits = digits, format = "g"),
    "\n"
  )
  cat("  AIC:", formatC(x$AIC, digits = digits, format = "g"), "\n")
  cat("  BIC:", formatC(x$BIC, digits = digits, format = "g"), "\n\n")
  invisible(x)
}

#' @rdname Rprobit-methods
#' @exportS3Method

vcov.Rprobit <- function(object, ...) {
  checkmate::assert_r6(object, classes = "Rprobit")
  if (is.null(object$theta_vcov)) {
    stop(
      "Field `theta_vcov` is unavailable; run `$fit()` first.",
      call. = FALSE
    )
  }
  coefficient_names <- names(coef(object))
  theta_b_index <- seq_len(object$mod$length_theta_b)
  theta_b_vcov <- object$theta_vcov[
    theta_b_index,
    theta_b_index,
    drop = FALSE
  ]
  coefficient_vcov <- object$mod$H_b %*%
    theta_b_vcov %*%
    t(object$mod$H_b)
  coefficient_vcov <- (coefficient_vcov + t(coefficient_vcov)) / 2
  dimnames(coefficient_vcov) <- list(coefficient_names, coefficient_names)
  coefficient_vcov
}
