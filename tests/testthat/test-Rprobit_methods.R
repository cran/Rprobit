test_that("AIC is computed from the fitted log-likelihood", {
  m <- mod$new(
    formula = choice ~ x | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("A", "B"),
    error_struc = "fixed"
  )
  object <- Rprobit$new(m, choice_data$new())

  expect_error(
    AIC(object),
    "Field `ll_fit` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  object$.__enclos_env__$private$.theta_fit <- 0.5
  object$.__enclos_env__$private$.ll_fit <- -10

  expect_equal(AIC(object), 22)
  expect_equal(AIC(object, k = 3), 23)
  expect_error(AIC(object, k = -1), "not >= 0", fixed = TRUE)
})

test_that("BIC uses the total number of choice occasions", {
  m <- mod$new(
    formula = choice ~ x | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("A", "B"),
    error_struc = "fixed"
  )
  empty_object <- Rprobit$new(m, choice_data$new())
  empty_object$.__enclos_env__$private$.theta_fit <- 0.5
  empty_object$.__enclos_env__$private$.ll_fit <- -10

  expect_error(
    BIC(empty_object),
    "BIC() requires at least one choice occasion.",
    fixed = TRUE
  )

  d <- choice_data$new(data.frame(choice = c("A", "B", "A")))
  object <- Rprobit$new(m, d)
  expect_error(
    BIC(object),
    "Field `ll_fit` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  object$.__enclos_env__$private$.theta_fit <- 0.5
  object$.__enclos_env__$private$.ll_fit <- -10

  expect_equal(BIC(object), 20 + log(3))
})

test_that("coef returns named regression coefficients", {
  m <- mod$new(
    formula = choice ~ x | z + 0 | 0,
    re = "z",
    alt = 3,
    alt_names = c("A", "B", "C"),
    error_struc = "fixed"
  )
  object <- Rprobit$new(m, choice_data$new())

  expect_error(
    coef(object),
    "Field `theta_fit` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  object$.__enclos_env__$private$.theta_fit <- c(0.5, -1, 2, 0.2, 0.3)
  coefficients <- coef(object)

  expect_identical(coefficients, c(x = 0.5, z_B = -1, z_C = 2))
  expect_type(coefficients, "double")
})

test_that("confint computes coefficient confidence intervals", {
  H_b <- matrix(
    c(
      1, 0,
      0, 1,
      1, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  m <- mod$new(
    formula = choice ~ x | z + 0 | 0,
    re = character(),
    alt = 3,
    alt_names = c("A", "B", "C"),
    error_struc = "fixed",
    H_b = H_b,
    f_b = matrix(0, nrow = 3, ncol = 1)
  )
  object <- Rprobit$new(m, choice_data$new())
  object$.__enclos_env__$private$.theta_fit <- c(1, 2)

  expect_error(
    confint(object),
    "Field `theta_vcov` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  theta_vcov <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)
  object$.__enclos_env__$private$.theta_vcov <- theta_vcov
  estimates <- c(x = 1, z_B = 2, z_C = 3)
  standard_errors <- sqrt(diag(H_b %*% theta_vcov %*% t(H_b)))
  quantiles <- stats::qnorm(c(0.025, 0.975))
  expected <- estimates + standard_errors %o% quantiles
  rownames(expected) <- names(estimates)
  colnames(expected) <- c("2.5 %", "97.5 %")

  expect_equal(confint(object), expected)

  probabilities_90 <- c(0.05, 0.95)
  expected_90 <- estimates +
    standard_errors %o% stats::qnorm(probabilities_90)
  rownames(expected_90) <- names(estimates)
  colnames(expected_90) <- c("5 %", "95 %")
  expect_equal(
    confint(object, parm = c("z_C", "x"), level = 0.9),
    expected_90[c("z_C", "x"), , drop = FALSE]
  )
  expect_equal(
    confint(object, parm = c(3, 1), level = 0.9),
    expected_90[c("z_C", "x"), , drop = FALSE]
  )
  expect_error(
    confint(object, level = 1),
    "`level` must be strictly between 0 and 1.",
    fixed = TRUE
  )
})

test_that("logLik returns the fitted value with model metadata", {
  m <- mod$new(
    formula = choice ~ x | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("A", "B"),
    error_struc = "fixed"
  )
  d <- choice_data$new(data.frame(choice = c("A", "B", "A")))
  object <- Rprobit$new(m, d)

  expect_error(
    logLik(object),
    "Field `ll_fit` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  object$.__enclos_env__$private$.ll_fit <- -10.5
  likelihood <- logLik(object)

  expect_s3_class(likelihood, "logLik")
  expect_identical(as.numeric(likelihood), -10.5)
  expect_identical(attr(likelihood, "df"), m$length_theta)
  expect_identical(attr(likelihood, "nobs"), sum(d$Tp))
})

test_that("plot returns ggplot coefficient distributions", {
  m <- mod$new(
    formula = choice ~ x + z | 0 | 0,
    re = c("x", "z"),
    alt = 2,
    alt_names = c("A", "B"),
    error_struc = "fixed"
  )
  object <- Rprobit$new(m, choice_data$new())

  expect_error(
    plot(object),
    "Field `theta_fit` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  object$.__enclos_env__$private$.theta_fit <- c(0.5, -0.5, 1, 0.5)

  density_plot <- plot(object, parm = "x")
  expect_s3_class(density_plot, "ggplot")
  expect_identical(density_plot$labels$x, "x")
  expect_identical(density_plot$labels$y, "Density")
  expect_identical(nrow(density_plot$data), 201L)
  expect_true(all(density_plot$data$distribution >= 0))

  cdf_plot <- plot(object, parm = 1, cdf = TRUE)
  expect_s3_class(cdf_plot, "ggplot")
  expect_identical(cdf_plot$labels$y, "Cumulative probability")
  expect_true(all(diff(cdf_plot$data$distribution) >= 0))

  joint_plot <- plot(object, parm = c("x", "z"))
  expect_s3_class(joint_plot, "ggplot")
  expect_identical(joint_plot$labels$x, "x")
  expect_identical(joint_plot$labels$y, "z")
  expect_equal(nrow(joint_plot$data), 51L^2)
  expect_true(all(joint_plot$data$distribution >= 0))

  degenerate_object <- Rprobit$new(m, choice_data$new())
  degenerate_object$.__enclos_env__$private$.theta_fit <-
    c(0.5, -0.5, 0, 0)
  covariance_error <- paste(
    "The covariance matrix for the selected random coefficient(s)",
    "must be positive definite."
  )
  expect_error(
    plot(degenerate_object, parm = "x"),
    covariance_error,
    fixed = TRUE
  )
  expect_error(
    plot(degenerate_object, parm = c("x", "z"), cdf = TRUE),
    covariance_error,
    fixed = TRUE
  )

  mixed_mod <- mod$new(
    formula = choice ~ x + z | 0 | 0,
    re = "z",
    alt = 2,
    alt_names = c("A", "B"),
    error_struc = "fixed"
  )
  mixed_object <- Rprobit$new(mixed_mod, choice_data$new())
  mixed_object$.__enclos_env__$private$.theta_fit <- c(0.5, -0.5, 0.5)
  default_plot <- plot(mixed_object)
  expect_identical(default_plot$labels$x, "z")
  expect_error(
    plot(mixed_object, parm = "x"),
    paste(
      "Argument `parm` may only select random coefficients; non-random",
      "coefficient(s): `x`."
    ),
    fixed = TRUE
  )

  expect_error(
    plot(object, parm = "unknown"),
    "`parm` contains unknown coefficient name(s): `unknown`.",
    fixed = TRUE
  )
  expect_error(plot(object, parm = 3), "not <= 2", fixed = TRUE)

  fixed_mod <- mod$new(
    formula = choice ~ x | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("A", "B"),
    error_struc = "fixed"
  )
  fixed_object <- Rprobit$new(fixed_mod, choice_data$new())
  fixed_object$.__enclos_env__$private$.theta_fit <- 0.5
  expect_error(
    plot(fixed_object),
    "Method `plot()` requires at least one random coefficient.",
    fixed = TRUE
  )
})

test_that("predict returns choice probabilities from the fitted model", {
  m <- mod$new(
    formula = choice ~ x | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("A", "B"),
    error_struc = "fixed"
  )
  df <- data.frame(
    id = c(1, 1),
    choice = c("A", "B"),
    x_A = c(1, 0),
    x_B = c(0, 1)
  )
  d <- choice_data$new(df = df, column_decider_id = "id")
  object <- Rprobit$new(m, d)

  expect_error(
    predict(object),
    "Field `theta_fit` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  object$.__enclos_env__$private$.theta_fit <- 1
  predictions <- predict(object)
  probability_A <- stats::pnorm(c(1, -1) / sqrt(2))
  expected_probabilities <- cbind(
    A = probability_A,
    B = 1 - probability_A
  )

  expect_s3_class(predictions, "data.frame")
  expect_identical(names(predictions), c("A", "B"))
  expect_equal(
    as.matrix(predictions),
    expected_probabilities,
    tolerance = 1e-8
  )
  expect_equal(rowSums(predictions), rep(1, 2))
  expect_equal(
    predict(object, type = "choice"),
    expected_probabilities[cbind(seq_len(nrow(df)), c(1, 2))]
  )

  new_data <- choice_data$new(
    df = df[2:1, ],
    column_decider_id = "id"
  )
  new_predictions <- predict(object, newdata = new_data)
  expect_equal(
    as.matrix(new_predictions[c("A", "B")]),
    expected_probabilities[2:1, ],
    tolerance = 1e-8
  )

  frequencies <- predict(object, type = "frequencies")
  expect_s3_class(frequencies, "table")
  expect_identical(
    frequencies,
    structure(c(A = 1L, B = 1L), dim = 2L, dimnames = list(
      predicted = c("A", "B")
    ), class = "table")
  )

  confusion_data <- df
  confusion_data$choice <- c("B", "A")
  confusion_data <- choice_data$new(
    df = confusion_data,
    column_decider_id = "id"
  )
  confusion <- predict(
    object,
    newdata = confusion_data,
    type = "confusion"
  )
  expect_s3_class(confusion, "table")
  expect_identical(
    confusion,
    structure(
      c(0L, 1L, 1L, 0L),
      dim = c(2L, 2L),
      dimnames = list(
        observed = c("A", "B"),
        predicted = c("A", "B")
      ),
      class = "table"
    )
  )

  expect_error(
    predict(object, newdata = df),
    "`newdata` must be `NULL` or a `choice_data` object.",
    fixed = TRUE
  )
})

test_that("residuals return observed indicators minus fitted probabilities", {
  m <- mod$new(
    formula = choice ~ x | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("A", "B"),
    error_struc = "fixed"
  )
  df <- data.frame(
    id = c(1, 1),
    choice = c("A", "B"),
    x_A = c(1, 0),
    x_B = c(0, 1)
  )
  d <- choice_data$new(df = df, column_decider_id = "id")
  object <- Rprobit$new(m, d)

  expect_error(
    residuals(object),
    "Field `theta_fit` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  object$.__enclos_env__$private$.theta_fit <- 1
  probability_A <- stats::pnorm(c(1, -1) / sqrt(2))
  probabilities <- cbind(
    A = probability_A,
    B = 1 - probability_A
  )
  observed <- rbind(
    c(A = 1, B = 0),
    c(A = 0, B = 1)
  )
  expected <- observed - probabilities
  residual_values <- residuals(object)

  expect_type(residual_values, "double")
  expect_equal(residual_values, expected, tolerance = 1e-8)
  expect_equal(rowSums(residual_values), c(0, 0), tolerance = 1e-8)
})

test_that("summary presents model results", {
  m <- mod$new(
    formula = choice ~ x | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("A", "B"),
    error_struc = "fixed"
  )
  df <- data.frame(
    id = c(1, 1),
    choice = c("A", "B"),
    x_A = c(1, 0),
    x_B = c(0, 1)
  )
  d <- choice_data$new(df = df, column_decider_id = "id")
  object <- Rprobit$new(m, d)
  object$.__enclos_env__$private$.theta_fit <- 1
  object$.__enclos_env__$private$.ll_fit <- -10

  expect_error(
    summary(object),
    "Field `theta_vcov` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  object$.__enclos_env__$private$.theta_vcov <- matrix(0.04)
  result <- summary(object)
  expected_coefficients <- matrix(
    c(1, 0.2, 5, 2 * stats::pnorm(5, lower.tail = FALSE)),
    nrow = 1,
    dimnames = list(
      "x",
      c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    )
  )
  probability_A <- stats::pnorm(c(1, -1) / sqrt(2))
  expected_predictions <- data.frame(
    A = probability_A,
    B = 1 - probability_A
  )
  expected_prediction_summary <- rbind(
    "Observed choice frequency" = c(A = 0.5, B = 0.5),
    "Mean predicted probability" = c(A = 0.5, B = 0.5)
  )

  expect_s3_class(result, "summary.Rprobit")
  expect_named(
    result,
    c(
      "formula", "coefficients", "covariance", "parameters", "residuals",
      "prediction_summary", "predictions", "logLik", "AIC", "BIC", "nobs",
      "N", "Tp", "alternatives"
    )
  )
  expect_equal(
    result$coefficients["x", , drop = FALSE],
    expected_coefficients
  )
  expect_equal(
    result$coefficients[
      c("Sigma[A,A]", "Sigma[B,A]", "Sigma[B,B]"),
      "Estimate"
    ],
    c("Sigma[A,A]" = 1, "Sigma[B,A]" = 0, "Sigma[B,B]" = 1)
  )
  expect_true(all(is.na(
    result$coefficients[
      c("Sigma[A,A]", "Sigma[B,A]", "Sigma[B,B]"),
      -1L
    ]
  )))
  expect_equal(result$covariance, matrix(0.04, dimnames = rep(list("x"), 2)))
  expect_equal(result$predictions, expected_predictions, tolerance = 1e-8)
  expect_equal(
    result$prediction_summary,
    expected_prediction_summary,
    tolerance = 1e-8
  )
  expect_equal(result$residuals, residuals(object))
  expect_identical(result$logLik, logLik(object))
  expect_equal(result$AIC, AIC(object))
  expect_equal(result$BIC, BIC(object))
  expect_identical(result$nobs, 2L)
  expect_identical(result$N, 1L)
  expect_identical(result$Tp, 2L)
  expect_identical(result$alternatives, c("A", "B"))
})

test_that("vcov returns the covariance matrix of the coefficients", {
  H_b <- matrix(
    c(
      1, 0,
      0, 1,
      1, 1
    ),
    nrow = 3,
    byrow = TRUE
  )
  m <- mod$new(
    formula = choice ~ x | z + 0 | 0,
    re = character(),
    alt = 3,
    alt_names = c("A", "B", "C"),
    error_struc = "fixed",
    H_b = H_b,
    f_b = matrix(0, nrow = 3, ncol = 1)
  )
  object <- Rprobit$new(m, choice_data$new())
  object$.__enclos_env__$private$.theta_fit <- c(1, 2)

  expect_error(
    vcov(object),
    "Field `theta_vcov` is unavailable; run `$fit()` first.",
    fixed = TRUE
  )

  theta_vcov <- matrix(c(0.04, 0.01, 0.01, 0.09), nrow = 2)
  object$.__enclos_env__$private$.theta_vcov <- theta_vcov
  expected <- H_b %*% theta_vcov %*% t(H_b)
  dimnames(expected) <- rep(list(c("x", "z_B", "z_C")), 2)

  expect_equal(vcov(object), expected)
  expect_identical(dimnames(vcov(object)), dimnames(expected))
})
