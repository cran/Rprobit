choice_data_wide <- data.frame(
  id = c(1, 1, 2),
  choice = c("bus", "car", "bus"),
  price_bus = c(1, 2, 3),
  price_car = c(4, 1, 0),
  income = c(10, 20, 30),
  time_bus = c(5, 6, 7),
  time_car = c(8, 9, 10)
)

choice_data_long <- data.frame(
  id = rep(c(1, 1, 2), each = 2),
  occasion = rep(c(1, 2, 1), each = 2),
  mode = rep(c("bus", "car"), 3),
  choice = c(1, 0, 0, 1, 1, 0),
  price = c(1, 4, 2, 1, 3, 0),
  income = rep(c(10, 20, 30), each = 2),
  time = c(5, 8, 6, 9, 7, 10)
)

test_that("choice_data object can be defined", {

  ### empty data
  empty <- choice_data$new()
  expect_identical(empty$alt_names, character())
  expect_equal(empty$choice_counts, as.table(numeric()))
  expect_equal(empty$choice_frequencies, as.table(numeric()))

  ### wide data
  d <- choice_data$new(
    df = choice_data_wide,
    column_decider_id = "id"
  )
  expect_s3_class(d, "choice_data")
  expected_wide <- cbind(choice_data_wide, occasion = c(1L, 2L, 1L))
  expect_identical(d$df, expected_wide)
  expect_identical(d$format, "wide")
  expect_identical(d$alt_names, c("bus", "car"))
  expect_null(d$theta_dgp)
  expect_identical(d$N, 2L)
  expect_identical(d$Tp, c(2L, 1L))
  expect_equal(d$choice_counts, as.table(c(bus = 2, car = 1)))
  expect_equal(
    d$choice_frequencies,
    as.table(c(bus = 2 / 3, car = 1 / 3))
  )
  expect_snapshot(print(d))

  df_subclass <- choice_data_wide
  class(df_subclass) <- c("choice_data_subclass", class(df_subclass))
  coerced <- choice_data$new(df = df_subclass)
  expect_identical(class(coerced$df), "data.frame")

  ### wide data with occasion identifiers
  d_occasion <- choice_data$new(
    df = cbind(choice_data_wide, occasion = c(1, 2, 1)),
    column_decider_id = "id",
    column_occasion_id = "occasion"
  )
  expect_identical(d_occasion$N, 2L)
  expect_identical(d_occasion$Tp, c(2L, 1L))

  ### automatic occasion identifiers are added to wide and long data
  automatic_occasion <- choice_data$new(
    df = choice_data_wide,
    column_decider_id = "id"
  )
  expect_identical(automatic_occasion$df$occasion, c(1L, 2L, 1L))
  long_without_occasion <- choice_data_long
  long_without_occasion$occasion <- NULL
  long_no_occasion <- choice_data$new(
    df = long_without_occasion,
    format = "long",
    column_decider_id = "id",
    column_alternative = "mode"
  )
  expect_identical(
    long_no_occasion$df$occasion,
    c(1L, 1L, 2L, 2L, 1L, 1L)
  )

  ### long data
  long <- choice_data$new(
    df = choice_data_long,
    format = "long",
    column_decider_id = "id",
    column_occasion_id = "occasion",
    column_alternative = "mode"
  )
  expect_identical(long$format, "long")
  expect_identical(long$alt_names, c("bus", "car"))
  expect_identical(long$N, 2L)
  expect_identical(long$Tp, c(2L, 1L))
  expect_equal(long$choice_counts, d$choice_counts)
  expect_equal(long$choice_frequencies, d$choice_frequencies)
  expect_snapshot(print(long))

  ### validation of constructor inputs
  expect_error(
    choice_data$new(
      df = choice_data_wide,
      format = "long"
    ),
    "`column_alternative` must be specified",
    fixed = TRUE
  )
  expect_warning(
    choice_data$new(
      df = choice_data_wide,
      column_alternative = "mode"
    ),
    "`column_alternative` is only used",
    fixed = TRUE
  )

  expect_error(
    choice_data$new(
      df = choice_data_wide,
      column_occasion_id = "id"
    ),
    "`column_occasion_id` requires `column_decider_id`.",
    fixed = TRUE
  )

  expect_error(
    choice_data$new(
      df = choice_data_long,
      format = "long",
      column_occasion_id = "occasion",
      column_alternative = "mode"
    ),
    "`column_occasion_id` requires `column_decider_id`.",
    fixed = TRUE
  )

  duplicate_occasion <- cbind(choice_data_wide, occasion = c(1, 1, 1))
  expect_error(
    choice_data$new(
      df = duplicate_occasion,
      column_decider_id = "id",
      column_occasion_id = "occasion"
    ),
    "Wide data must have unique decider and occasion identifier rows.",
    fixed = TRUE
  )

  d <- choice_data$new(df = choice_data_wide)
  expect_error(
    d$choice_counts <- c(bus = 1),
    "Field `choice_counts` is read-only.",
    fixed = TRUE
  )
  expect_error(
    d$alt_names <- c("bus", "car"),
    "Field `alt_names` is read-only.",
    fixed = TRUE
  )
  expect_error(
    d$choice_frequencies <- c(bus = 1),
    "Field `choice_frequencies` is read-only.",
    fixed = TRUE
  )
  expect_error(
    d$theta_dgp <- numeric(),
    "Field `theta_dgp` is read-only.",
    fixed = TRUE
  )

  ### validation of long data
  bad_long <- choice_data_long
  bad_long$mode[2] <- "bus"
  expect_error(
    choice_data$new(
      df = bad_long,
      format = "long",
      column_decider_id = "id",
      column_occasion_id = "occasion",
      column_alternative = "mode"
    ),
    "Each long-data choice set must contain each alternative once.",
    fixed = TRUE
  )

  label_long <- choice_data_long
  label_long$choice <- rep(c("bus", "car", "bus"), each = 2)
  expect_error(
    choice_data$new(
      df = label_long,
      format = "long",
      column_decider_id = "id",
      column_occasion_id = "occasion",
      column_alternative = "mode"
    ),
    "Each long-data choice set must identify exactly one choice.",
    fixed = TRUE
  )

  missing_choice_long <- choice_data_long
  missing_choice_long$choice[2] <- NA
  expect_warning(
    long_miss <- choice_data$new(
      df = missing_choice_long,
      format = "long",
      column_decider_id = "id",
      column_occasion_id = "occasion",
      column_alternative = "mode"
    ),
    "Choice column contains missing values.",
    fixed = TRUE
  )
  expect_snapshot(print(long_miss))

  no_choice_long <- choice_data_long
  no_choice_long$choice <- NULL
  expect_warning(
    choice_data$new(
      df = no_choice_long,
      format = "long",
      column_decider_id = "id",
      column_occasion_id = "occasion",
      column_alternative = "mode"
    ),
    "Choice column contains missing values.",
    fixed = TRUE
  )

  long_no_id <- choice_data_long
  long_no_id$id <- NULL
  no_id <- choice_data$new(
    df = long_no_id,
    format = "long",
    column_alternative = "mode"
  )
  expect_identical(no_id$N, 3L)
  expect_identical(no_id$Tp, rep(1L, 3))
  expect_equal(no_id$choice_counts, as.table(c(bus = 2, car = 1)))

  ### validation of design inputs
  m <- mod$new(
    formula = choice ~ missing | 0 | 0,
    alt = 2,
    alt_names = c("bus", "car")
  )
  expect_error(
    d$build_design(m),
    "Required covariate column(s) are missing",
    fixed = TRUE
  )

  m <- mod$new(
    formula = choice ~ price | income | time,
    alt = 2,
    alt_names = c("bus", "car")
  )
  non_numeric <- choice_data_wide
  non_numeric$price_bus <- "not numeric"
  expect_error(
    choice_data$new(df = non_numeric)$build_design(m),
    "Covariate column(s) must be numeric: `price_bus`.",
    fixed = TRUE
  )
  missing_value <- choice_data_wide
  missing_value$price_car[1] <- NA
  expect_warning(
    missing_data <- choice_data$new(df = missing_value),
    "Covariate column(s) contain missing values: `price_car`.",
    fixed = TRUE
  )
  expect_no_warning(missing_data$build_design(m))
})

test_that("choice_data object can simulate data", {

  ### simulation
  m <- mod$new(
    formula = choice ~ quality + price | income | time,
    re = "quality",
    alt = 2,
    alt_names = c("bus", "car"),
    ref_alt = "bus"
  )
  set.seed(1)
  d <- choice_data$new()
  theta <- rep(0.1, m$length_theta)
  simulated <- d$simulate(
    mod = m,
    theta = theta,
    N = 2,
    Tp = c(1, 2)
  )
  expect_s3_class(simulated, "data.frame")
  expect_identical(d$df, simulated)
  expect_identical(d$N, 2L)
  expect_identical(d$Tp, c(1L, 2L))
  expect_identical(d$theta_dgp, theta)
  expect_true(all(simulated$choice %in% m$alt_names))
  expect_identical(simulated$occasion, c(1L, 1L, 2L))
  expect_true(all(c("income", "time_bus", "time_car") %in% names(simulated)))
  expect_length(d$build_design(m), 2L)

  ### simulation standard deviation
  zero_sd <- choice_data$new()
  zero_simulated <- zero_sd$simulate(
    mod = m,
    theta = theta,
    N = 2,
    Tp = 1,
    sd = 0
  )
  covariate_columns <- setdiff(
    names(zero_simulated),
    c("id", "occasion", "choice")
  )
  expect_true(all(vapply(
    zero_simulated[covariate_columns],
    function(x) all(x == 0),
    logical(1)
  )))

  ### validation of simulation inputs
  original <- choice_data$new(df = choice_data_wide)
  original_df <- original$df
  original_theta <- original$theta_dgp
  set.seed(2)
  rng_state <- .Random.seed
  expect_error(
    original$simulate(
      mod = m,
      theta = numeric(),
      N = 2,
      Tp = 1
    ),
    "theta"
  )
  expect_error(
    original$simulate(
      mod = m,
      theta = theta,
      N = 2,
      Tp = 1,
      sd = -1
    ),
    "sd"
  )
  expect_identical(.Random.seed, rng_state)
  expect_identical(original$df, original_df)
  expect_identical(original$theta_dgp, original_theta)

  ### resimulate choices for existing covariates
  m_choice <- mod$new(
    formula = choice ~ price | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("bus", "car"),
    ref_alt = "bus",
    H_L = matrix(0, 3, 0),
    f_L = matrix(0, 3, 1)
  )
  choice_covariates <- data.frame(
    id = c(1, 1, 2),
    choice = c("bus", "bus", "bus"),
    price_bus = c(0, 5, 2),
    price_car = c(10, 1, 3)
  )
  choice_sim <- choice_data$new()
  resimulated <- choice_sim$redraw_choices(
    df = choice_covariates,
    mod = m_choice,
    theta = 1,
    column_decider_id = "id"
  )
  expect_identical(resimulated$choice, c("car", "bus", "car"))
  expect_identical(
    resimulated[names(choice_covariates)[-2]],
    choice_covariates[-2]
  )
  expect_identical(choice_sim$df, resimulated)
  expect_identical(choice_sim$N, 2L)
  expect_identical(choice_sim$Tp, c(2L, 1L))
  expect_identical(choice_sim$theta_dgp, 1)

  resimulated$choice <- "bus"
  expect_identical(
    choice_sim$redraw_choices(resimulated, mod = m_choice)$choice,
    c("car", "bus", "car")
  )
  expect_identical(choice_sim$theta_dgp, 1)

  no_theta <- choice_data$new()
  expect_error(
    no_theta$redraw_choices(
      df = choice_covariates,
      mod = m_choice,
      column_decider_id = "id"
    ),
    "Argument `theta` is required when `theta_dgp` is unavailable.",
    fixed = TRUE
  )
  expect_equal(no_theta$df, data.frame())
})

test_that("choice_data object builds design matrices", {

  ### wide data design
  m <- mod$new(
    formula = choice ~ price | income | time,
    re = character(),
    alt = 2,
    alt_names = c("bus", "car"),
    ref_alt = "bus"
  )
  d <- choice_data$new(
    df = choice_data_wide,
    column_decider_id = "id"
  )
  design <- d$build_design(m)
  expect_identical(length(design), 2L)
  expect_identical(design[[1]]$y, c(1L, 2L))
  expect_equal(
    design[[1]]$X[[1]],
    matrix(
      c(1, 4, 0, 10, 0, 1, 5, 0, 0, 8),
      nrow = 2,
      dimnames = list(
        c("bus", "car"),
        c("price", "income_car", "ASC_car", "time_bus", "time_car")
      )
    )
  )

  differenced <- d$build_design(m, difference = TRUE)
  expect_identical(rownames(differenced[[1]]$X[[1]]), "car - bus")
  expect_equal(differenced[[1]]$X[[1]]["car - bus", "price"], 3)
  expect_identical(rownames(differenced[[1]]$X[[2]]), "bus - car")
  expect_equal(differenced[[1]]$X[[2]]["bus - car", "price"], 1)
  expect_identical(rownames(design[[1]]$X[[1]]), c("bus", "car"))
  expect_identical(
    rownames(d$build_design(m)[[1]]$X[[1]]),
    c("bus", "car")
  )

  ### long data design
  long <- choice_data$new(
    df = choice_data_long,
    format = "long",
    column_decider_id = "id",
    column_occasion_id = "occasion",
    column_alternative = "mode"
  )
  expect_equal(long$build_design(m), design)
  expect_equal(long$build_design(m, difference = TRUE), differenced)

  ### formula reuse
  d <- choice_data$new(
    df = choice_data_wide,
    column_decider_id = "id"
  )
  m <- mod$new(
    formula = choice ~ price | 0 | 0,
    alt = 2,
    alt_names = c("bus", "car"),
    ref_alt = "bus"
  )
  design_price <- d$build_design(m)
  expect_identical(colnames(design_price[[1]]$X[[1]]), "price")

  expect_warning(
    m$formula <- choice ~ price | income | 0,
    "Setting `formula` adapted",
    fixed = TRUE
  )
  design_price_income <- d$build_design(m)
  expect_identical(
    colnames(design_price_income[[1]]$X[[1]]),
    c("price", "income_car", "ASC_car")
  )
})

test_that("choice_data object can plot choice frequencies", {

  ### wide data plot
  d <- choice_data$new(
    df = choice_data_wide,
    column_decider_id = "id"
  )
  p <- d$plot()
  expect_s3_class(p, "ggplot")
  expect_identical(as.character(p$data$alternative), c("bus", "car"))
  expect_identical(p$data$value, c(2, 1))
  expect_identical(p$labels$x, "Alternative")
  expect_identical(p$labels$y, "Count")
  y_scale <- p$scales$get_scales("y")
  expect_equal(y_scale$breaks, 0:2)

  ### plot breaks
  large <- choice_data$new(
    df = data.frame(choice = c(rep("bus", 25), rep("car", 2)))
  )
  large_breaks <- large$plot()$scales$get_scales("y")$breaks
  expect_true(all(large_breaks == as.integer(large_breaks)))
  expect_true(length(large_breaks) < 26L)
  expect_true(25L %in% large_breaks)

  ### long data plot
  long <- choice_data$new(
    df = choice_data_long,
    format = "long",
    column_decider_id = "id",
    column_occasion_id = "occasion",
    column_alternative = "mode"
  )
  expect_equal(long$plot()$data$value, c(2, 1))

  ### validation of plot inputs
  expect_error(
    choice_data$new()$plot(),
    "No choices available to plot.",
    fixed = TRUE
  )
})
