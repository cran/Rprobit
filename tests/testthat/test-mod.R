test_that("mod object can be defined", {

  m <- mod$new(
    formula = choice ~ price | 1 | time,
    re = "ASC",
    alt = 3,
    corr_re = FALSE,
    error_struc = "full"
  )
  expect_s3_class(m, "mod")
  expect_identical(m$formula, choice ~ price | 1 | time)
  expect_identical(m$alt_names, c("A", "B", "C"))
  expect_identical(m$alt, 3L)
  expect_identical(m$ref_alt, "A")
  expect_identical(m$re, "ASC")
  expect_identical(m$b_count, c(4L, 2L))
  expect_identical(
    m$var_names,
    list(
      fixed = c("price", "time_A", "time_B", "time_C"),
      random = c("ASC_B", "ASC_C")
    )
  )
  expect_identical(m$length_theta_b, 6L)
  expect_identical(m$length_theta_O, 2L)
  expect_identical(m$length_theta_L, 2L)
  expect_identical(m$length_theta, 10L)
  expect_identical(unname(m$H_b), unname(diag(6)))
  expect_identical(unname(m$f_b), unname(matrix(0, 6, 1)))
  expect_identical(
    unname(m$H_O),
    unname(matrix(c(1, 0, 0, 0, 0, 1), 3, 2))
  )
  expect_identical(unname(m$f_O), unname(matrix(0, 3, 1)))
  expect_identical(unname(m$H_L), unname(diag(6)[, 5:6]))
  expect_identical(
    unname(m$f_L),
    unname(matrix(c(0, 0, 0, 1, 0, 0), 6, 1))
  )
  expect_identical(
    rownames(m$H_b),
    c(
      "b[1]:price", "b[2]:time_A", "b[3]:time_B",
      "b[4]:time_C", "b[5]:ASC_B", "b[6]:ASC_C"
    )
  )
  expect_identical(
    rownames(m$f_b),
    c("price", "time_A", "time_B", "time_C", "ASC_B", "ASC_C")
  )
  expect_identical(colnames(m$H_b), paste0("theta_", 1:6))
  expect_identical(colnames(m$f_b), " ")
  expect_identical(
    rownames(m$H_O),
    c("ASC_B_:_ASC_B", "ASC_C_:_ASC_B", "ASC_C_:_ASC_C")
  )
  expect_identical(
    rownames(m$f_O),
    c("ASC_B_:_ASC_B", "ASC_C_:_ASC_B", "ASC_C_:_ASC_C")
  )
  expect_identical(colnames(m$H_O), paste0("theta_", 7:8))
  expect_identical(colnames(m$f_O), " ")
  expect_identical(
    rownames(m$H_L),
    c("A_:_A", "B_:_A", "C_:_A", "B_:_B", "C_:_B", "C_:_C")
  )
  expect_identical(
    rownames(m$f_L),
    c("A_:_A", "B_:_A", "C_:_A", "B_:_B", "C_:_B", "C_:_C")
  )
  expect_identical(colnames(m$H_L), paste0("theta_", 9:10))
  expect_identical(colnames(m$f_L), " ")

  expect_snapshot(print(m))
  expect_snapshot(m$print())

  ### validation of constructor inputs
  expect_error(
    mod$new(
      formula = choice ~ 0 | 0 | 0,
      re = character(),
      alt = 2,
      alt_names = c("1", "2"),
      ref_alt = "1",
      corr_re = FALSE,
      error_struc = "full",
      matrix(1, 1, 1)
    ),
    "All `...` arguments must be named.",
    fixed = TRUE
  )
  expect_error(
    mod$new(HX = matrix(1, 1, 1)),
    "The following `...` argument(s) are unknown: `HX`.",
    fixed = TRUE
  )
  expect_error(
    mod$new(),
    "`formula` must be specified.",
    fixed = TRUE
  )
  expect_error(
    mod$new(
      formula = choice ~ A | B + 0 | 0,
      alt = 3,
      alt_names = c("1", "2", "3"),
      re = "B",
      H_b = matrix(0, 2, 0),
      f_b = matrix(0, 3, 1),
      H_O = matrix(0, 3, 0),
      H_L = matrix(0, 3, 0)
    ),
    "`H_b` must have 3 row(s)",
    fixed = TRUE
  )
  expect_error(
    mod$new(
      formula = choice ~ A | 0 | 0,
      re = character(),
      alt = 2,
      alt_names = c("1", "2"),
      H_b = matrix(0, 1, 0),
      H_O = matrix(0, 0, 0),
      H_L = matrix(0, 3, 0),
      f_L = matrix(0, 6, 1)
    ),
    "`f_L` must have 3 row(s)",
    fixed = TRUE
  )

})

test_that("mod object accepts custom restrictions", {

  expect_no_warning(
    m <- mod$new(
      formula = choice ~ A | B + 0 | 0,
      alt = 3,
      alt_names = c("1", "2", "3"),
      re = "B",
      error_struc = "fixed",
      H_b = matrix(c(1, 0, 1, 0, 1, 1), 3, 2),
      f_b = matrix(c(1, 2, 3), 3, 1),
      H_O = matrix(1, 3, 1),
      f_O = matrix(0, 3, 1),
      H_L = matrix(1, 6, 2),
      f_L = matrix(0, 6, 1)
    )
  )
  expect_s3_class(m, "mod")
  expect_identical(m$formula, choice ~ A | B + 0 | 0)
  expect_identical(m$alt_names, c("1", "2", "3"))
  expect_identical(m$re, "B")
  expect_identical(m$b_count, c(1L, 2L))
  expect_identical(
    m$var_names,
    list(fixed = "A", random = c("B_2", "B_3"))
  )
  expect_identical(m$length_theta_b, 2L)
  expect_identical(m$length_theta_O, 1L)
  expect_identical(m$length_theta_L, 2L)
  expect_identical(m$length_theta, 5L)
  expect_identical(
    unname(m$H_b),
    unname(matrix(c(1, 0, 1, 0, 1, 1), 3, 2))
  )
  expect_identical(unname(m$f_b), unname(matrix(c(1, 2, 3), 3, 1)))
  expect_identical(unname(m$H_O), unname(matrix(1, 3, 1)))
  expect_identical(unname(m$f_O), unname(matrix(0, 3, 1)))
  expect_identical(unname(m$H_L), unname(matrix(1, 6, 2)))
  expect_identical(unname(m$f_L), unname(matrix(0, 6, 1)))

})

test_that("mod object handles reference alternative", {

  m <- mod$new(
    formula = choice ~ price | income | time,
    re = c("ASC", "time"),
    alt = 3,
    alt_names = c("train", "bus", "car"),
    ref_alt = "bus"
  )
  expect_identical(m$ref_alt, "bus")
  expect_identical(
    m$var_names,
    list(
      fixed = c("price", "income_train", "income_car"),
      random = c(
        "ASC_train", "ASC_car", "time_train", "time_bus", "time_car"
      )
    )
  )
  m$alt_names <- c("walk", "bus", "bike")
  expect_identical(m$ref_alt, "bus")
  expect_identical(
    m$var_names,
    list(
      fixed = c("price", "income_walk", "income_bike"),
      random = c(
        "ASC_walk", "ASC_bike", "time_walk", "time_bus", "time_bike"
      )
    )
  )
  expect_warning(
    m$ref_alt <- "bike",
    "Setting `ref_alt` adapted",
    fixed = TRUE
  )
  expect_identical(m$ref_alt, "bike")
  expect_identical(
    m$var_names,
    list(
      fixed = c("price", "income_walk", "income_bus"),
      random = c(
        "ASC_walk", "ASC_bus", "time_walk", "time_bus", "time_bike"
      )
    )
  )
  expect_error(
    m$ref_alt <- "car",
    "`ref_alt` is not in `alt_names`.",
    fixed = TRUE
  )
  expect_error(
    m$alt_names <- c("walk", "bus", "car"),
    "`ref_alt` is not in `alt_names`.",
    fixed = TRUE
  )

})

test_that("mod object sets covariance restriction defaults", {

  m_correlated <- mod$new(
    formula = choice ~ A1 + A2 + A3 | 0 | 0,
    re = c("A1", "A2", "A3"),
    alt = 2,
    alt_names = c("1", "2"),
    corr_re = TRUE
  )
  expect_s3_class(m_correlated, "mod")
  expect_identical(m_correlated$formula, choice ~ A1 + A2 + A3 | 0 | 0)
  expect_identical(m_correlated$alt_names, c("1", "2"))
  expect_identical(m_correlated$re, c("A1", "A2", "A3"))
  expect_identical(m_correlated$length_theta_b, 3L)
  expect_identical(m_correlated$length_theta_O, 6L)
  expect_identical(m_correlated$length_theta_L, 0L)
  expect_identical(m_correlated$length_theta, 9L)
  expect_identical(unname(m_correlated$H_b), unname(diag(3)))
  expect_identical(unname(m_correlated$f_b), unname(matrix(0, 3, 1)))
  expect_identical(unname(m_correlated$H_O), unname(diag(6)))
  expect_identical(unname(m_correlated$f_O), unname(matrix(0, 6, 1)))
  expect_identical(unname(m_correlated$H_L), unname(matrix(0, 3, 0)))
  expect_identical(unname(m_correlated$f_L), unname(matrix(c(0, 0, 1), 3, 1)))

  m_diag <- mod$new(
    formula = choice ~ A | 0 | 0,
    re = character(),
    alt = 4,
    alt_names = as.character(1:4),
    error_struc = "diag"
  )
  expect_s3_class(m_diag, "mod")
  expect_identical(m_diag$formula, choice ~ A | 0 | 0)
  expect_identical(m_diag$alt_names, as.character(1:4))
  expect_identical(m_diag$re, character())
  expect_identical(m_diag$length_theta_b, 1L)
  expect_identical(m_diag$length_theta_O, 0L)
  expect_identical(m_diag$length_theta_L, 3L)
  expect_identical(m_diag$length_theta, 4L)
  expect_identical(unname(m_diag$H_b), unname(matrix(1, 1, 1)))
  expect_identical(unname(m_diag$f_b), unname(matrix(0, 1, 1)))
  expect_identical(unname(m_diag$H_O), unname(matrix(0, 0, 0)))
  expect_identical(unname(m_diag$f_O), unname(matrix(0, 0, 1)))
  expect_identical(unname(m_diag$H_L), unname(diag(10)[, c(5, 8, 10)]))
  expect_identical(unname(m_diag$H_L[c(5, 8, 10), ]), diag(3))
  expect_identical(unname(m_diag$f_L), unname(matrix(c(1, rep(0, 9)), 10, 1)))

  m_full <- mod$new(
    formula = choice ~ A | 0 | 0,
    re = character(),
    alt = 3,
    alt_names = c("1", "2", "3"),
    error_struc = "full"
  )
  expect_s3_class(m_full, "mod")
  expect_identical(m_full$formula, choice ~ A | 0 | 0)
  expect_identical(m_full$alt_names, c("1", "2", "3"))
  expect_identical(m_full$re, character())
  expect_identical(m_full$length_theta_b, 1L)
  expect_identical(m_full$length_theta_O, 0L)
  expect_identical(m_full$length_theta_L, 2L)
  expect_identical(m_full$length_theta, 3L)
  expect_identical(unname(m_full$H_b), unname(matrix(1, 1, 1)))
  expect_identical(unname(m_full$f_b), unname(matrix(0, 1, 1)))
  expect_identical(unname(m_full$H_O), unname(matrix(0, 0, 0)))
  expect_identical(unname(m_full$f_O), unname(matrix(0, 0, 1)))
  expect_identical(unname(m_full$H_L), unname(diag(6)[, 5:6]))
  expect_identical(
    unname(m_full$f_L),
    unname(matrix(c(0, 0, 0, 1, 0, 0), 6, 1))
  )

  m_ident <- mod$new(
    formula = choice ~ A | 0 | 0,
    re = character(),
    alt = 3,
    alt_names = c("1", "2", "3"),
    error_struc = "ident"
  )
  expect_s3_class(m_ident, "mod")
  expect_identical(m_ident$formula, choice ~ A | 0 | 0)
  expect_identical(m_ident$alt_names, c("1", "2", "3"))
  expect_identical(m_ident$re, character())
  expect_identical(m_ident$length_theta_b, 1L)
  expect_identical(m_ident$length_theta_O, 0L)
  expect_identical(m_ident$length_theta_L, 1L)
  expect_identical(m_ident$length_theta, 2L)
  expect_identical(unname(m_ident$H_b), unname(matrix(1, 1, 1)))
  expect_identical(unname(m_ident$f_b), unname(matrix(0, 1, 1)))
  expect_identical(unname(m_ident$H_O), unname(matrix(0, 0, 0)))
  expect_identical(unname(m_ident$f_O), unname(matrix(0, 0, 1)))
  expect_identical(
    unname(m_ident$H_L),
    unname(matrix(c(1, 0, 0, 1, 0, 1), 6, 1))
  )
  expect_identical(
    unname(m_ident$H_L[c(1, 4, 6), , drop = FALSE]),
    matrix(1, 3, 1)
  )
  expect_identical(unname(m_ident$f_L), unname(matrix(0, 6, 1)))

  m_fixed <- mod$new(
    formula = choice ~ A | 0 | 0,
    re = character(),
    alt = 3,
    alt_names = c("1", "2", "3"),
    error_struc = "fixed"
  )
  expect_s3_class(m_fixed, "mod")
  expect_identical(m_fixed$formula, choice ~ A | 0 | 0)
  expect_identical(m_fixed$alt_names, c("1", "2", "3"))
  expect_identical(m_fixed$re, character())
  expect_identical(m_fixed$length_theta_b, 1L)
  expect_identical(m_fixed$length_theta_O, 0L)
  expect_identical(m_fixed$length_theta_L, 0L)
  expect_identical(m_fixed$length_theta, 1L)
  expect_identical(unname(m_fixed$H_b), unname(matrix(1, 1, 1)))
  expect_identical(unname(m_fixed$f_b), unname(matrix(0, 1, 1)))
  expect_identical(unname(m_fixed$H_O), unname(matrix(0, 0, 0)))
  expect_identical(unname(m_fixed$f_O), unname(matrix(0, 0, 1)))
  expect_identical(unname(m_fixed$H_L), unname(matrix(0, 6, 0)))
  expect_identical(
    unname(m_fixed$f_L),
    matrix(c(1, 0, 0, 1, 0, 1), 6, 1)
  )

  m_toggle <- mod$new(
    formula = choice ~ price | 1 | time,
    re = "ASC",
    alt = 3,
    alt_names = c("train", "bus", "car"),
    ref_alt = "train",
    corr_re = FALSE,
    error_struc = "diag"
  )
  expect_warning(
    m_toggle$re <- character(),
    "Setting `re` adapted `H_O`, `f_O`"
  )
  expect_identical(unname(m_toggle$H_O), unname(matrix(0, 0, 0)))
  expect_identical(unname(m_toggle$f_O), unname(matrix(0, 0, 1)))
  expect_identical(m_toggle$length_theta_O, 0L)

  expect_warning(
    m_toggle$re <- "ASC",
    "Setting `re` adapted `H_O`, `f_O`"
  )
  expect_identical(unname(m_toggle$H_O), unname(diag(3)[, c(1, 3)]))
  expect_identical(unname(m_toggle$f_O), unname(matrix(0, 3, 1)))
  expect_identical(m_toggle$length_theta_O, 2L)

  expect_warning(
    m_toggle$alt_names <- c("train", "bus"),
    "Setting `alt_names` adapted"
  )
  expect_identical(unname(m_toggle$H_L), unname(matrix(0, 3, 0)))
  expect_identical(
    unname(m_toggle$f_L),
    unname(matrix(c(0, 0, 1), 3, 1))
  )
  expect_identical(m_toggle$length_theta_L, 0L)

  expect_warning(
    m_toggle$alt_names <- c("train", "bus", "car"),
    "Setting `alt_names` adapted"
  )
  expect_identical(unname(m_toggle$H_L), unname(diag(6)[, 5:6]))
  expect_identical(
    unname(m_toggle$f_L),
    unname(matrix(c(0, 0, 0, 1, 0, 0), 6, 1))
  )
  expect_identical(m_toggle$length_theta_L, 2L)

})

test_that("mod object keeps model components conformable", {

  mod_test <- mod$new(
    formula = choice ~ 0 | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("1", "2")
  )

  ### change number of alternatives
  expect_warning(
    mod_test$alt_names <- c("1", "2", "3"),
    "Setting `alt_names` adapted `alt`, `H_L`, `f_L`"
  )
  expect_identical(mod_test$alt_names, c("1", "2", "3"))
  expect_identical(mod_test$alt, 3L)
  expect_identical(unname(mod_test$H_L), unname(diag(6)[, 5:6]))
  expect_identical(
    unname(mod_test$f_L),
    unname(matrix(c(0, 0, 0, 1, 0, 0), 6, 1))
  )
  expect_identical(mod_test$length_theta_L, 2L)
  expect_identical(mod_test$length_theta, 2L)

  ### change formula and random coefficients
  expect_warning(
    mod_test$formula <- choice ~ A1 + A2 | B + 0 | 0,
    "Setting `formula` adapted `H_b`, `f_b`"
  )
  expect_warning(
    mod_test$re <- c("A1", "A2"),
    "Setting `re` adapted `H_O`, `f_O`"
  )
  expect_identical(unname(mod_test$H_b), unname(matrix(0, 4, 0)))
  expect_identical(unname(mod_test$f_b), unname(matrix(0, 4, 1)))
  expect_identical(unname(mod_test$H_O), unname(diag(3)[, c(1, 3)]))
  expect_identical(unname(mod_test$f_O), unname(matrix(0, 3, 1)))
  expect_identical(mod_test$length_theta_b, 0L)
  expect_identical(mod_test$length_theta_O, 2L)
  expect_identical(mod_test$length_theta, 4L)
  expect_identical(mod_test$b_count, c(2L, 2L))

  ### change b restrictions
  mod_test$H_b <- matrix(1, nrow = 4, ncol = 2)
  expect_identical(unname(mod_test$H_b), unname(matrix(1, 4, 2)))
  expect_identical(mod_test$length_theta_b, 2L)
  expect_identical(mod_test$length_theta, 6L)
  expect_identical(unname(mod_test$f_b), unname(matrix(0, 4, 1)))

  mod_test$f_b <- matrix(1, nrow = 4, ncol = 1)
  expect_identical(unname(mod_test$f_b), unname(matrix(1, 4, 1)))
  expect_identical(mod_test$length_theta_b, 2L)
  expect_identical(mod_test$length_theta, 6L)
  expect_identical(unname(mod_test$H_b), unname(matrix(1, 4, 2)))

  mod_test$f_b <- 1:4
  expect_identical(unname(mod_test$f_b), unname(matrix(as.numeric(1:4), 4, 1)))
  expect_identical(mod_test$length_theta_b, 2L)
  expect_identical(mod_test$length_theta, 6L)

  ### change Omega restrictions
  mod_test$H_O <- matrix(1, nrow = 3, ncol = 1)
  expect_identical(unname(mod_test$H_O), unname(matrix(1, 3, 1)))
  expect_identical(mod_test$length_theta_O, 1L)
  expect_identical(mod_test$length_theta, 5L)
  expect_identical(unname(mod_test$f_O), unname(matrix(0, 3, 1)))

  expect_error(
    mod_test$f_O <- matrix(1, nrow = 6, ncol = 1),
    "`f_O` must have 3 row(s)",
    fixed = TRUE
  )
  mod_test$f_O <- 1:3
  expect_identical(unname(mod_test$f_O), unname(matrix(as.numeric(1:3), 3, 1)))
  expect_identical(mod_test$length_theta_O, 1L)
  expect_identical(mod_test$length_theta, 5L)
  expect_identical(unname(mod_test$H_O), unname(matrix(1, 3, 1)))

  ### change Sigma restrictions
  mod_test$H_L <- matrix(1, nrow = 6, ncol = 2)
  expect_identical(unname(mod_test$H_L), unname(matrix(1, 6, 2)))
  expect_identical(mod_test$length_theta_L, 2L)
  expect_identical(mod_test$length_theta, 5L)
  expect_identical(
    unname(mod_test$f_L),
    unname(matrix(c(0, 0, 0, 1, 0, 0), 6, 1))
  )

  expect_warning(
    mod_test$alt_names <- c("1", "2", "3", "4"),
    "Setting `alt_names` adapted `H_b`, `f_b`, `alt`, `H_L`, `f_L`"
  )
  expect_identical(mod_test$alt_names, c("1", "2", "3", "4"))
  expect_identical(mod_test$alt, 4L)
  expect_identical(mod_test$b_count, c(3L, 2L))
  expect_identical(mod_test$length_theta_L, 5L)
  expect_identical(mod_test$length_theta, 8L)
  expect_identical(
    unname(mod_test$H_L),
    unname(diag(10)[, 6:10])
  )
  expect_identical(
    unname(mod_test$f_L),
    unname(matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0), 10, 1))
  )

  mod_quiet <- mod$new(
    formula = choice ~ 0 | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("1", "2")
  )
  mod_quiet$suppress_warnings <- TRUE
  expect_no_warning(mod_quiet$alt_names <- c("1", "2", "3"))
  expect_identical(mod_quiet$alt, 3L)

  mod_test$f_L <- 1:10
  expect_identical(
    unname(mod_test$f_L),
    unname(matrix(as.numeric(1:10), 10, 1))
  )
  expect_identical(mod_test$length_theta_L, 5L)
  expect_identical(mod_test$length_theta, 8L)

  ### check edge cases
  expect_error(mod_test$length_theta_b <- 1, "read-only", fixed = TRUE)
  expect_error(mod_test$length_theta_O <- 1, "read-only", fixed = TRUE)
  expect_error(mod_test$length_theta_L <- 1, "read-only", fixed = TRUE)
  expect_error(mod_test$length_theta <- 1, "read-only", fixed = TRUE)
  expect_error(mod_test$alt <- 2, "read-only", fixed = TRUE)
  expect_error(mod_test$var_names <- list(), "read-only", fixed = TRUE)
  expect_warning(
    mod_test$ref_alt <- "3",
    "Setting `ref_alt` adapted",
    fixed = TRUE
  )
  expect_identical(mod_test$ref_alt, "3")
  expect_error(mod_test$ref_alt <- 2, "string", fixed = TRUE)
  expect_error(
    mod_test$ref_alt <- "5",
    "`ref_alt` is not in `alt_names`.",
    fixed = TRUE
  )
  expect_error(mod_test$b_count <- c(1, 1), "read-only", fixed = TRUE)
  expect_error(
    mod_test$H_L <- matrix(0, nrow = 2, ncol = 1),
    "`H_L` must have 10 row(s)",
    fixed = TRUE
  )
  expect_error(
    mod_test$f_b <- matrix(1, nrow = sum(mod_test$b_count), ncol = 2),
    "`f_b` must have one column.",
    fixed = TRUE
  )
  expect_error(
    mod_test$f_O <- matrix(1, nrow = 3, ncol = 2),
    "`f_O` must have one column.",
    fixed = TRUE
  )
  expect_error(
    mod_test$f_L <- matrix(1, nrow = 10, ncol = 2),
    "`f_L` must have one column.",
    fixed = TRUE
  )
  expect_error(
    mod_test$f_b <- "not numeric",
    "`f_b` must be a matrix or numeric vector.",
    fixed = TRUE
  )
  expect_error(
    mod_test$f_O <- "not numeric",
    "`f_O` must be a matrix or numeric vector.",
    fixed = TRUE
  )
  expect_error(
    mod_test$f_L <- "not numeric",
    "`f_L` must be a matrix or numeric vector.",
    fixed = TRUE
  )
  expect_error(
    mod_test$H_b <- matrix(
      "not numeric",
      nrow = sum(mod_test$b_count),
      ncol = 1
    ),
    "numeric"
  )
  expect_error(
    mod_test$H_O <- matrix(NA_real_, nrow = 3, ncol = 1),
    "missing"
  )
  expect_error(
    mod_test$H_L <- matrix(Inf, nrow = 10, ncol = 1),
    "`H_L` must contain only finite values.",
    fixed = TRUE
  )
  expect_error(
    mod$new(
      formula = choice ~ A | 0 | 0,
      re = "Z",
      alt = 2,
      alt_names = c("1", "2")
    ),
    "The following random coefficient variable(s) are unknown",
    fixed = TRUE
  )
  expect_error(
    mod$new(
      formula = choice ~ A | 0 | 0,
      re = character(),
      alt = 2,
      alt_names = c("1", "2"),
      error_struc = "bad"
    )
  )
  expect_error(
    mod$new(
      formula = choice ~ A | 0 | 0,
      re = character(),
      alt = 2,
      alt_names = c("1", "2"),
      corr_re = NA
    )
  )
  expect_error(
    mod$new(
      formula = choice ~ A | 0 | 0,
      re = character(),
      alt = 2,
      alt_names = "1"
    ),
    "Argument `alt_names` must have length `alt`.",
    fixed = TRUE
  )
  expect_error(
    mod$new(
      formula = choice ~ A | 0 | 0,
      re = character(),
      alt = 2,
      alt_names = c("1", "2"),
      ref_alt = 3
    ),
    "string",
    fixed = TRUE
  )
  expect_error(
    mod$new(
      formula = choice ~ A | 0 | 0,
      re = character(),
      alt = 2,
      alt_names = c("1", "2"),
      ref_alt = "3"
    ),
    "`ref_alt` is not in `alt_names`.",
    fixed = TRUE
  )

})

test_that("mod object builds model parameters from theta", {

  ### default mod object
  m_default <- mod$new(
    formula = choice ~ 0 | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("1", "2")
  )
  parameter <- m_default$build_parameter(numeric(0))
  expect_equal(unname(parameter$b), unname(m_default$f_b))
  expect_equal(unname(parameter$b), unname(matrix(0, 0, 1)))

  ### custom mod object
  m_custom <- mod$new(
    formula = choice ~ A1 + A2 | B + 0 | 0,
    re = c("A1", "A2"),
    alt = 3,
    alt_names = c("1", "2", "3")
  )
  theta <- seq_len(m_custom$length_theta)
  parameter <- m_custom$build_parameter(theta)
  expect_named(
    parameter,
    c("b", "O", "Omega_chol", "Omega", "L", "Sigma_chol", "Sigma")
  )
  expect_equal(unname(parameter$b), unname(matrix(1:4, 4, 1)))
  expect_identical(rownames(parameter$b), c("B_2", "B_3", "A1", "A2"))
  expect_identical(colnames(parameter$b), " ")
  expect_equal(unname(parameter$O), unname(matrix(c(5, 0, 6), 3, 1)))
  expect_identical(
    rownames(parameter$O),
    c("A1_:_A1", "A2_:_A1", "A2_:_A2")
  )
  expect_identical(colnames(parameter$O), " ")
  expect_equal(
    unname(parameter$Omega_chol),
    matrix(c(5, 0, 0, 6), 2, 2)
  )
  expect_identical(rownames(parameter$Omega_chol), c("A1", "A2"))
  expect_identical(colnames(parameter$Omega_chol), c("A1", "A2"))
  expect_equal(unname(parameter$Omega), unname(diag(c(25, 36))))
  expect_identical(rownames(parameter$Omega), c("A1", "A2"))
  expect_identical(colnames(parameter$Omega), c("A1", "A2"))
  expect_equal(unname(parameter$L), unname(matrix(c(0, 0, 0, 1, 7, 8), 6, 1)))
  expect_identical(
    rownames(parameter$L),
    c("1_:_1", "2_:_1", "3_:_1", "2_:_2", "3_:_2", "3_:_3")
  )
  expect_identical(colnames(parameter$L), " ")
  expect_equal(
    unname(parameter$Sigma_chol),
    matrix(c(0, 0, 0, 0, 1, 7, 0, 0, 8), 3, 3)
  )
  expect_identical(rownames(parameter$Sigma_chol), c("1", "2", "3"))
  expect_identical(colnames(parameter$Sigma_chol), c("1", "2", "3"))
  expect_equal(
    unname(parameter$Sigma),
    matrix(c(0, 0, 0, 0, 1, 7, 0, 7, 113), 3, 3)
  )
  expect_identical(rownames(parameter$Sigma), c("1", "2", "3"))
  expect_identical(colnames(parameter$Sigma), c("1", "2", "3"))

  ### mod object with correlated random effects
  m_correlated <- mod$new(
    formula = choice ~ A1 + A2 + A3 | 0 | 0,
    re = c("A1", "A2", "A3"),
    alt = 2,
    alt_names = c("1", "2"),
    corr_re = TRUE
  )
  theta <- seq_len(m_correlated$length_theta)
  parameter <- m_correlated$build_parameter(theta)
  expect_equal(unname(parameter$b), unname(matrix(1:3, 3, 1)))
  expect_equal(unname(parameter$O), unname(matrix(4:9, 6, 1)))
  expect_equal(
    unname(parameter$Omega_chol),
    matrix(c(4, 5, 6, 0, 7, 8, 0, 0, 9), 3, 3)
  )
  expect_equal(
    unname(parameter$Omega),
    unname(parameter$Omega_chol %*% t(parameter$Omega_chol))
  )
  expect_equal(unname(parameter$L), unname(matrix(c(0, 0, 1), 3, 1)))
  expect_equal(
    unname(parameter$Sigma_chol),
    unname(matrix(c(0, 0, 0, 1), 2, 2))
  )
  expect_equal(unname(parameter$Sigma), unname(matrix(c(0, 0, 0, 1), 2, 2)))

  ### mod object with fixed error structure
  m_fixed <- mod$new(
    formula = choice ~ A | 0 | 0,
    re = character(),
    alt = 3,
    alt_names = c("1", "2", "3"),
    error_struc = "fixed"
  )
  parameter <- m_fixed$build_parameter(1)
  expect_equal(unname(parameter$O), unname(matrix(0, 0, 1)))
  expect_equal(unname(parameter$Omega_chol), unname(matrix(0, 0, 0)))
  expect_equal(unname(parameter$Omega), unname(matrix(0, 0, 0)))
  expect_equal(
    unname(parameter$L),
    matrix(c(1, 0, 0, 1, 0, 1), 6, 1)
  )
  expect_equal(
    unname(parameter$Sigma_chol),
    diag(3)
  )
  expect_equal(unname(parameter$Sigma), unname(diag(3)))
  expect_error(m_fixed$build_parameter(numeric(0)), "length")
})

test_that("mod object builds theta from model parameters", {

  ### theta can be uniquely reproduced
  H_b <- matrix(
    c(1, 0, 2, 0, 1, -1, 1, 1, 0, -1, 0, 1, 0, 2, 1),
    nrow = 5, byrow = TRUE
  )
  f_b <- matrix(c(0.5, -0.25, 1, 0.75, -1.5), ncol = 1)
  H_O <- matrix(
    c(1, 0, 0.5, 1, 0, 1),
    nrow = 3, byrow = TRUE
  )
  f_O <- matrix(c(1, 0.1, 0.5), ncol = 1)
  H_L <- matrix(
    c(1, 0, 0, 0.2, 1, 0, 0, 0.5, 1, 0, 1, 0, 1, -0.5, 0.25, 0, 0, 1),
    nrow = 6, byrow = TRUE
  )
  f_L <- matrix(c(1, -0.1, 0.2, 0.8, -0.3, 1.1), ncol = 1)
  m <- mod$new(
    formula = choice ~ A1 + A2 + A3 | B + 0 | 0,
    alt = 3,
    alt_names = c("1", "2", "3"),
    re = c("A1", "A2"),
    H_b = H_b,
    f_b = f_b,
    H_O = H_O,
    f_O = f_O,
    H_L = H_L,
    f_L = f_L
  )
  theta <- c(0.4, -0.8, 1.2, 0.3, 0.7, 0.5, -0.4, 0.6)
  parameter <- m$build_parameter(theta)
  expect_no_warning(
    theta2 <- m$build_theta(
      b = parameter$b,
      Omega = parameter$Omega,
      Sigma = parameter$Sigma
    )
  )
  expect_equal(theta2, theta)
  expect_error(
    m$build_theta(
      b = parameter$b,
      Omega = diag(3),
      Sigma = parameter$Sigma
    ),
    paste(
      "`Omega` has dimension 3 x 3, but `mod` expects a 2 x 2",
      "covariance matrix to match the random coefficient count."
    ),
    fixed = TRUE
  )
  expect_error(
    m$build_theta(
      b = parameter$b,
      Omega = matrix(1, 2, 2),
      Sigma = parameter$Sigma
    ),
    "`Omega` must be positive definite.",
    fixed = TRUE
  )

  m_without_re <- mod$new(
    formula = choice ~ A | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("1", "2")
  )
  expect_error(
    m_without_re$build_theta(
      b = 0,
      Sigma = matrix(1, 2, 2)
    ),
    "`Sigma` must be positive definite.",
    fixed = TRUE
  )

  ### theta needs to be approximated
  m <- mod$new(
    formula = choice ~ A1 + A2 | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("1", "2"),
    H_b = matrix(1, nrow = 2, ncol = 1)
  )
  expect_warning(
    theta <- m$build_theta(
      b = c(1, 3),
      Sigma = diag(2)
    ),
    "least-squares approximation",
    fixed = TRUE
  )
  expect_equal(theta, 2)
  m$suppress_warnings <- TRUE
  expect_no_warning(
    theta <- m$build_theta(
      b = c(1, 3),
      Sigma = diag(2)
    )
  )
  expect_equal(theta, 2)

})

test_that("mod object can make theta unique", {
  m <- mod$new(
    formula = choice ~ A + B | 0 | 0,
    re = c("A", "B"),
    alt = 3,
    corr_re = TRUE
  )
  theta <- c(0.5, -0.5, -2, 0.5, -3, 0.4, -0.7)
  theta_unique <- c(0.5, -0.5, 2, -0.5, 3, 0.4, 0.7)
  expect_equal(m$make_theta_unique(theta), theta_unique)
  expect_equal(m$make_theta_unique(theta_unique), theta_unique)
})

test_that("mod object checks model identifiability", {

  ### this model is identified
  m_default <- mod$new(
    formula = choice ~ 0 | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("1", "2"),
    corr_re = FALSE,
    error_struc = "full"
  )
  expect_true(m_default$check_identifiability())

  ### the diagonal default includes a scale restriction
  m_diag <- mod$new(
    formula = choice ~ 0 | 0 | 0,
    re = character(),
    alt = 3,
    alt_names = c("1", "2", "3"),
    error_struc = "diag"
  )
  expect_true(m_diag$check_identifiability())

  ### this model is not identified (beta_3 = beta_1 + beta_2, all estimated)
  m_linear <- mod$new(
    formula = choice ~ A1 + A2 + A3 | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("1", "2"),
    H_b = matrix(c(1, 0, 1, 0, 1, 1, 0, 0, 0), nrow = 3, ncol = 3),
    f_b = matrix(c(0, 0, 0), ncol = 1),
    H_O = matrix(0, 0, 0),
    f_O = matrix(0, 0, 0),
    H_L = matrix(0, nrow = 3, ncol = 0),
    f_L = matrix(c(0, 0, 1), ncol = 1)
  )
  expect_warning(
    expect_false(m_linear$check_identifiability()),
    "The parameter set may not be identified; check the specification.",
    fixed = TRUE
  )
  m_linear$suppress_warnings <- TRUE
  expect_no_warning(expect_false(m_linear$check_identifiability()))

  ### not estimating beta_3 makes the model identified
  m_linear2 <- m_linear$clone()
  m_linear2$H_b <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 3, ncol = 2)
  expect_true(m_linear2$check_identifiability())

  ### this model is not identified (no scale restriction)
  m_scale <- mod$new(
    formula = choice ~ A1 + A2 | 0 | 0,
    re = character(),
    alt = 2,
    alt_names = c("1", "2"),
    H_b = diag(2),
    f_b = matrix(c(0, 0), ncol = 1),
    H_O = matrix(0, 0, 0),
    f_O = matrix(0, 0, 0),
    H_L = diag(3),
    f_L = matrix(c(0, 0, 0), ncol = 1)
  )
  expect_warning(
    expect_false(m_scale$check_identifiability()),
    "The parameter set may not be identified; check the specification.",
    fixed = TRUE
  )

  ### fixing one beta entry still leaves too many parameters in Sigma
  m_scale2 <- m_scale$clone()
  m_scale2$H_b <- diag(2)[, -1, drop = FALSE]
  m_scale2$f_b <- matrix(c(1, 0), ncol = 1)
  expect_warning(
    expect_false(m_scale2$check_identifiability()),
    "The parameter set may not be identified; check the specification.",
    fixed = TRUE
  )

  ### fixing one Sigma entry makes the model identified
  m_scale3 <- m_scale$clone()
  m_scale3$H_L <- matrix(0, nrow = 3, ncol = 0)
  m_scale3$f_L <- matrix(c(0, 0, 1), ncol = 1)
  expect_true(m_scale3$check_identifiability())

})
