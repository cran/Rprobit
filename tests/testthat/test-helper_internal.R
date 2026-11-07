test_that("read_formula parses model formulas correctly", {
  cases <- list(
    list(
      formula = choice ~ A,
      expected = list(
        choice_variable = "choice",
        covariate_types = list("A", character(0), character(0)),
        ASC = TRUE
      )
    ),
    list(
      formula = choice ~ A1 + A2 | B + 0 | C,
      expected = list(
        choice_variable = "choice",
        covariate_types = list(c("A1", "A2"), "B", "C"),
        ASC = FALSE
      )
    ),
    list(
      formula = C ~ A | B,
      expected = list(
        choice_variable = "C",
        covariate_types = list("A", "B", character(0)),
        ASC = TRUE
      )
    ),
    list(
      formula = choice ~ A | 1 | C,
      expected = list(
        choice_variable = "choice",
        covariate_types = list("A", character(0), "C"),
        ASC = TRUE
      )
    ),
    list(
      formula = choice ~ 0 | 0 | 0,
      expected = list(
        choice_variable = "choice",
        covariate_types = list(character(0), character(0), character(0)),
        ASC = FALSE
      )
    )
  )
  for (case in cases) expect_equal(read_formula(case$formula), case$expected)
  expect_error(
    read_formula(choice ~ A | 0 + 1 | C),
    "Terms `0` and `1` cannot both occur in the second formula part.",
    fixed = TRUE
  )
  expect_error(
    read_formula("choice ~ A"),
    "Argument `formula` must be a formula.",
    fixed = TRUE
  )
  expect_error(
    read_formula(~A),
    "Argument `formula` must be a two-sided formula.",
    fixed = TRUE
  )
  expect_error(
    read_formula(choice ~ A | B | C | D),
    "Argument `formula` must have at most three covariate parts.",
    fixed = TRUE
  )
})

test_that("vech_chol_to_chol reconstructs Cholesky root", {
  expect_equal(vech_chol_to_chol(numeric(0)), matrix(0, 0, 0))
  expect_equal(vech_chol_to_chol(2), matrix(2, 1, 1))
  expect_equal(
    vech_chol_to_chol(matrix(1:6, 6, 1)),
    matrix(c(1, 2, 3, 0, 4, 5, 0, 0, 6), 3, 3)
  )
  expect_error(
    vech_chol_to_chol(1:2),
    "Length of `vech_chol` is not triangular.",
    fixed = TRUE
  )
})

test_that("check_covariance_matrix validates covariance matrices", {

  ### valid covariance matrices
  expect_true(check_covariance_matrix(matrix(1, nrow = 1)))
  expect_true(check_covariance_matrix(diag(c(1, 2, 3))))
  expect_true(check_covariance_matrix(matrix(c(1, 0.5, 0.5, 1), nrow = 2)))
  almost_symmetric <- matrix(c(1, 0.5, 0.5 + 1e-10, 1), nrow = 2)
  expect_true(check_covariance_matrix(almost_symmetric, tolerance = 1e-8))
  expect_true(check_covariance_matrix(matrix(c(1, 1, 1, 1), nrow = 2)))
  expect_true(check_covariance_matrix(diag(2), dim = 2))
  expect_true(check_covariance_matrix(matrix(0, 0, 0), dim = 0))

  ### invalid covariance matrices
  invalid_cases <- list(
    list("not numeric", "The value must be numeric."),
    list(c(1, 2, 3), "The value must be a matrix."),
    list(matrix(1:6, nrow = 2, ncol = 3), "The matrix must be square."),
    list(
      matrix(c(1, NA, NA, 1), nrow = 2),
      "The matrix must not contain missing values."
    ),
    list(
      matrix(c(1, Inf, Inf, 1), nrow = 2),
      "The matrix must contain only finite values."
    ),
    list(
      matrix(c(1, 0.5, 0.7, 1), nrow = 2),
      "The matrix must be symmetric."
    ),
    list(
      matrix(c(1, 2, 2, 1), nrow = 2),
      "The matrix must have only positive eigenvalues."
    ),
    list(diag(2), "The matrix must have dimension 3.", list(dim = 3)),
    list(
      matrix(0, 0, 0), "The matrix must have dimension 1.",
      list(dim = 1)
    )
  )
  for (case in invalid_cases) {
    args <- c(list(x = case[[1]]), if (length(case) > 2) case[[3]] else list())
    expect_no_warning(expect_false(do.call(check_covariance_matrix, args)))
    expect_warning(
      expect_false(
        do.call(check_covariance_matrix, c(args, list(warning = TRUE)))
      ),
      case[[2]],
      fixed = TRUE
    )
  }
})

test_that("adjust_covariance_matrix stabilizes covariance estimates", {
  covariance <- matrix(c(2, 0.5, 0.5, 1), nrow = 2)
  adjusted <- adjust_covariance_matrix(covariance)
  expect_equal(adjusted$mat, covariance)
  expect_equal(adjusted$inverse, solve(covariance))

  nonsymmetric <- matrix(c(1, 0.4, 0.2, 1), nrow = 2)
  adjusted <- adjust_covariance_matrix(nonsymmetric)
  expect_equal(adjusted$mat, (nonsymmetric + t(nonsymmetric)) / 2)

  indefinite <- matrix(c(1, 2, 2, 1), nrow = 2)
  expect_warning(
    adjusted <- adjust_covariance_matrix(indefinite, label = "test matrix"),
    "The fitted test matrix was not positive definite and was adjusted.",
    fixed = TRUE
  )
  expect_true(check_covariance_matrix(adjusted$mat))
  expect_equal(adjusted$mat %*% adjusted$inverse, diag(2), tolerance = 1e-8)

  nonfinite <- matrix(c(1, Inf, Inf, 1), nrow = 2)
  expect_warning(
    adjusted <- adjust_covariance_matrix(nonfinite),
    paste0(
      "The fitted covariance matrix contained non-finite values and was ",
      "replaced by an identity matrix."
    ),
    fixed = TRUE
  )
  expect_equal(adjusted$mat, diag(2))
  expect_equal(adjusted$inverse, diag(2))

  empty <- matrix(numeric(), nrow = 0, ncol = 0)
  expect_equal(
    adjust_covariance_matrix(empty),
    list(mat = empty, inverse = empty)
  )
  expect_error(
    adjust_covariance_matrix(matrix(1:6, nrow = 2)),
    "`x` must be a square matrix.",
    fixed = TRUE
  )
  expect_error(
    adjust_covariance_matrix(diag(2), tolerance = 0),
    "`tolerance` must be positive.",
    fixed = TRUE
  )
})

test_that("complete_long_data_frame completes ModeCanada choice sets", {
  skip_if_not_installed("mlogit")
  utils::data("ModeCanada", package = "mlogit", envir = environment())

  alt_names <- unique(as.character(ModeCanada$alt))
  completed <- complete_long_data_frame(
    df = ModeCanada,
    column_choice = "choice",
    column_decider_id = "case",
    column_occasion_id = NULL,
    alt_names = alt_names,
    column_alternative = "alt"
  )

  expected_rows <- length(unique(ModeCanada$case)) * length(alt_names)
  expect_equal(nrow(completed), expected_rows)
  expect_identical(names(completed), c(names(ModeCanada), "occasion"))
  expect_identical(
    anyDuplicated(completed[c("case", "occasion", "alt")]),
    0L
  )
  expect_true(all(table(completed$case) == length(alt_names)))
  expect_true(all(tapply(completed$choice, completed$case, sum) == 1))
})
