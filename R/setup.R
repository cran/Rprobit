#' Setup for Estimation of Multinomial Probit Models
#'
#' @description
#' This function is usually the starting point for working with the `{Rprobit}`
#' package. It helps to construct an [`Rprobit`] object that is used to fit a
#' multinomial probit model to choice data.
#'
#' - The probit model is specified via `formula`, see details.
#'
#' - If no data is provided (`df = NULL`), data is simulated.
#'
#' - If data is provided, `setup()` checks that design matrices can be built
#'   from `df` for the specified model.
#'
#' - Estimation can be performed in place via the method `Rprobit$fit()` or
#'   without modifying the original object via [`fit()`].
#'
#' - Methods can be used to analyze the model fit, see [`Rprobit-methods`].
#'
#' @inheritParams read_formula
#'
#' @param df \[`data.frame` | `NULL`\]\cr
#' Empirical choice data.
#'
#' Can be `NULL`, in which case data is simulated.
#'
#' @param re \[`character()`\]\cr
#' Covariates with random coefficients.
#'
#' @param format \[`character(1)`\]\cr
#' Format of `df`, either `"wide"` or `"long"`.
#'
#' @param column_decider_id \[`character(1)` | `NULL`\]\cr
#' Decider identifier column in `df`.
#'
#' @param column_occasion_id \[`character(1)` | `NULL`\]\cr
#' Choice occasion identifier column in `df`.
#'
#' @param column_choice \[`character(1)` | `NULL`\]\cr
#' Choice column in `df`.
#'
#' If `NULL`, the left-hand side variable of `formula` is used.
#'
#' @param column_alternative \[`character(1)` | `NULL`\]\cr
#' Alternative column. Required for `format = "long"`.
#'
#' @param alt_names \[`character()` | `NULL`\]\cr
#' Alternative labels.
#'
#' Must be supplied if `df = NULL`; otherwise derived from `df`.
#'
#' @param ref_alt \[`character(1)` | `NULL`\]\cr
#' Reference alternative.
#'
#' @param corr_re \[`logical(1)`\]\cr
#' Whether random effects are correlated, see details
#'
#' @param error_struc \[`character(1)`\]\cr
#' Error structure, see details. One of:
#' - `"diag"`, diagonal matrix `Sigma`. The first Cholesky diagonal entry is
#'   fixed to one, and the remaining diagonal entries are free.
#' - `"fixed"`, `Sigma` fixed to identity matrix.
#' - `"ident"`, `Sigma` proportional to identity matrix. The first Cholesky
#'   diagonal entry shares one common free parameter with all remaining
#'   diagonal entries.
#' - `"full"` (default), fully flexible matrix `Sigma`. The first Cholesky
#'   column is fixed to zero for level identification, `L[2,2]` is fixed to one
#'   for scale identification, and all remaining Cholesky entries are free.
#'
#' @param N \[`integer(1)`\]\cr
#' Number of deciders used when simulating data.
#'
#' @param Tp \[`integer()`\]\cr
#' Number of choice occasions per decider used when simulating data.
#'
#' @param sd \[`numeric(1)`\]\cr
#' Standard deviation for the centered normal distribution used when simulating
#' covariates.
#'
#' @return
#' An [`Rprobit`] object.
#'
#' @inheritSection read_formula Specifying the model formula
#' @inheritSection mod Common covariance structures
#' @inheritSection Rprobit Fitting workflows
#'
#' @examples
#' ### empirical choice data
#' if (requireNamespace("mlogit", quietly = TRUE)) {
#'   data("Train", package = "mlogit")
#'   r_train <- setup(
#'     formula = choice ~ price + time | 0 | 0,
#'     df = Train,
#'     column_decider_id = "id",
#'     column_occasion_id = "choiceid"
#'   )
#'   r_train
#' }
#'
#' ### simulated choice data
#' r_sim <- setup(
#'   formula = choice ~ price + time | income | 0,
#'   alt_names = c("train", "bus", "car"),
#'   N = 20,
#'   Tp = 2
#' )
#' r_sim
#'
#' @export

setup <- function(
    formula,
    df = NULL,
    re = character(),
    format = "wide",
    column_decider_id = NULL,
    column_occasion_id = NULL,
    column_choice = NULL,
    column_alternative = NULL,
    alt_names = NULL,
    ref_alt = NULL,
    corr_re = FALSE,
    error_struc = "full",
    N = 100L,
    Tp = 1L,
    sd = 1
  ) {
  ### Check that a model formula was provided.
  if (missing(formula)) {
    stop("Argument `formula` must be specified.", call. = FALSE)
  }
  ### Use the formula left-hand side as choice column if none was supplied.
  if (is.null(column_choice)) {
    column_choice <- read_formula_lhs(formula)
  }
  if (is.null(df)) {
    ### Require alternative labels when data must be simulated.
    if (is.null(alt_names)) {
      stop(
        "Argument `alt_names` is required when `df = NULL`.",
        call. = FALSE
      )
    }
    ### Initialize empty choice data for simulation.
    cd <- choice_data$new()
  } else {
    ### Initialize empirical choice data from the supplied data frame.
    cd <- choice_data$new(
      df = df,
      format = format,
      column_decider_id = column_decider_id,
      column_occasion_id = column_occasion_id,
      column_choice = column_choice,
      column_alternative = column_alternative
    )
  }
  ### Validate supplied alternative labels.
  if (!is.null(alt_names)) {
    checkmate::assert_character(
      alt_names, any.missing = FALSE, min.len = 2, unique = TRUE
    )
    alt_names <- as.character(alt_names)
  }
  ### Derive or check alternative labels against empirical data.
  if (!is.null(df)) {
    data_alt_names <- cd$alt_names
    if (is.null(alt_names)) {
      alt_names <- data_alt_names
    } else {
      compatible_alt_names <- all(alt_names %in% data_alt_names) ||
        all(data_alt_names %in% alt_names)
      if (!compatible_alt_names) {
        stop(
          "Argument `alt_names` is incompatible with `df`.",
          call. = FALSE
        )
      }
    }
  }
  ### Use the first alternative as reference if none was supplied.
  if (is.null(ref_alt)) {
    ref_alt <- alt_names[1]
  }
  ### Create the model object.
  m <- mod$new(
    formula = formula,
    re = re,
    alt = length(alt_names),
    alt_names = alt_names,
    ref_alt = ref_alt,
    corr_re = corr_re,
    error_struc = error_struc
  )
  if (is.null(df)) {
    ### Simulate choice data from the model specification.
    cd$simulate(mod = m, N = N, Tp = Tp, sd = sd)
  } else {
    ### Validate that design matrices can be built from empirical data.
    cd$build_design(mod = m)
  }
  ### Combine the model and choice data into an Rprobit object.
  Rprobit$new(mod = m, choice_data = cd)
}
