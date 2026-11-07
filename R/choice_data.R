#' Choice Data
#'
#' @description
#' An object of the `choice_data` R6 class stores the data for choice modeling
#' and can interpret `data.frame`s in wide or long format.
#'
#' \strong{Wide data} contain one row per choice occasion. Alternative-varying
#' covariates are supplied in columns with names of the form
#' `covariate_alternative`, for example `price_bus` and `price_car`;
#' alternative-constant covariates are supplied as ordinary columns. The
#' argument `column_choice` names the column that contains the chosen
#' alternative label. The argument `column_decider_id` optionally names the
#' decider identifier column. If it is `NULL`, each row is treated as one
#' decider with one choice occasion. The argument `column_occasion_id`
#' optionally names a choice occasion identifier column. The argument
#' `column_alternative` is not used for wide data.
#'
#' \strong{Long data} contain one row per alternative and choice occasion.
#' Alternative-varying covariates are supplied as ordinary columns, one value
#' per alternative row; alternative-constant covariates are supplied as columns
#' that are constant within each choice set. The argument `column_alternative`
#' is required and names the column that contains the alternative labels. The
#' argument `column_choice` names the column that identifies the chosen row in
#' each choice set, encoded by a logical indicator or a `0`/`1` numeric
#' indicator. Choice sets are identified by `column_decider_id` and
#' `column_occasion_id` when both are supplied, by consecutive complete
#' alternative blocks within each decider when only `column_decider_id` is
#' supplied, and by consecutive complete alternative blocks when neither is
#' supplied.
#'
#' For a given [`mod`] object,
#'
#' - choice data can be simulated via the method `$simulate()`,
#' - choices can be redrawn via the method `$redraw_choices()`,
#' - design matrices are created via the method `$build_design()`.
#'
#' @param df \[`data.frame`\]\cr
#' Choice data.
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
#' If both `column_decider_id` and `column_occasion_id` are `NULL` and
#' `format = "wide"`, every row is treated as one decider with one choice
#' occasion.
#'
#' If `column_occasion_id` is `NULL` and `format = "long"`, rows are
#' interpreted as consecutive complete choice sets within each decider.
#'
#' @param column_choice \[`character(1)`\]\cr
#' Choice column in `df`.
#'
#' @param column_alternative \[`character(1)` | `NULL`\]\cr
#' Alternative column. Required for `format = "long"`.
#'
#' @param mod \[`mod`\]\cr
#' Model object.
#'
#' @param theta \[`numeric()`\]\cr
#' Free parameter vector.
#'
#' @param N \[`integer(1)`\]\cr
#' Number of deciders.
#'
#' @param Tp \[`integer()`\]\cr
#' Number of choice occasions per decider. Either scalar or length `N`.
#'
#' @param sd \[`numeric(1)`\]\cr
#' Standard deviation of the centered normal distribution used to draw
#' covariates.
#'
#' @param difference \[`logical(1)`\]\cr
#' If `TRUE`, subtract the row of the chosen alternative from each design
#' matrix and drop that row.
#'
#' @param check_choice \[`logical(1)`\]\cr
#' If `TRUE`, warn about missing or invalid choices.
#'
#' @examples
#' ### read choice data
#' if (requireNamespace("mlogit", quietly = TRUE)) {
#'
#'   ### wide data
#'   data("Train", package = "mlogit")
#'   (train_data <- choice_data$new(
#'     df = Train,
#'     column_decider_id = "id",
#'     column_occasion_id = "choiceid"
#'   ))
#'   train_data$choice_frequencies
#'
#'   train_mod <- mod$new(
#'     formula = choice ~ price + time | 0 | 0,
#'     alt = length(train_data$alt_names),
#'     alt_names = train_data$alt_names
#'   )
#'   train_design <- train_data$build_design(train_mod)
#'   train_design[[1]]
#'
#'   ### long data
#'   data("ModeCanada", package = "mlogit")
#'   suppressWarnings(mode_data <- choice_data$new(
#'     df = subset(ModeCanada, noalt == 3),
#'     format = "long",
#'     column_decider_id = "case",
#'     column_alternative = "alt"
#'   ))
#'   mode_data$choice_counts
#'
#'   mode_mod <- mod$new(
#'     formula = choice ~ cost + ivt | income | 0,
#'     alt = length(mode_data$alt_names),
#'     alt_names = mode_data$alt_names
#'   )
#'   mode_design <- mode_data$build_design(mode_mod, difference = TRUE)
#'   mode_design[[1]] # no information about 'air' available
#' }
#'
#' ### simulate choice data
#' sim_mod <- mod$new(
#'   formula = choice ~ price + time | income | 0,
#'   re = "time",
#'   alt = 3,
#'   alt_names = c("train", "bus", "car")
#' )
#' sim_data <- choice_data$new()
#' df <- sim_data$simulate(mod = sim_mod, N = 100, sd = 2)
#' head(df)
#' plot(sim_data)
#'
#' @export

choice_data <- R6::R6Class(
  classname = "choice_data",
  public = list(

    #' @description
    #' Create a new `choice_data` object.
    #'
    #' @return
    #' A new `choice_data` object.
    initialize = function(
      df = data.frame(),
      format = "wide",
      column_decider_id = NULL,
      column_occasion_id = NULL,
      column_choice = "choice",
      column_alternative = NULL
    ) {
      checkmate::assert_data_frame(df)
      df <- as.data.frame(df)
      checkmate::assert_choice(format, choices = c("wide", "long"))
      checkmate::assert_string(
        column_decider_id,
        min.chars = 1,
        null.ok = TRUE
      )
      checkmate::assert_string(
        column_occasion_id,
        min.chars = 1,
        null.ok = TRUE
      )
      checkmate::assert_string(column_choice, min.chars = 1)
      checkmate::assert_string(
        column_alternative,
        min.chars = 1,
        null.ok = TRUE
      )
      if (is.null(column_decider_id) && !is.null(column_occasion_id)) {
        stop(
          "Argument `column_occasion_id` requires `column_decider_id`.",
          call. = FALSE
        )
      }
      if (format == "wide" && !is.null(column_alternative)) {
        warning(
          "Argument `column_alternative` is only used with `format = \"long\"`.",
          call. = FALSE
        )
        column_alternative <- NULL
      }
      if (format == "long" && is.null(column_alternative)) {
        stop(
          "Argument `column_alternative` must be specified for long data.",
          call. = FALSE
        )
      }
      if (nrow(df) > 0L) {
        # add column, if not there
        missing_choice <- setdiff(
          c(column_choice),
          names(df)
        )

        if (length(missing_choice) >0L) {
          df$choice <- NA
          column_choice = "choice"
        }

        missing_columns <- setdiff(
          c(
            if (!is.null(column_decider_id)) column_decider_id,
            if (!is.null(column_occasion_id)) column_occasion_id,
            column_choice,
            if (format == "long") column_alternative
          ),
          names(df)
        )
        if (length(missing_columns) > 0L) {
          stop(
            sprintf(
              "Required column(s) are missing: %s.",
              paste(sprintf("`%s`", missing_columns), collapse = ", ")
            ),
            call. = FALSE
          )
        }
        if (format == "wide" && is.null(column_occasion_id)) {
          occasion_column <- "occasion"
          if (!occasion_column %in% names(df)) {
            if (is.null(column_decider_id)) {
              df[[occasion_column]] <- 1L
            } else {
              ids <- df[[column_decider_id]]
              df[[occasion_column]] <- ave(seq_along(ids), ids, FUN = seq_along)
            }
          }
          if (!is.null(column_decider_id)) {
            column_occasion_id <- occasion_column
          }
        }
        if (!is.null(column_decider_id) &&
            any(is.na(df[[column_decider_id]]))) {
          stop(
            "Decider id column must not contain missing values.",
            call. = FALSE
          )
        }
        if (!is.null(column_occasion_id) &&
            any(is.na(df[[column_occasion_id]]))) {
          stop(
            "Occasion id column must not contain missing values.",
            call. = FALSE
          )
        }
        if (format == "long") {
          alt_names <- unique(df[[column_alternative]])
          # complete data set
          df <- complete_long_data_frame(
            df,
            column_choice,
            column_decider_id,
            column_occasion_id,
            alt_names,
            column_alternative
          )
          if (
            is.null(column_occasion_id) && !is.null(column_decider_id)
            && "occasion" %in% names(df)
          ) {
            column_occasion_id <- "occasion"
          }
        }

        if (
          format == "wide" &&
            !is.null(column_decider_id) &&
            !is.null(column_occasion_id)
        ) {
          key <- df[c(column_decider_id, column_occasion_id)]
          if (anyDuplicated(key) > 0L) {
            stop(
              paste(
                "Wide data must have unique decider and occasion",
                "identifier rows."
              ),
              call. = FALSE
            )
          }
        }
        if (any(is.na(df[[column_choice]]))) {
          warning(
            "Choice column contains missing values. Correct the data frame, if this is not wanted.",
            call. = FALSE
          )
        }
        if (format == "long" && any(is.na(df[[column_alternative]]))) {
          stop(
            "Alternative column must not contain missing values.",
            call. = FALSE
          )
        }
        if (format == "long") {
          choice_sets <- private$.group_long_choice_sets(
            df = df,
            column_decider_id = column_decider_id,
            column_occasion_id = column_occasion_id,
            column_alternative = column_alternative
          )
          invisible(vapply(choice_sets$rows, function(rows) {
            private$.extract_long_choice(
              df = df,
              rows = rows,
              column_choice = column_choice,
              column_alternative = column_alternative
            )
          }, character(1)))
        }
        covariate_columns <- setdiff(
          names(df),
          c(
            column_decider_id,
            column_occasion_id,
            column_choice,
            column_alternative
          )
        )
        has_missing <- covariate_columns[vapply(
          df[covariate_columns],
          anyNA,
          logical(1)
        )]
        if (length(has_missing) > 0L) {
          warning(
            sprintf(
              "Covariate column(s) contain missing values: %s.",
              paste(sprintf("`%s`", has_missing), collapse = ", ")
            ),
            call. = FALSE
          )
        }
      }
      private$.df <- df
      private$.format <- format
      private$.column_decider_id <- column_decider_id
      private$.column_occasion_id <- column_occasion_id
      private$.column_choice <- column_choice
      private$.column_alternative <- column_alternative
      private$.update_metadata()
      invisible(self)
    },

    #' @description
    #' Simulate choice data. Covariates are drawn from a centered normal
    #' distribution.
    #'
    #' @return
    #' The simulated choice `data.frame`.
    #'
    #' It is also stored in the object together with `theta` in the fields
    #' `df` and `theta_dgp`.
    simulate = function(
      mod,
      theta = stats::rnorm(mod$length_theta),
      N = 100L,
      Tp = 1L,
      sd = 1
    ) {
      checkmate::assert_r6(mod, classes = "mod")
      checkmate::assert_numeric(
        theta,
        any.missing = FALSE,
        len = mod$length_theta
      )
      checkmate::assert_integerish(N, len = 1, lower = 1)
      checkmate::assert_integerish(Tp, lower = 1, any.missing = FALSE)
      checkmate::assert_number(sd, lower = 0, finite = TRUE)
      N <- as.integer(N)
      if (length(Tp) == 1L) {
        Tp <- rep(as.integer(Tp), N)
      } else if (length(Tp) == N) {
        Tp <- as.integer(Tp)
      } else {
        stop(
          "Argument `Tp` must be scalar or have length `N`.",
          call. = FALSE
        )
      }
      covariate_types <- read_formula_rhs(mod$formula)
      alt_names <- mod$alt_names
      id_column <- private$.column_decider_id
      if (is.null(id_column)) {
        id_column <- "id"
      }
      occasion_column <- private$.column_occasion_id

      rows <- sum(Tp)
      df <- data.frame(id_value = rep(seq_len(N), times = Tp))
      names(df) <- id_column
      if (is.null(occasion_column)) {
        occasion_column <- "occasion"
      }
      df[[occasion_column]] <- unlist(lapply(Tp, seq_len), use.names = FALSE)

      df[[private$.column_choice]] <- alt_names[1]
      for (var in covariate_types[[1]]) {
        for (alt_name in alt_names) {
          df[[paste0(var, "_", alt_name)]] <- stats::rnorm(rows, sd = sd)
        }
      }
      for (var in covariate_types[[2]]) {
        df[[var]] <- stats::rnorm(rows, sd = sd)
      }
      for (var in covariate_types[[3]]) {
        for (alt_name in alt_names) {
          df[[paste0(var, "_", alt_name)]] <- stats::rnorm(rows, sd = sd)
        }
      }
      self$redraw_choices(
        df = df,
        mod = mod,
        theta = theta,
        column_decider_id = id_column,
        column_occasion_id = occasion_column,
        column_choice = private$.column_choice
      )
    },

    #' @description
    #' Redraw choices for existing wide-format covariates.
    #'
    #' The `df` argument must be in wide format. It is recommended to use the
    #' output of `simulate()` and modify the covariates as desired.
    #'
    #' @return
    #' The choice `data.frame` with redrawn choices.
    #'
    #' It is also stored in the object together with `theta` in the fields
    #' `df` and `theta_dgp`.
    redraw_choices = function(
      df,
      mod,
      theta = self$theta_dgp,
      column_decider_id = NULL,
      column_occasion_id = NULL,
      column_choice = "choice"
    ) {
      checkmate::assert_r6(mod, classes = "mod")
      if (is.null(theta)) {
        stop(
          "Argument `theta` is required when `theta_dgp` is unavailable.",
          call. = FALSE
        )
      }
      checkmate::assert_numeric(
        theta,
        any.missing = FALSE,
        len = mod$length_theta
      )
      checkmate::assert_data_frame(df)
      df <- as.data.frame(df)
      if (nrow(df) == 0L) {
        stop("Argument `df` must contain at least one row.", call. = FALSE)
      }
      if (missing(column_decider_id)) {
        column_decider_id <- private$.column_decider_id
      }
      if (missing(column_occasion_id)) {
        column_occasion_id <- private$.column_occasion_id
      }
      if (missing(column_choice)) {
        column_choice <- private$.column_choice
      }
      checkmate::assert_string(
        column_decider_id,
        min.chars = 1,
        null.ok = TRUE
      )
      checkmate::assert_string(
        column_occasion_id,
        min.chars = 1,
        null.ok = TRUE
      )
      checkmate::assert_string(column_choice, min.chars = 1)
      if (is.null(column_decider_id) && !is.null(column_occasion_id)) {
        stop(
          "Argument `column_occasion_id` requires `column_decider_id`.",
          call. = FALSE
        )
      }
      missing_columns <- setdiff(
        c(
          if (!is.null(column_decider_id)) column_decider_id,
          if (!is.null(column_occasion_id)) column_occasion_id
        ),
        names(df)
      )
      if (length(missing_columns) > 0L) {
        stop(
          sprintf(
            "Required column(s) are missing: %s.",
            paste(sprintf("`%s`", missing_columns), collapse = ", ")
          ),
          call. = FALSE
        )
      }
      df[[column_choice]] <- mod$alt_names[1]
      data_sim <- choice_data$new(
        df = df,
        format = "wide",
        column_decider_id = column_decider_id,
        column_occasion_id = column_occasion_id,
        column_choice = column_choice
      )
      design <- data_sim$build_design(mod = mod, difference = FALSE)
      theta <- mod$make_theta_unique(theta)
      parameter <- mod$build_parameter(theta)
      b <- as.numeric(parameter$b)
      O <- parameter$Omega_chol
      L <- parameter$Sigma_chol
      random_count <- length(mod$var_names$random)
      random_index <- seq_along(b) > length(mod$var_names$fixed)
      if (!all(dim(O) == c(random_count, random_count))) {
        stop(
          "Argument `Omega` must match the random coefficient count.",
          call. = FALSE
        )
      }
      if (!all(dim(L) == c(mod$alt, mod$alt))) {
        stop(
          "Argument `Sigma` must match the number of alternatives.",
          call. = FALSE
        )
      }
      choice_rows <- if (is.null(column_decider_id)) {
        seq_len(nrow(df))
      } else {
        id_values <- unique(df[[column_decider_id]])
        unlist(
          lapply(
            id_values,
            function(id_value) which(df[[column_decider_id]] == id_value)
          ),
          use.names = FALSE
        )
      }
      row <- 1L
      for (n in seq_along(design)) {
        b_n <- b
        if (random_count > 0L) {
          b_n[random_index] <- b_n[random_index] +
            as.vector(O %*% stats::rnorm(random_count))
        }

        for (t in seq_along(design[[n]]$X)) {
          Xnt <- design[[n]]$X[[t]]
          error <- as.vector(L %*% stats::rnorm(mod$alt))
          V <- as.vector(Xnt %*% b_n) + error
          choice_index <- which.max(V)
          df[[column_choice]][choice_rows[row]] <- mod$alt_names[choice_index]
          row <- row + 1L
        }
      }
      private$.df <- df
      private$.format <- "wide"
      private$.column_decider_id <- column_decider_id
      private$.column_occasion_id <- column_occasion_id
      private$.column_choice <- column_choice
      private$.column_alternative <- NULL
      private$.theta_dgp <- theta
      private$.update_metadata()
      private$.df
    },

    #' @description
    #' Build model design matrices.
    #'
    #' @return
    #' A nested `list`. Each decider element contains `X`, a list of design
    #' matrices, and `y`, the chosen alternative indices in `mod$alt_names`
    #' order.
    build_design = function(mod, difference = FALSE, check_choice = TRUE) {
      checkmate::assert_r6(mod, classes = "mod")
      checkmate::assert_flag(difference)
      checkmate::assert_flag(check_choice)
      covariate_types <- read_formula_rhs(mod$formula)
      ASC <- read_formula_asc(mod$formula)
      alt_names <- mod$alt_names
      if (private$.format == "wide") {
        df <- private$.df
      } else {
        df <- private$.long_to_wide(
          df = private$.df,
          column_choice = private$.column_choice,
          covariate_types = covariate_types,
          column_decider_id = private$.column_decider_id,
          column_occasion_id = private$.column_occasion_id,
          alt_names = alt_names,
          column_alternative = private$.column_alternative
        )
      }
      if (nrow(df) > 0L) {
        missing_columns <- setdiff(private$.column_choice, names(df))
        if (length(missing_columns) > 0L) {
          stop("The choice column is missing.", call. = FALSE)
        }
        if (isTRUE(check_choice) && any(is.na(df[[private$.column_choice]]))) {
          warning(
            "Choice column contains missing values.",
            call. = FALSE
          )
        }
        choice_index <- match(
          as.character(df[[private$.column_choice]]), alt_names
        )
        if (isTRUE(check_choice) && any(is.na(choice_index))) {
          warning(
            "Choice column contains values outside `mod$alt_names`.",
            call. = FALSE
          )
        }
        needed <- unique(c(
          unlist(
            lapply(
              covariate_types[[1]],
              function(var) paste0(var, "_", alt_names)
            ),
            use.names = FALSE
          ),
          covariate_types[[2]],
          unlist(
            lapply(
              covariate_types[[3]],
              function(var) paste0(var, "_", alt_names)
            ),
            use.names = FALSE
          )
        ))
        missing_columns <- setdiff(needed, names(df))
        if (length(missing_columns) > 0L) {
          stop(
            sprintf(
              "Required covariate column(s) are missing: %s.",
              paste(sprintf("`%s`", missing_columns), collapse = ", ")
            ),
            call. = FALSE
          )
        }
        non_numeric <- needed[!vapply(df[needed], is.numeric, logical(1))]
        if (length(non_numeric) > 0L) {
          stop(
            sprintf(
              "Covariate column(s) must be numeric: %s.",
              paste(sprintf("`%s`", non_numeric), collapse = ", ")
            ),
            call. = FALSE
          )
        }
      }
      ref_index <- match(mod$ref_alt, alt_names)
      kept_alt_names <- alt_names[-ref_index]
      DN <- diag(mod$alt)[, -ref_index, drop = FALSE]
      choice_index <- match(
        as.character(df[[private$.column_choice]]), alt_names
      )
      id_values <- if (is.null(private$.column_decider_id)) {
        seq_len(nrow(df))
      } else {
        unique(df[[private$.column_decider_id]])
      }
      var_names <- unlist(mod$var_names, use.names = FALSE)
      out <- vector("list", length(id_values))
      Tp <- private$.Tp
      tot_choices <- 0

      for (n in seq_along(id_values)) {
        rows <- if (is.null(private$.column_decider_id)) {
          n
        } else {
          which(df[[private$.column_decider_id]] == id_values[n])
        }
        Xn <- vector("list", length(rows)-sum(is.na(choice_index[rows])))
        yn <- choice_index[rows]
        anz_t = 0 # count the number of decisions where y is not NA.
        for (t in seq_along(rows)) {
          row_index <- rows[t]
          fixed_blocks <- list()
          random_blocks <- list()
          for (var in covariate_types[[1]]) {
            columns <- paste0(var, "_", alt_names)
            values <- as.numeric(df[row_index, columns, drop = TRUE])
            block <- matrix(values, nrow = mod$alt, ncol = 1)
            colnames(block) <- var
            rownames(block) <- alt_names
            if (var %in% mod$re) {
              random_blocks[[length(random_blocks) + 1L]] <- block
            } else {
              fixed_blocks[[length(fixed_blocks) + 1L]] <- block
            }
          }
          for (var in covariate_types[[2]]) {
            block <- DN * df[[var]][row_index]
            colnames(block) <- paste0(var, "_", kept_alt_names)
            rownames(block) <- alt_names
            if (var %in% mod$re) {
              random_blocks[[length(random_blocks) + 1L]] <- block
            } else {
              fixed_blocks[[length(fixed_blocks) + 1L]] <- block
            }
          }
          if (ASC) {
            block <- DN
            colnames(block) <- paste0("ASC_", kept_alt_names)
            rownames(block) <- alt_names
            if ("ASC" %in% mod$re) {
              random_blocks[[length(random_blocks) + 1L]] <- block
            } else {
              fixed_blocks[[length(fixed_blocks) + 1L]] <- block
            }
          }
          for (var in covariate_types[[3]]) {
            columns <- paste0(var, "_", alt_names)
            values <- as.numeric(df[row_index, columns, drop = TRUE])
            block <- diag(values)
            colnames(block) <- paste0(var, "_", alt_names)
            rownames(block) <- alt_names
            if (var %in% mod$re) {
              random_blocks[[length(random_blocks) + 1L]] <- block
            } else {
              fixed_blocks[[length(fixed_blocks) + 1L]] <- block
            }
          }
          blocks <- c(fixed_blocks, random_blocks)
          X <- if (length(blocks) == 0L) {
            matrix(0, nrow = mod$alt, ncol = 0)
          } else {
            do.call(cbind, blocks)
          }
          rownames(X) <- alt_names
          colnames(X) <- var_names

          chosen_NA <- FALSE

          if (isTRUE(difference)) {
            chosen_index <- choice_index[row_index]

            if (is.na(chosen_index)){
              chosen_index = 1
              chosen_NA = TRUE
            } else {
              anz_t = anz_t + 1
            }
            chosen_row <- X[chosen_index, , drop = FALSE]
            X <- X - matrix(1, nrow(X), 1) %*% chosen_row
            X <- X[-chosen_index, , drop = FALSE]
            rownames(X) <- paste(
              alt_names[-chosen_index],
              alt_names[chosen_index],
              sep = " - "
            )
          } else {
            anz_t = anz_t +1
          }
          if (chosen_NA == FALSE){
            Xn[[anz_t]] <- X
            yn[anz_t] <- choice_index[row_index]
          }
        }
        Tp[n] <- anz_t
        tot_choices <- tot_choices + length(rows)-anz_t

        out[[n]] <- list(X = Xn, y = yn[1:anz_t])
      }
      private$.Tp <- as.integer(Tp)
      if (tot_choices >0L){
        warning(
            sprintf(
              "Deleted %d decisions due to NAs.", tot_choices
            ),
            call. = FALSE
        )
      }
      out
    },

    #' @description
    #' Print a summary of the choice data.
    #'
    #' @param what \[`character()`\]\cr
    #' Parts to print, any subset of `"header"`, `"data_format"`,
    #' `"dimensions"`, `"alternatives"`, and `"covariates"`.
    print = function(what = c(
      "header", "data_format", "dimensions", "alternatives", "covariates"
    )) {
      choices <- c(
        "header", "data_format", "dimensions", "alternatives", "covariates"
      )
      checkmate::assert_character(
        what, any.missing = FALSE, unique = TRUE
      )
      unknown <- setdiff(what, choices)
      if (length(unknown) > 0L) {
        stop(
          sprintf(
            "Argument `what` contains unknown value(s): %s.",
            paste(sprintf("`%s`", unknown), collapse = ", ")
          ),
          call. = FALSE
        )
      }
      if ("header" %in% what) {
        cat("<choice_data object>\n\n")
      }
      if ("data_format" %in% what) {
        cat(sprintf("Data format: %s\n", self$format))
      }
      if ("dimensions" %in% what) {
        cat(sprintf("Deciders: %d\n", self$N))
        if (length(self$Tp) > 0L) {
          cat(sprintf(
            "Choice occasions: from %d to %d\n",
            min(self$Tp),
            max(self$Tp)
          ))
        }
      }
      if ("alternatives" %in% what) {
        alternative_labels <- names(private$.choice_counts())
        cat(sprintf(
          "Alternatives: %s\n",
          if (length(alternative_labels) == 0L) {
            "none"
          } else {
            paste(alternative_labels, collapse = ", ")
          }
        ))
      }
      if ("covariates" %in% what) {
        covariates <- private$.covariate_constant_varying()
        cat(sprintf(
          "Alternative-constant covariates: %s\n",
          if (length(covariates$constant) == 0L) {
            "none"
          } else {
            paste(covariates$constant, collapse = ", ")
          }
        ))
        cat(sprintf(
          "Alternative-varying covariates: %s\n",
          if (length(covariates$varying) == 0L) {
            "none"
          } else {
            paste(covariates$varying, collapse = ", ")
          }
        ))
      }
      invisible(self)
    },

    #' @description
    #' Plot choice frequencies.
    #'
    #' @return
    #' A [ggplot2::ggplot()] object.
    plot = function() {
      values <- private$.choice_counts()
      if (length(values) == 0L) {
        stop("No choices available to plot.", call. = FALSE)
      }
      max_count <- max(values)
      if (max_count <= 10L) {
        count_breaks <- seq.int(0L, max_count, by = 1L)
      } else {
        count_breaks <- pretty(c(0L, max_count))
        count_breaks <- unique(as.integer(round(count_breaks)))
        count_breaks <- count_breaks[
          count_breaks >= 0L & count_breaks <= max_count
        ]
        count_breaks <- sort(unique(c(0L, count_breaks, max_count)))
      }
      plot_data <- data.frame(
        alternative = factor(names(values), levels = names(values)),
        value = as.numeric(values)
      )
      plot <- ggplot2::ggplot(
        data = plot_data,
        mapping = ggplot2::aes(
          x = !!as.name("alternative"),
          y = !!as.name("value")
        )
      )
      plot <- plot + ggplot2::geom_col()
      plot <- plot + ggplot2::scale_y_continuous(breaks = count_breaks)
      plot + ggplot2::labs(
        x = "Alternative",
        y = "Count",
        title = "Choice distribution"
      )
    }
  ),

  private = list(
    .df = data.frame(),
    .format = "wide",
    .column_decider_id = NULL,
    .column_occasion_id = NULL,
    .column_choice = "choice",
    .column_alternative = NULL,
    .N = 0L,
    .Tp = integer(),
    .id_values = NULL,
    .theta_dgp = NULL,

    # Identify alternative-constant and alternative-varying covariates.
    .covariate_constant_varying = function() {
      columns <- setdiff(
        names(private$.df),
        c(
          private$.column_decider_id,
          private$.column_occasion_id,
          private$.column_choice,
          private$.column_alternative
        )
      )
      covariates <- list(constant = character(), varying = character())
      if (length(columns) == 0L || nrow(private$.df) == 0L) {
        return(covariates)
      }
      if (private$.format == "wide") {
        alt_names <- names(private$.choice_counts())
        if (length(alt_names) == 0L) {
          covariates$constant <- columns
          return(covariates)
        }
        varying_columns <- character()
        for (column in columns) {
          for (alt_name in alt_names) {
            suffix <- paste0("_", alt_name)
            if (endsWith(column, suffix)) {
              base <- substr(column, 1L, nchar(column) - nchar(suffix))
              if (nzchar(base)) {
                candidate_columns <- paste0(base, "_", alt_names)
                if (all(candidate_columns %in% columns)) {
                  covariates$varying <- c(covariates$varying, base)
                  varying_columns <- c(varying_columns, candidate_columns)
                }
              }
            }
          }
        }
        covariates$varying <- unique(covariates$varying)
        covariates$constant <- setdiff(columns, unique(varying_columns))
        return(covariates)
      }
      choice_sets <- private$.group_long_choice_sets(
        df = private$.df,
        column_decider_id = private$.column_decider_id,
        column_occasion_id = private$.column_occasion_id,
        column_alternative = private$.column_alternative
      )
      is_constant <- vapply(
        columns,
        function(column) all(vapply(
          choice_sets$rows,
          function(rows) length(unique(na.omit(private$.df[[column]][rows]))) <= 1L,
          logical(1)
        )),
        logical(1)
      )
      covariates$constant <- columns[is_constant]
      covariates$varying <- columns[!is_constant]
      covariates
    },

    # Group long-format rows into choice sets.
    .group_long_choice_sets = function(
      df,
      column_decider_id,
      column_occasion_id,
      column_alternative,
      alt_names = NULL
    ) {
      checkmate::assert_data_frame(df)
      checkmate::assert_string(
        column_decider_id,
        min.chars = 1,
        null.ok = TRUE
      )
      checkmate::assert_string(
        column_occasion_id,
        min.chars = 1,
        null.ok = TRUE
      )
      checkmate::assert_string(column_alternative, min.chars = 1)
      columns <- c(
        if (!is.null(column_decider_id)) column_decider_id,
        if (!is.null(column_occasion_id)) column_occasion_id,
        column_alternative
      )
      missing_columns <- setdiff(columns, names(df))
      if (length(missing_columns) > 0L) {
        stop(
            sprintf(
              "Required long-data column(s) are missing: %s.",
            paste(sprintf("`%s`", missing_columns), collapse = ", ")
          ),
          call. = FALSE
        )
      }
      if (any(is.na(df[[column_alternative]]))) {
        stop(
          "Alternative column must not contain missing values.",
          call. = FALSE
        )
      }
      if (is.null(alt_names)) {
        alt_names <- unique(as.character(df[[column_alternative]]))
      } else {
        checkmate::assert_character(
          alt_names,
          any.missing = FALSE,
          min.len = 1,
          unique = TRUE
        )
        alt_names <- as.character(alt_names)
      }
      if (nrow(df) == 0L || length(alt_names) == 0L) {
        return(list(rows = list(), alt_names = alt_names))
      }
      if (!is.null(column_occasion_id)) {
        key_columns <- c(
          if (!is.null(column_decider_id)) column_decider_id,
          column_occasion_id
        )
        parts <- lapply(df[key_columns], function(x) {
          x <- as.character(x)
          paste0(nchar(x), ":", x)
        })
        key <- do.call(paste, c(parts, sep = "\r"))
        groups <- match(key, unique(key))
      } else if (is.null(column_decider_id)) {
        if (nrow(df) %% length(alt_names) != 0L) {
          stop(
            "Long data rows must form complete alternative blocks.",
            call. = FALSE
          )
        }
        groups <- rep(
          seq_len(nrow(df) / length(alt_names)),
          each = length(alt_names)
        )
      } else {
        groups <- integer(nrow(df))
        next_group <- 1L
        for (id_value in unique(df[[column_decider_id]])) {
          rows <- which(df[[column_decider_id]] == id_value)
          if (length(rows) %% length(alt_names) != 0L) {
            stop(
              "Long data rows must form complete alternative blocks.",
              call. = FALSE
            )
          }
          row_groups <- rep(
            seq_len(length(rows) / length(alt_names)),
            each = length(alt_names)
          )
          groups[rows] <- row_groups + next_group - 1L
          next_group <- max(groups[rows]) + 1L
        }
      }
      group_rows <- lapply(
        unique(groups),
        function(group) which(groups == group)
      )
      for (rows in group_rows) {
        alternatives <- as.character(df[[column_alternative]][rows])
        if (length(alternatives) != length(alt_names) ||
            anyDuplicated(alternatives) > 0L ||
            !setequal(alternatives, alt_names)) {
          stop(
            "Each long-data choice set must contain each alternative once.",
            call. = FALSE
          )
        }
      }
      list(rows = group_rows, alt_names = alt_names)
    },

    # Update metadata from the stored data.
    .update_metadata = function() {
      if (nrow(private$.df) == 0L) {
        private$.N <- 0L
        private$.Tp <- integer()
        private$.id_values <- NULL
        return(invisible(NULL))
      }
      if (private$.format == "wide") {
        if (is.null(private$.column_decider_id)) {
          private$.N <- nrow(private$.df)
          private$.Tp <- rep(1L, private$.N)
          private$.id_values <- seq_len(private$.N)
        } else {
          id_values <- unique(private$.df[[private$.column_decider_id]])
          private$.N <- length(id_values)
          if (is.null(private$.column_occasion_id)) {
            id_index <- match(
              private$.df[[private$.column_decider_id]],
              id_values
            )
            private$.Tp <- as.integer(tabulate(id_index, nbins = private$.N))
          } else {
            private$.Tp <- vapply(
              id_values,
              function(id_value) {
                rows <- private$.df[[private$.column_decider_id]] == id_value
                length(unique(private$.df[[private$.column_occasion_id]][rows]))
              },
              integer(1)
            )
          }
          private$.id_values <- id_values
        }
        return(invisible(NULL))
      }
      choice_sets <- private$.group_long_choice_sets(
        df = private$.df,
        column_decider_id = private$.column_decider_id,
        column_occasion_id = private$.column_occasion_id,
        column_alternative = private$.column_alternative
      )
      group_rows <- choice_sets$rows
      if (length(group_rows) == 0L) {
        private$.N <- 0L
        private$.Tp <- integer()
        private$.id_values <- NULL
        return(invisible(NULL))
      }
      if (is.null(private$.column_decider_id)) {
        private$.N <- length(group_rows)
        private$.Tp <- rep(1L, private$.N)
        private$.id_values <- seq_len(private$.N)
      } else {
        id_values <- unique(private$.df[[private$.column_decider_id]])
        group_id_index <- match(
          vapply(
            group_rows,
            function(rows) {
              as.character(private$.df[[private$.column_decider_id]][rows[1]])
            },
            character(1)
          ),
          as.character(id_values)
        )
        private$.N <- length(id_values)
        private$.Tp <- as.integer(
          tabulate(group_id_index, nbins = private$.N)
        )
        private$.id_values <- id_values
      }
      invisible(NULL)
    },

    # Convert long choice data to wide format.
    .long_to_wide = function(
      df,
      column_choice,
      covariate_types,
      column_decider_id,
      column_occasion_id,
      alt_names,
      column_alternative
    ) {
      checkmate::assert_string(
        column_alternative,
        min.chars = 1,
        null.ok = TRUE
      )
      checkmate::assert_string(
        column_occasion_id,
        min.chars = 1,
        null.ok = TRUE
      )
      missing_columns <- setdiff(column_alternative, names(df))
      if (length(missing_columns) > 0L) {
        stop(
            sprintf(
              "Required alternative column(s) are missing: %s.",
            paste(sprintf("`%s`", missing_columns), collapse = ", ")
          ),
          call. = FALSE
        )
      }
      long_columns <- unique(c(
        if (!is.null(column_decider_id)) column_decider_id,
        if (!is.null(column_occasion_id)) column_occasion_id,
        column_choice,
        column_alternative,
        unlist(covariate_types)
      ))
      missing_columns <- setdiff(long_columns, names(df))
      if (length(missing_columns) > 0L) {
        stop(
            sprintf(
              "Required long-data column(s) are missing: %s.",
            paste(sprintf("`%s`", missing_columns), collapse = ", ")
          ),
          call. = FALSE
        )
      }
      choice_sets <- private$.group_long_choice_sets(
        df = df,
        column_decider_id = column_decider_id,
        column_occasion_id = column_occasion_id,
        column_alternative = column_alternative,
        alt_names = alt_names
      )
      alt_names <- choice_sets$alt_names
      group_rows <- choice_sets$rows
      columns <- c(
        if (!is.null(column_decider_id)) column_decider_id,
        if (!is.null(column_occasion_id)) column_occasion_id,
        column_choice,
        covariate_types[[2]],
        unlist(
          lapply(
            covariate_types[[1]],
            function(var) paste0(var, "_", alt_names)
          ),
          use.names = FALSE
        ),
        unlist(
          lapply(
            covariate_types[[3]],
            function(var) paste0(var, "_", alt_names)
          ),
          use.names = FALSE
        )
      )
      wide_rows <- vector("list", length(group_rows))
      for (group in seq_along(group_rows)) {
        rows <- group_rows[[group]]
        alternatives <- as.character(df[[column_alternative]][rows])
        alt_rows <- rows[match(alt_names, alternatives)]
        row <- list()
        if (!is.null(column_decider_id)) {
          row[[column_decider_id]] <- df[[column_decider_id]][rows[1]]
        }
        if (!is.null(column_occasion_id)) {
          row[[column_occasion_id]] <- df[[column_occasion_id]][rows[1]]
        }
        row[[column_choice]] <- private$.extract_long_choice(
          df = df,
          rows = rows,
          column_choice = column_choice,
          column_alternative = column_alternative
        )
        for (var in covariate_types[[2]]) {
          values <- df[[var]][rows]
          #if (any(is.na(values))) {
          #  row[[var]] <- NA
          #} else {
            unique_values <- unique(na.omit(values))
            if (length(unique_values) != 1L) {
              stop(
                sprintf(
                  paste(
                    "Long covariate `%s` must be constant within each",
                    "choice set."
                  ),
                  var
                ),
                call. = FALSE
              )
            }
            row[[var]] <- unique_values[1]
          #}
        }
        for (var in covariate_types[[1]]) {
          row[paste0(var, "_", alt_names)] <- as.list(df[[var]][alt_rows])
        }
        for (var in covariate_types[[3]]) {
          row[paste0(var, "_", alt_names)] <- as.list(df[[var]][alt_rows])
        }
        wide_rows[[group]] <- as.data.frame(
          row[columns],
          stringsAsFactors = FALSE
        )
      }
      wide <- do.call(rbind, wide_rows)
      rownames(wide) <- NULL
      wide
    },

    # Extract the chosen alternative from long-format choice set.
    .extract_long_choice = function(
      df,
      rows,
      column_choice,
      column_alternative
    ) {
      choices <- df[[column_choice]][rows]
      alternatives <- as.character(df[[column_alternative]][rows])
      if (any(is.na(choices))) {
       # warning(
       #   "Choice column contains missing values.",
       #   call. = FALSE
       # )
        return(c("NA"))
      }
      if (is.logical(choices)) {
        index <- which(!is.na(choices) & choices)
        if (length(index) == 1L) {
          return(alternatives[index])
        }
      }
      if (is.numeric(choices)) {
        values <- choices[!is.na(choices)]
        if (length(values) > 0L && all(values %in% c(0, 1))) {
          index <- which(!is.na(choices) & choices == 1)
          if (length(index) == 1L) {
            return(alternatives[index])
          }
        }
      }
      stop(
        "Each long-data choice set must identify exactly one choice.",
        call. = FALSE
      )
    },

    # Count choices in the stored data.
    .choice_counts = function() {
      if (nrow(private$.df) == 0L) {
        return(numeric())
      }
      if (private$.format == "wide") {
        choice_values <- as.character(private$.df[[private$.column_choice]])
        choice_values <- choice_values[!is.na(choice_values)]
        levels <- sort(unique(choice_values))
      } else {
        levels <- unique(as.character(
          private$.df[[private$.column_alternative]]
        ))
        choice_sets <- private$.group_long_choice_sets(
          df = private$.df,
          column_decider_id = private$.column_decider_id,
          column_occasion_id = private$.column_occasion_id,
          column_alternative = private$.column_alternative
        )
        choice_values <- vapply(choice_sets$rows, function(rows) {
          private$.extract_long_choice(
            df = private$.df,
            rows = rows,
            column_choice = private$.column_choice,
            column_alternative = private$.column_alternative
          )
        }, character(1))
      }
      choice_values <- choice_values[!is.na(choice_values)]
      counts <- table(factor(choice_values, levels = levels))
      stats::setNames(as.numeric(counts), names(counts))
    }
  ),

  active = list(

    #' @field df (read only) \[`data.frame`\]\cr
    #' Choice data.
    df = function(value) {
      if (missing(value)) {
        private$.df
      } else {
        stop("Field `df` is read-only.", call. = FALSE)
      }
    },

    #' @field format (read only) \[`character(1)`\]\cr
    #' Data format, either `"wide"` or `"long"`.
    format = function(value) {
      if (missing(value)) {
        private$.format
      } else {
        stop("Field `format` is read-only.", call. = FALSE)
      }
    },

    #' @field alt_names (read only) \[`character()`\]\cr
    #' Alternative labels.
    alt_names = function(value) {
      if (missing(value)) {
        labels <- names(private$.choice_counts())
        if (is.null(labels)) {
          character()
        } else {
          labels
        }
      } else {
        stop("Field `alt_names` is read-only.", call. = FALSE)
      }
    },

    #' @field N (read only) \[`integer(1)`\]\cr
    #' Number of deciders.
    N = function(value) {
      if (missing(value)) {
        private$.N
      } else {
        stop("Field `N` is read-only.", call. = FALSE)
      }
    },

    #' @field Tp (read only) \[`integer()`\]\cr
    #' Number of choice occasions per decider.
    Tp = function(value) {
      if (missing(value)) {
        private$.Tp
      } else {
        stop("Field `Tp` is read-only.", call. = FALSE)
      }
    },

    #' @field theta_dgp (read only) \[`numeric()` | `NULL`\]\cr
    #' Free parameter vector used for simulating the current data. `NULL` if
    #' the object does not contain simulated data.
    theta_dgp = function(value) {
      if (missing(value)) {
        private$.theta_dgp
      } else {
        stop("Field `theta_dgp` is read-only.", call. = FALSE)
      }
    },

    #' @field choice_counts (read only) \[`table`\]\cr
    #' Absolute choice counts.
    choice_counts = function(value) {
      if (missing(value)) {
        as.table(private$.choice_counts())
      } else {
        stop("Field `choice_counts` is read-only.", call. = FALSE)
      }
    },

    #' @field choice_frequencies (read only) \[`table`\]\cr
    #' Relative choice frequencies.
    choice_frequencies = function(value) {
      if (missing(value)) {
        count_table <- self$choice_counts
        counts <- as.numeric(count_table)
        names(counts) <- names(count_table)
        if (length(counts) == 0L) {
          return(as.table(numeric()))
        }
        as.table(counts / sum(counts))
      } else {
        stop("Field `choice_frequencies` is read-only.", call. = FALSE)
      }
    }
  )
)
