#' Model Specification
#'
#' @description
#' An object of the `mod` R6 class defines the probit model structure and
#' stores linear restrictions on the parameter blocks \eqn{b}, \eqn{\Omega}, and
#' \eqn{\Sigma} from the utility system
#' \deqn{
#' U_n = X_n \begin{bmatrix} b_f \\ \beta_{r,n} \end{bmatrix} + \epsilon_n,
#' \quad \beta_{r,n} \sim N(b_r, \Omega),
#' \quad b = \begin{bmatrix} b_f \\ b_r \end{bmatrix},
#' \quad \epsilon_n \sim N(0, \Sigma),
#' }
#' which are represented via a matrix \eqn{H} and a column vector \eqn{f} for
#' each parameter block as follows:
#' \deqn{b = H_b \theta_b + f_b,}
#' \deqn{O = \mathrm{vech}(\mathrm{chol}(\Omega)) = H_O \theta_O + f_O,}
#' \deqn{L = \mathrm{vech}(\mathrm{chol}(\Sigma)) = H_L \theta_L + f_L.}
#'
#' - A `mod` object can be created via
#'
#'   ```r
#'   m <- mod$new(
#'     formula = choice ~ price | 1 | time,  # model formula
#'     re = c("price", "ASC"),               # random effects
#'     alt = 3,                              # number of alternatives
#'     alt_names = c("train", "bus", "car"), # alternative labels (optional)
#'     ref_alt = "train",                    # reference alternative
#'     corr_re = TRUE,                       # random effects are correlated
#'     error_struc = "diag"                  # diagonal error covariance matrix
#'   )
#'   ```
#'   with the inputs
#'
#'   - `formula`, a symbolic description of the choice model, see details,
#'   - `re`, the covariates to be connected to random coefficients,
#'   - `alt`, the number of choice alternatives,
#'   - `alt_names`, labels for the choice alternatives,
#'   - `ref_alt`, the reference alternative,
#'   - and `corr_re` and `error_struc` can be used for common \eqn{\Omega} and
#'     \eqn{\Sigma} structures, see details.
#'
#' - The source of truth for the model structure is `formula`, `re`,
#'   `alt_names`, and `ref_alt`. These can be changed via their active bindings.
#'   The dependent quantities `alt`, `b_count`, `var_names`, `length_theta_b`,
#'   `length_theta_O`, `length_theta_L`, and `length_theta` are read-only and
#'   are updated internally, together with the dimensions of the \eqn{H}
#'   matrices and \eqn{f} vectors. The \eqn{H} matrices and \eqn{f} vectors can
#'   be accessed (for example, `m$f_L` for the vector \eqn{f_L}) and modified
#'   within the dimensions implied by `formula`, `re`, `alt_names`, and
#'   `ref_alt`.
#'
#' - Warning messages can be suppressed by the active binding
#'   `m$suppress_warnings`. Its default is `FALSE`.
#'
#' - The vector \eqn{\theta = \begin{bmatrix} \theta_b & \theta_O & \theta_L
#'   \end{bmatrix}^\top} over which the model is fitted to data can be
#'   transformed to the parameters \eqn{b, \Omega, \Sigma} via
#'   `m$build_parameter()`. The method `m$build_theta()` reverses this task.
#'
#' - Restrictions on the model parameters are used to identify a model with
#'   respect to level and scale. By default, if `error_struc = "full"`, the
#'   first column and row of \eqn{\Sigma} is fixed to zero and the
#'   \eqn{(2,2)}-entry \eqn{\Sigma_{2,2}} is fixed to one. The method
#'   `m$check_identifiability()` checks whether a specification seems
#'   identifiable.
#'
#' @inheritSection read_formula Specifying the model formula
#'
#' @section Common covariance structures:
#' The covariance matrices \eqn{\Omega} and \eqn{\Sigma} are parameterized via
#' the lower-triangular Cholesky factors \eqn{O} and \eqn{L}.
#'
#' For the \strong{random-effect covariance matrix \eqn{\Omega}}, the argument
#' `corr_re` controls whether correlations between random effects are included:
#'
#' - If `corr_re = FALSE`, \eqn{\Omega} is diagonal.
#' - If `corr_re = TRUE`, \eqn{\Omega} is fully flexible.
#'
#' For the \strong{error covariance matrix \eqn{\Sigma}}, `error_struc` selects
#' one of the following common structures:
#'
#' - `"full"` creates a flexible covariance matrix. The entry `L[2,2]` of the
#'   Cholesky factor is fixed to one for scale identification, the first
#'   Cholesky column is fixed to zero, and all remaining Cholesky entries are
#'   free.
#' - `"diag"` creates a diagonal covariance matrix. The first Cholesky diagonal
#'   entry is fixed to one, and the remaining diagonal entries are free.
#' - `"ident"` creates a covariance matrix proportional to the identity. All
#'   Cholesky diagonal entries share one common free parameter.
#' - `"fixed"` fixes the covariance matrix to the identity, i.e. all Cholesky
#'   diagonal entries are fixed to one.
#'
#' @param formula \[`formula`\]\cr
#' A symbolic description of the choice model.
#'
#' @param re \[`character()`\]\cr
#' Covariates with random coefficients.
#'
#' @param alt \[`integer(1)`\]\cr
#' Number of choice alternatives.
#'
#' @param alt_names \[`character()`\]\cr
#' Alternative labels. If supplied, its length must match `alt`.
#'
#' @param ref_alt \[`character(1)`\]\cr
#' Reference alternative. Must match one of `alt_names`.
#'
#' @param corr_re \[`logical(1)`\]\cr
#' Whether random effects are correlated.
#'
#' @param error_struc \[`character(1)`\]\cr
#' Error structure, one of:
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
#' @param ... Optional restriction matrices or fixed offset vectors named
#' `H_b`, `f_b`, `H_O`, `f_O`, `H_L`, and `f_L`. These must have the row
#' dimensions implied by `formula`, `re`, `alt_names`, and `ref_alt`; all
#' user-supplied dimnames are ignored and removed on assignment.
#'
#' @param theta \[`numeric`\]\cr
#' A numeric vector containing the free model parameters in the order
#' \eqn{\theta = \begin{bmatrix} \theta_b & \theta_O & \theta_L
#' \end{bmatrix}^\top}.
#'
#' @param b \[`numeric` | `matrix`\]\cr
#' Coefficient vector.
#'
#' @param Omega \[`matrix` | `NULL`\]\cr
#' Random-effect covariance matrix.
#'
#' Can be `NULL` if there are no random effects.
#'
#' @param Sigma \[`matrix`\]\cr
#' Error-term covariance matrix.
#'
#' @param tolerance \[`numeric(1)`\]\cr
#' A positive numeric value defining the check tolerance.
#'
#' @examples
#' # define model structure:
#' m <- mod$new(
#'   formula = choice ~ price | 1 | time,
#'   re = c("price", "ASC"),
#'   alt = 3,
#'   alt_names = c("train", "bus", "car"),
#'   ref_alt = "train",
#'   corr_re = TRUE,
#'   error_struc = "diag"
#' )
#' print(m)
#'
#' # transform parameter vector to b, Omega, and Sigma:
#' theta <- runif(m$length_theta)
#' parameter <- m$build_parameter(theta)
#' parameter$b
#' parameter$Omega
#' parameter$Sigma
#'
#' # reverse this operation
#' theta2 <- m$build_theta(
#'   b = parameter$b,
#'   Omega = parameter$Omega,
#'   Sigma = parameter$Sigma
#' )
#' all.equal(theta, theta2)
#'
#' # is the model specification identified?
#' m$check_identifiability()
#' m$H_L <- diag(6)[, 4:6]; m$f_L[4, ] <- 0 # removed scale identification
#' m$check_identifiability()
#'
#' @export

mod <- R6::R6Class(
  classname = "mod",
  public = list(

    #' @description
    #' Create a new `mod` object.
    #'
    #' @return
    #' A new `mod` object.
    initialize = function(
      formula,
      re = character(),
      alt = 2,
      alt_names = LETTERS[seq_len(alt)],
      ref_alt = alt_names[1],
      corr_re = FALSE,
      error_struc = "full",
      ...
    ) {

      ### read ... arguments
      restrictions <- list(...)
      restriction_names <- names(restrictions)
      allowed_restrictions <- c("H_b", "f_b", "H_O", "f_O", "H_L", "f_L")
      if (
        length(restrictions) > 0 &&
          (is.null(restriction_names) || any(restriction_names == ""))
      ) {
        stop("All `...` arguments must be named.", call. = FALSE)
      }
      unknown_restrictions <- setdiff(restriction_names, allowed_restrictions)
      if (length(unknown_restrictions) > 0) {
        stop(
          paste0(
            "The following `...` argument(s) are unknown: ",
            paste(sprintf("`%s`", unknown_restrictions), collapse = ", "),
            "."
          ),
          call. = FALSE
        )
      }

      ### check inputs
      if (missing(formula)) {
        stop("Argument `formula` must be specified.", call. = FALSE)
      }
      checkmate::assert_int(alt, lower = 2)
      alt <- as.integer(alt)
      if (missing(alt_names)) {
        alt_names <- LETTERS[seq_len(alt)]
      }
      covariate_types <- read_formula_rhs(formula)
      ASC <- read_formula_asc(formula)
      checkmate::assert_character(
        alt_names, any.missing = FALSE, unique = TRUE
      )
      alt_names <- as.character(alt_names)
      if (length(alt_names) != alt) {
        stop("Argument `alt_names` must have length `alt`.", call. = FALSE)
      }
      checkmate::assert_string(ref_alt, min.chars = 1)
      ref_alt_index <- match(ref_alt, alt_names)
      if (is.na(ref_alt_index)) {
        stop("Argument `ref_alt` is not in `alt_names`.", call. = FALSE)
      }
      ref_alt_index <- as.integer(ref_alt_index)
      checkmate::assert_character(re, any.missing = FALSE)
      re <- unique(re)
      private$.check_re(covariate_types, ASC, re)
      checkmate::assert_flag(corr_re)
      checkmate::assert_choice(
        error_struc, choices = c("diag", "fixed", "ident", "full")
      )

      ### define (and sync) H/f (w/o warnings)
      old_suppress_warnings <- private$.suppress_warnings
      private$.suppress_warnings <- TRUE
      on.exit(
        private$.suppress_warnings <- old_suppress_warnings,
        add = TRUE
      )
      private$.formula <- formula
      private$.covariate_types <- covariate_types
      private$.ASC <- ASC
      private$.re <- re
      private$.alt <- alt
      private$.alt_names <- alt_names
      private$.ref_alt <- ref_alt
      private$.ref_alt_index <- ref_alt_index
      b_count <- self$b_count
      Hf_b <- private$.setup_b(sum(b_count))
      Hf_Omega <- private$.setup_Omega(b_count[2], corr_re)
      Hf_Sigma <- private$.setup_Sigma(self$alt, error_struc)
      self$H_b <- Hf_b$H
      self$f_b <- Hf_b$f
      self$H_O <- Hf_Omega$H
      self$f_O <- Hf_Omega$f
      self$H_L <- Hf_Sigma$H
      self$f_L <- Hf_Sigma$f
      for (restriction in intersect(allowed_restrictions, restriction_names)) {
        self[[restriction]] <- restrictions[[restriction]]
      }
    },

    #' @description
    #' Transform a free parameter vector to model parameters.
    #'
    #' @return
    #' A named list with the entries `b`, `O`, `Omega_chol`, `Omega`, `L`,
    #' `Sigma_chol`, and `Sigma`.
    build_parameter = function(theta = stats::runif(self$length_theta)) {
      checkmate::assert_numeric(
        theta, any.missing = FALSE, len = self$length_theta
      )
      var_names <- self$var_names
      H_b <- self$H_b
      f_b <- self$f_b
      if (self$length_theta_b > 0) {
        b <- H_b %*% theta[seq_len(self$length_theta_b)] + f_b
      } else {
        b <- f_b
      }
      rownames(b) <- unlist(var_names, use.names = FALSE)
      colnames(b) <- " "
      H_O <- self$H_O
      f_O <- self$f_O
      if (self$length_theta_O > 0) {
        omega_index <- self$length_theta_b + seq_len(self$length_theta_O)
        O <- H_O %*% theta[omega_index] + f_O
      } else {
        O <- f_O
      }
      random_names <- var_names$random
      rownames(O) <- private$.vech_labels(random_names)
      colnames(O) <- " "
      Omega_chol <- vech_chol_to_chol(O)
      dimnames(Omega_chol) <- list(random_names, random_names)
      Omega <- Omega_chol %*% t(Omega_chol)
      dimnames(Omega) <- list(random_names, random_names)
      H_L <- self$H_L
      f_L <- self$f_L
      if (self$length_theta_L > 0) {
        sigma_index <- self$length_theta_b +
          self$length_theta_O +
          seq_len(self$length_theta_L)
        L <- H_L %*% theta[sigma_index] + f_L
      } else {
        L <- f_L
      }
      alt_names <- self$alt_names
      rownames(L) <- private$.vech_labels(alt_names)
      colnames(L) <- " "
      Sigma_chol <- vech_chol_to_chol(L)
      dimnames(Sigma_chol) <- list(alt_names, alt_names)
      Sigma <- Sigma_chol %*% t(Sigma_chol)
      dimnames(Sigma) <- list(alt_names, alt_names)
      list(
        "b" = b,
        "O" = O,
        "Omega_chol" = Omega_chol,
        "Omega" = Omega,
        "L" = L,
        "Sigma_chol" = Sigma_chol,
        "Sigma" = Sigma
      )
    },

    #' @description
    #' Transform model parameters to a free parameter vector.
    #'
    #' If the supplied parameters cannot be represented exactly under the
    #' current restrictions, the returned vector is a least-squares
    #' approximation.
    #'
    #' @return
    #' A numeric vector `theta`.
    build_theta = function(b, Omega = NULL, Sigma) {
      b_count <- self$b_count
      if (is.matrix(b)) {
        if (ncol(b) != 1) {
          stop("Argument `b` must have one column.", call. = FALSE)
        }
      } else if (is.numeric(b)) {
        b <- matrix(b, ncol = 1)
      } else {
        stop(
          "Argument `b` must be a matrix or numeric vector.",
          call. = FALSE
        )
      }
      checkmate::assert_matrix(
        b, mode = "numeric", any.missing = FALSE, nrows = sum(b_count)
      )
      size <- b_count[2]
      if (is.null(Omega)) {
        if (size > 0) {
          stop(
            "Argument `Omega` is required when `mod` has random effects.",
            call. = FALSE
          )
        }
        Omega <- matrix(0, 0, 0)
      } else {
        if (!check_covariance_matrix(Omega)) {
          stop("Argument `Omega` must be a covariance matrix.", call. = FALSE)
        }
        if (nrow(Omega) != size) {
          stop(
            "Argument `Omega` has dimension ", nrow(Omega), " x ", ncol(Omega),
            ", but `mod` expects a ", size, " x ", size,
            " covariance matrix to match the random coefficient count.",
            call. = FALSE
          )
        }
      }
      tolerance <- sqrt(.Machine$double.eps)
      Omega_chol <- matrix(0, size, size)
      if (size > 0) {
        if (any(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values <=
                tolerance)) {
          stop("Argument `Omega` must be positive definite.", call. = FALSE)
        }
        Omega_chol <- t(chol(Omega))
      }
      O <- matrix(Omega_chol[lower.tri(Omega_chol, diag = TRUE)], ncol = 1)
      size <- self$alt
      if (!check_covariance_matrix(Sigma)) {
        stop("Argument `Sigma` must be a covariance matrix.", call. = FALSE)
      }
      if (nrow(Sigma) != size) {
        stop(
          "Argument `Sigma` has dimension ", nrow(Sigma), " x ", ncol(Sigma),
          ", but `mod` expects a ", size, " x ", size,
          " covariance matrix to match `mod$alt_names`.",
          call. = FALSE
        )
      }
      Sigma_chol <- matrix(0, size, size)
      if (size > 0) {
        positive <- diag(Sigma) > tolerance
        if (any(!positive) &&
            any(Sigma[!positive, , drop = FALSE] != 0)) {
          stop("Argument `Sigma` must be positive definite.", call. = FALSE)
        }
        if (any(positive)) {
          Sigma_positive <- Sigma[positive, positive, drop = FALSE]
          Sigma_chol_positive <- tryCatch(
            t(chol(Sigma_positive)),
            error = function(error) NULL
          )
          if (is.null(Sigma_chol_positive)) {
            stop("Argument `Sigma` must be positive definite.", call. = FALSE)
          }
          Sigma_chol[positive, positive] <- Sigma_chol_positive
        }
      }
      L <- matrix(Sigma_chol[lower.tri(Sigma_chol, diag = TRUE)], ncol = 1)
      theta <- numeric()
      used_approximation <- FALSE
      theta_blocks <- list(
        list(H = self$H_b, y = b - self$f_b),
        list(H = self$H_O, y = O - self$f_O),
        list(H = self$H_L, y = L - self$f_L)
      )
      for (block in theta_blocks) {
        H <- block$H
        y <- block$y
        if (ncol(H) > 0) {
          if (nrow(H) == 0) {
            theta <- c(theta, rep(0, ncol(H)))
          } else {
            theta_block <- stats::lm.fit(x = H, y = y)$coefficients
            theta_block[is.na(theta_block)] <- 0
            theta_block <- as.vector(theta_block)
            fitted <- H %*% theta_block
            used_approximation <- used_approximation ||
              !isTRUE(all.equal(
                as.vector(fitted), as.vector(y), check.attributes = FALSE
              ))
            theta <- c(theta, theta_block)
          }
        } else if (nrow(H) > 0) {
          used_approximation <- used_approximation || max(abs(y)) > 0
        }
      }
      if (used_approximation && !self$suppress_warnings) {
        warning(
          paste0(
            "The supplied parameters cannot be represented exactly by `mod`; ",
            "`theta` is a least-squares approximation."
          ),
          call. = FALSE
        )
      }
      theta
    },

    #' @description
    #' Replace arbitrary Cholesky entry signs with a unique representation
    #' using positive diagonal entries.
    #'
    #' @return
    #' A numeric vector representing the same model parameters as `theta`.
    make_theta_unique = function(theta) {
      checkmate::assert_numeric(
        theta, any.missing = FALSE, finite = TRUE, len = self$length_theta
      )
      parameter <- self$build_parameter(theta)
      unique_theta <- tryCatch(
        suppressWarnings(self$build_theta(
          b = parameter$b,
          Omega = parameter$Omega,
          Sigma = parameter$Sigma
        )),
        error = function(error) NULL
      )
      # keep theta if positive-root representation cannot be constructed
      if (is.null(unique_theta)) return(theta)
      # verify that custom restrictions did not alter any model parameter
      unique_parameter <- self$build_parameter(unique_theta)
      same_parameters <- isTRUE(all.equal(
        c(parameter$b, parameter$Omega, parameter$Sigma),
        c(
          unique_parameter$b,
          unique_parameter$Omega,
          unique_parameter$Sigma
        ),
        check.attributes = FALSE
      ))
      if (same_parameters) unique_theta else theta
    },

    #' @description
    #' Check whether the model specification seems identifiable.
    #'
    #' @return
    #' `TRUE` if the numerical rank check indicates that the model is
    #' identifiable, otherwise `FALSE`.
    check_identifiability = function(tolerance = 1e-6) {
      checkmate::assert_number(tolerance, finite = TRUE)
      if (tolerance <= 0) {
        stop("Argument `tolerance` must be positive.", call. = FALSE)
      }
      npar <- self$length_theta
      if (npar == 0) return(TRUE)
      theta <- stats::rnorm(npar)
      parameter <- self$build_parameter(theta)
      b <- parameter$b
      norm_b <- norm(b)
      Omega <- parameter$Omega
      Sigma <- parameter$Sigma
      alt <- self$alt
      delta_1 <- diag(alt)
      delta_1[, 1] <- -1
      M_1 <- matrix(delta_1[seq.int(2L, alt), ], nrow = alt - 1L)
      Sigmac <- M_1 %*% Sigma %*% t(M_1)
      if (norm_b > 0) {
        b <- b / norm_b
        Omega <- Omega / norm_b^2
        Sigmac <- Sigmac / norm_b^2
      } else {
        norm_Sigmac <- norm(Sigmac)
        if (norm_Sigmac > 0) {
          Omega <- Omega / norm_Sigmac
          Sigmac <- Sigmac / norm_Sigmac
        }
      }
      vech_0 <- c(as.vector(b), as.vector(Omega), as.vector(Sigmac))
      dims <- length(vech_0)
      systems <- matrix(0, dims, npar)
      for (i in seq_len(npar)) {
        theta_shift <- theta
        theta_shift[i] <- theta_shift[i] + tolerance
        parameter <- self$build_parameter(theta_shift)
        b <- parameter$b
        norm_b <- norm(b)
        Omega <- parameter$Omega
        Sigma <- parameter$Sigma
        Sigmac <- M_1 %*% Sigma %*% t(M_1)
        if (norm_b > 0) {
          b <- b / norm_b
          Omega <- Omega / norm_b^2
          Sigmac <- Sigmac / norm_b^2
        } else {
          norm_Sigmac <- norm(Sigmac)
          if (norm_Sigmac > 0) {
            Omega <- Omega / norm_Sigmac
            Sigmac <- Sigmac / norm_Sigmac
          }
        }
        vech_j <- c(as.vector(b), as.vector(Omega), as.vector(Sigmac))
        systems[, i] <- (vech_j - vech_0) / tolerance
      }
      if (qr(systems, tol = tolerance)$rank < npar) {
        if (!self$suppress_warnings) {
          warning(
            "The parameter set may not be identified; check the specification.",
            call. = FALSE
          )
        }
        return(FALSE)
      }
      return(TRUE)
    },

    #' @description
    #' Print a summary of the model specification.
    #'
    #' @param what \[`character()`\]\cr
    #' Parts to print, any subset of `"header"`, `"formula"`, `"re"`,
    #' `"alt_names"`, `"b_count"`, `"length_theta"`, `"utility_equations"`,
    #' and `"restrictions"`.
    #'
    #' @return
    #' The `mod` object, invisibly. Called for its console output side effect.
    print = function(what = c(
      "header", "formula", "re", "alt_names", "b_count", "length_theta",
      "utility_equations", "restrictions"
    )) {
      choices <- c(
        "header", "formula", "re", "alt_names", "b_count", "length_theta",
        "utility_equations", "restrictions"
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
        cat("<mod object>\n\n")
      }
      if ("formula" %in% what) {
        cat(sprintf("Formula: %s\n", deparse1(self$formula)))
      }
      if ("re" %in% what) {
        cat(sprintf(
          "Random effects: %s\n",
          if (length(self$re) > 0) paste(self$re, collapse = ", ") else "none"
        ))
      }
      if ("alt_names" %in% what) {
        cat(
          sprintf(
            "Alternative names: %s\n",
            paste(self$alt_names, collapse = ", ")
          )
        )
      }
      if ("b_count" %in% what) {
        b_count <- self$b_count
        cat(sprintf(
          "Coefficients: %d fixed, %d random\n",
          b_count[1],
          b_count[2]
        ))
      }
      if ("length_theta" %in% what) {
        cat(
          sprintf(
            "Free parameters total: %d\n",
            self$length_theta
          )
        )
      }
      if (
        any(
          c("formula", "re", "alt_names", "b_count", "length_theta") %in% what
        ) &&
          any(c("utility_equations", "restrictions") %in% what)
      ) {
        cat("\n")
      }
      if ("utility_equations" %in% what) {
        private$.print_utility_equations()
      }
      if ("restrictions" %in% what) {
        private$.print_restrictions()
      }
      invisible(self)
    }

  ),
  private = list(

    .formula = NULL,
    .covariate_types = list(character(), character(), character()),
    .ASC = TRUE,
    .re = character(),
    .alt = 2L,
    .alt_names = character(),
    .ref_alt = character(),
    .ref_alt_index = 1L,
    .length_theta_b = 0,
    .length_theta_O = 0,
    .length_theta_L = 0,
    .f_b = matrix(0, 0, 0),
    .H_b = matrix(0, 0, 0),
    .H_O = matrix(0, 0, 0),
    .f_O = matrix(0, 0, 0),
    .H_L = matrix(0, 0, 0),
    .f_L = matrix(0, 0, 0),
    .suppress_warnings = FALSE,

    # Print utility equations using the current model specification.
    .print_utility_equations = function() {
      b_names <- unlist(self$var_names, use.names = FALSE)
      b_index <- setNames(seq_along(b_names), b_names)
      cat("Utility equations:\n\n")
      for (alt_index in seq_along(self$alt_names)) {
        alt_name <- self$alt_names[alt_index]
        term_index <- integer()
        term_text <- character()
        for (var in private$.covariate_types[[1]]) {
          index <- b_index[[var]]
          covariate_name <- paste0(var, "_", alt_name)
          term_index <- c(term_index, index)
          term_text <- c(
            term_text,
            sprintf("b[%d] * %s", index, covariate_name)
          )
        }
        if (!identical(alt_index, private$.ref_alt_index)) {
          for (var in private$.covariate_types[[2]]) {
            coefficient_name <- paste0(var, "_", alt_name)
            index <- b_index[[coefficient_name]]
            term_index <- c(term_index, index)
            term_text <- c(term_text, sprintf("b[%d] * %s", index, var))
          }
          if (private$.ASC) {
            coefficient_name <- paste0("ASC_", alt_name)
            index <- b_index[[coefficient_name]]
            term_index <- c(term_index, index)
            term_text <- c(term_text, sprintf("b[%d]", index))
          }
        }
        for (var in private$.covariate_types[[3]]) {
          coefficient_name <- paste0(var, "_", alt_name)
          index <- b_index[[coefficient_name]]
          term_index <- c(term_index, index)
          term_text <- c(
            term_text,
            sprintf("b[%d] * %s", index, coefficient_name)
          )
        }
        terms <- term_text[order(term_index)]
        terms <- c(terms, sprintf("eps[%d]", alt_index))
        cat(sprintf("  U_%s = %s\n", alt_name, paste(terms, collapse = " + ")))
      }
      cat("\n")
    },

    # Print linear restriction matrices.
    .print_restrictions = function() {
      restriction_blocks <- list(
        list(
          parameter_name = "b",
          definition_full = "b = H_b theta_b + f_b",
          definition_H = "b = H_b theta_b",
          definition_f = "b = f_b",
          theta_name = "theta_b",
          H_name = "H_b",
          f_name = "f_b",
          H = self$H_b,
          f = self$f_b,
          n_theta = self$length_theta_b
        ),
        list(
          parameter_name = "O",
          definition_full = "O = vech(chol(Omega)) = H_O theta_O + f_O",
          definition_H = "O = vech(chol(Omega)) = H_O theta_O",
          definition_f = "O = vech(chol(Omega)) = f_O",
          theta_name = "theta_O",
          H_name = "H_O",
          f_name = "f_O",
          H = self$H_O,
          f = self$f_O,
          n_theta = self$length_theta_O
        ),
        list(
          parameter_name = "L",
          definition_full = "L = vech(chol(Sigma)) = H_L theta_L + f_L",
          definition_H = "L = vech(chol(Sigma)) = H_L theta_L",
          definition_f = "L = vech(chol(Sigma)) = f_L",
          theta_name = "theta_L",
          H_name = "H_L",
          f_name = "f_L",
          H = self$H_L,
          f = self$f_L,
          n_theta = self$length_theta_L
        )
      )
      cat("Parameter restrictions:\n\n")
      for (i in seq_along(restriction_blocks)) {
        block <- restriction_blocks[[i]]
        H <- block$H
        f <- block$f
        H_empty <- nrow(H) == 0 || ncol(H) == 0
        f_empty <- nrow(f) == 0 || ncol(f) == 0
        f_zero <- !f_empty && isTRUE(all(f == 0))
        show_f <- !f_empty && !f_zero
        if (H_empty && !show_f) {
          cat(sprintf("  No parameter %s\n", block$parameter_name))
        } else {
          definition <- block$definition_full
          if (H_empty) {
            definition <- block$definition_f
          } else if (!show_f) {
            definition <- block$definition_H
          }
          cat(sprintf("  %s\n", definition))
          cat(sprintf("  length(%s): %d\n", block$theta_name, block$n_theta))
          if (H_empty) {
            cat("\n")
          }
          if (!H_empty) {
            cat(sprintf("  %s (%d x %d) =\n", block$H_name, nrow(H), ncol(H)))
            print(H)
          }
          if (show_f) {
            cat(sprintf("  %s (%d x %d) =\n", block$f_name, nrow(f), ncol(f)))
            print(f)
          }
          cat("\n")
        }
        if (i < length(restriction_blocks)) {
          cat("\n")
        }
      }
    },

    # Warn about component adaptations.
    .warn_sync = function(source, adapted) {
      adapted <- setdiff(unique(adapted), source)
      if (self$suppress_warnings || length(adapted) == 0) {
        return(invisible(NULL))
      }
      adapted <- paste(sprintf("`%s`", adapted), collapse = ", ")
      warning(
        sprintf(
          "Setting `%s` adapted %s to keep model components conformable.",
          source,
          adapted
        ),
        call. = FALSE
      )
    },

    # Check random coefficient names against the formula.
    .check_re = function(covariate_types, ASC, re) {
      allowed_re <- unique(c(unlist(covariate_types), "ASC"))
      unknown_re <- setdiff(re, allowed_re)
      if (length(unknown_re) > 0) {
        stop(
          sprintf(
            "The following random coefficient variable(s) are unknown: %s.",
            paste(sprintf("`%s`", unknown_re), collapse = ", ")
          ),
          call. = FALSE
        )
      }
      if (!ASC && "ASC" %in% re) {
        stop(
          "Argument `ASC` cannot be random when ASCs are not included.",
          call. = FALSE
        )
      }
      invisible(TRUE)
    },

    # Return row labels for covariance-related vech vectors.
    .vech_labels = function(labels) {
      label_matrix <- outer(labels, labels, paste, sep = "_:_")
      label_matrix[lower.tri(label_matrix, diag = TRUE)]
    },

    # Set initial b restrictions.
    .setup_b = function(length_b) {
      list(
        H = diag(length_b),
        f = matrix(0, length_b, 1)
      )
    },

    # Set initial Omega restrictions.
    .setup_Omega = function(re_count, corr_re) {
      rows <- re_count * (re_count + 1) / 2
      if (corr_re) {
        H <- diag(rows)
      } else {
        index <- integer()
        if (re_count > 0) {
          index <- cumsum(c(1L, rev(seq_len(re_count - 1L)) + 1L))
        }
        H <- matrix(0, rows, length(index))
        if (length(index) > 0) {
          H[cbind(index, seq_along(index))] <- 1
        }
      }
      list(
        H = H,
        f = matrix(0, rows, 1)
      )
    },

    # Set initial Sigma restrictions.
    .setup_Sigma = function(alt, error_struc) {
      rows <- alt * (alt + 1) / 2
      diag_index <- integer()
      if (alt > 0) {
        diag_index <- cumsum(c(1L, rev(seq_len(alt - 1L)) + 1L))
      }
      fixed <- matrix(0, rows, 1)
      if (error_struc == "full") {
        index <- setdiff(seq_len(rows), seq_len(alt + 1L))
        fixed[alt + 1L, 1] <- 1
        H <- matrix(0, rows, length(index))
        if (length(index) > 0) {
          H[cbind(index, seq_along(index))] <- 1
        }
      }
      if (error_struc == "diag") {
        index <- diag_index[-1L]
        fixed[diag_index[1L], 1] <- 1
        H <- matrix(0, rows, length(index))
        if (length(index) > 0) {
          H[cbind(index, seq_along(index))] <- 1
        }
      }
      if (error_struc == "ident") {
        H <- matrix(0, rows, 1)
        H[diag_index, 1] <- 1
      }
      if (error_struc == "fixed") {
        H <- matrix(0, rows, 0)
        fixed[diag_index, 1] <- 1
      }
      list(H = H, f = fixed)
    },

    # Resize H matrix.
    .resize_H = function(x, rows, cols) {
      out <- matrix(0, rows, cols)
      if (rows > 0 && cols > 0 && nrow(x) > 0 && ncol(x) > 0) {
        row_index <- seq_len(min(rows, nrow(x)))
        col_index <- seq_len(min(cols, ncol(x)))
        out[row_index, col_index] <- x[row_index, col_index]
      }
      out
    },

    # Resize f vector.
    .resize_f = function(x, rows) {
      cols <- ncol(x)
      if (cols == 0) {
        cols <- 1
      }
      private$.resize_H(x, rows, cols)
    },

    # Check that row restrictions match the model specification.
    .assert_rows = function(value, rows, name) {
      if (nrow(value) != rows) {
        stop(
          sprintf(
            "Argument `%s` must have %d row(s) to match the model specification.",
            name, rows
          ),
          call. = FALSE
        )
      }
      invisible(TRUE)
    },

    # Synchronize dimensions implied by source-of-truth bindings.
    .sync_spec = function(
      source, sync_b = TRUE, sync_Sigma = TRUE, adapted = character()
    ) {
      old_H_b <- private$.H_b
      old_f_b <- private$.f_b
      old_H_O <- private$.H_O
      old_f_O <- private$.f_O
      old_H_L <- private$.H_L
      old_f_L <- private$.f_L
      if (sync_b) {
        b_count <- self$b_count
        rows_b <- sum(b_count)
        rows_Omega <- as.integer(
          b_count[2] * (b_count[2] + 1) / 2
        )
        private$.H_b <- private$.resize_H(
          private$.H_b,
          rows_b,
          ncol(private$.H_b)
        )
        private$.f_b <- private$.resize_f(private$.f_b, rows_b)
        if (nrow(old_H_O) != rows_Omega) {
          Hf_Omega <- private$.setup_Omega(b_count[2], corr_re = FALSE)
          private$.H_O <- Hf_Omega$H
          private$.f_O <- Hf_Omega$f
        } else {
          private$.H_O <- private$.resize_H(
            private$.H_O,
            rows_Omega,
            ncol(private$.H_O)
          )
          private$.f_O <- private$.resize_f(private$.f_O, rows_Omega)
        }
        private$.length_theta_b <- ncol(private$.H_b)
        private$.length_theta_O <- ncol(private$.H_O)
      }
      if (sync_Sigma) {
        alt <- self$alt
        rows <- alt * (alt + 1) / 2
        if (nrow(old_H_L) != rows) {
          Hf_Sigma <- private$.setup_Sigma(alt, error_struc = "full")
          private$.H_L <- Hf_Sigma$H
          private$.f_L <- Hf_Sigma$f
        } else {
          private$.H_L <- private$.resize_H(
            private$.H_L,
            rows,
            ncol(private$.H_L)
          )
          private$.f_L <- private$.resize_f(private$.f_L, rows)
        }
        private$.length_theta_L <- ncol(private$.H_L)
      }
      private$.warn_sync(
        source,
        c(
          if (!identical(private$.H_b, old_H_b)) "H_b",
          if (!identical(private$.f_b, old_f_b)) "f_b",
          if (!identical(private$.H_O, old_H_O)) "H_O",
          if (!identical(private$.f_O, old_f_O)) "f_O",
          adapted,
          if (!identical(private$.H_L, old_H_L)) "H_L",
          if (!identical(private$.f_L, old_f_L)) "f_L"
        )
      )
    }

  ),
  active = list(

    #' @field formula \[`formula`\]\cr
    #' Model formula. Writable; changing it adapts dependent restrictions.
    formula = function(value) {
      if (missing(value)) {
        private$.formula
      } else {
        covariate_types <- read_formula_rhs(value)
        ASC <- read_formula_asc(value)
        private$.check_re(
          covariate_types,
          ASC,
          private$.re
        )
        private$.formula <- value
        private$.covariate_types <- covariate_types
        private$.ASC <- ASC
        private$.sync_spec("formula", sync_b = TRUE, sync_Sigma = FALSE)
      }
    },

    #' @field re \[`character()`\]\cr
    #' Covariates with random coefficients. Writable; changing it adapts
    #' dependent restrictions.
    re = function(value) {
      if (missing(value)) {
        private$.re
      } else {
        checkmate::assert_character(value, any.missing = FALSE)
        value <- unique(value)
        private$.check_re(private$.covariate_types, private$.ASC, value)
        private$.re <- value
        private$.sync_spec("re", sync_b = TRUE, sync_Sigma = FALSE)
      }
    },

    #' @field alt_names \[`character()`\]\cr
    #' Alternative labels, at least 2. Writable; changing them adapts dependent
    #' restrictions.
    alt_names = function(value) {
      if (missing(value)) {
        private$.alt_names
      } else {
        old_alt <- private$.alt
        checkmate::assert_character(
          value, any.missing = FALSE, min.len = 2, unique = TRUE
        )
        value <- as.character(value)
        ref_alt_index <- match(private$.ref_alt, value)
        if (is.na(ref_alt_index)) {
          stop("Argument `ref_alt` is not in `alt_names`.", call. = FALSE)
        }
        ref_alt_index <- as.integer(ref_alt_index)
        private$.alt_names <- value
        private$.alt <- as.integer(length(private$.alt_names))
        private$.ref_alt_index <- ref_alt_index
        private$.sync_spec(
          "alt_names",
          sync_b = TRUE,
          sync_Sigma = TRUE,
          adapted = if (!identical(private$.alt, old_alt)) "alt"
        )
      }
    },

    #' @field ref_alt \[`character(1)`\]\cr
    #' Reference alternative. Writable; changing it adapts dependent
    #' restrictions.
    ref_alt = function(value) {
      if (missing(value)) {
        private$.ref_alt
      } else {
        old_var_names <- self$var_names
        old_b_labels <- unlist(old_var_names, use.names = FALSE)
        old_random_labels <- old_var_names$random
        old_O_labels <- private$.vech_labels(old_random_labels)
        old_H_b <- private$.H_b
        old_f_b <- private$.f_b
        old_H_O <- private$.H_O
        old_f_O <- private$.f_O
        checkmate::assert_string(value, min.chars = 1)
        ref_alt_index <- match(value, private$.alt_names)
        if (is.na(ref_alt_index)) {
          stop("Argument `ref_alt` is not in `alt_names`.", call. = FALSE)
        }
        private$.ref_alt <- value
        private$.ref_alt_index <- as.integer(ref_alt_index)
        new_var_names <- self$var_names
        new_b_labels <- unlist(new_var_names, use.names = FALSE)
        new_random_labels <- new_var_names$random
        new_O_labels <- private$.vech_labels(new_random_labels)
        cols <- ncol(old_H_b)
        H_b <- matrix(0, length(new_b_labels), cols)
        if (
          length(old_b_labels) > 0 &&
            length(new_b_labels) > 0 &&
            nrow(old_H_b) > 0 &&
            cols > 0
        ) {
          row_index <- match(new_b_labels, old_b_labels)
          keep <- !is.na(row_index)
          if (any(keep)) {
            H_b[keep, ] <- old_H_b[row_index[keep], , drop = FALSE]
          }
        }
        private$.H_b <- H_b
        if (ncol(private$.H_b) > 0) {
          keep <- colSums(abs(private$.H_b)) > 0
          private$.H_b <- private$.H_b[, keep, drop = FALSE]
        }
        new_free_b_labels <- new_b_labels[
          is.na(match(new_b_labels, old_b_labels))
        ]
        if (length(new_free_b_labels) > 0) {
          free <- matrix(0, length(new_b_labels), length(new_free_b_labels))
          row_index <- match(new_free_b_labels, new_b_labels)
          free[cbind(row_index, seq_along(new_free_b_labels))] <- 1
          private$.H_b <- cbind(private$.H_b, free)
        }
        cols <- ncol(old_f_b)
        if (cols == 0) {
          cols <- 1
        }
        f_b <- matrix(0, length(new_b_labels), cols)
        if (
          length(old_b_labels) > 0 &&
            length(new_b_labels) > 0 &&
            nrow(old_f_b) > 0 &&
            ncol(old_f_b) > 0
        ) {
          row_index <- match(new_b_labels, old_b_labels)
          keep <- !is.na(row_index)
          if (any(keep)) {
            f_b[keep, ] <- old_f_b[row_index[keep], , drop = FALSE]
          }
        }
        private$.f_b <- f_b
        old_diag_labels <- paste(
          old_random_labels, old_random_labels, sep = "_:_"
        )
        old_offdiag_labels <- setdiff(old_O_labels, old_diag_labels)
        old_diag_index <- match(old_diag_labels, old_O_labels)
        old_diag_index <- old_diag_index[!is.na(old_diag_index)]
        old_offdiag_index <- match(old_offdiag_labels, old_O_labels)
        old_offdiag_index <- old_offdiag_index[!is.na(old_offdiag_index)]
        diag_free <- length(old_diag_index) > 0 &&
          any(rowSums(abs(old_H_O[old_diag_index, , drop = FALSE])) > 0)
        offdiag_free <- length(old_offdiag_index) > 0 &&
          any(rowSums(abs(old_H_O[old_offdiag_index, , drop = FALSE])) > 0)
        new_diag_labels <- paste(
          new_random_labels, new_random_labels, sep = "_:_"
        )
        new_offdiag_labels <- setdiff(new_O_labels, new_diag_labels)
        free_O_labels <- c(
          if (diag_free) new_diag_labels,
          if (offdiag_free) new_offdiag_labels
        )
        cols <- ncol(old_H_O)
        H_O <- matrix(0, length(new_O_labels), cols)
        if (
          length(old_O_labels) > 0 &&
            length(new_O_labels) > 0 &&
            nrow(old_H_O) > 0 &&
            cols > 0
        ) {
          row_index <- match(new_O_labels, old_O_labels)
          keep <- !is.na(row_index)
          if (any(keep)) {
            H_O[keep, ] <- old_H_O[row_index[keep], , drop = FALSE]
          }
        }
        private$.H_O <- H_O
        if (ncol(private$.H_O) > 0) {
          keep <- colSums(abs(private$.H_O)) > 0
          private$.H_O <- private$.H_O[, keep, drop = FALSE]
        }
        new_free_O_labels <- new_O_labels[
          is.na(match(new_O_labels, old_O_labels)) &
            new_O_labels %in% free_O_labels
        ]
        if (length(new_free_O_labels) > 0) {
          free <- matrix(0, length(new_O_labels), length(new_free_O_labels))
          row_index <- match(new_free_O_labels, new_O_labels)
          free[cbind(row_index, seq_along(new_free_O_labels))] <- 1
          private$.H_O <- cbind(private$.H_O, free)
        }
        cols <- ncol(old_f_O)
        if (cols == 0) {
          cols <- 1
        }
        f_O <- matrix(0, length(new_O_labels), cols)
        if (
          length(old_O_labels) > 0 &&
            length(new_O_labels) > 0 &&
            nrow(old_f_O) > 0 &&
            ncol(old_f_O) > 0
        ) {
          row_index <- match(new_O_labels, old_O_labels)
          keep <- !is.na(row_index)
          if (any(keep)) {
            f_O[keep, ] <- old_f_O[row_index[keep], , drop = FALSE]
          }
        }
        private$.f_O <- f_O
        private$.length_theta_b <- ncol(private$.H_b)
        private$.length_theta_O <- ncol(private$.H_O)
        private$.warn_sync(
          "ref_alt",
          c(
            if (!identical(private$.H_b, old_H_b)) "H_b",
            if (!identical(private$.f_b, old_f_b)) "f_b",
            if (!identical(private$.H_O, old_H_O)) "H_O",
            if (!identical(private$.f_O, old_f_O)) "f_O"
          )
        )
      }
    },

    #' @field suppress_warnings \[`logical(1)`\]\cr
    #' Whether warnings are suppressed. Writable; defaults to `FALSE`.
    suppress_warnings = function(value) {
      if (missing(value)) {
        private$.suppress_warnings
      } else {
        checkmate::assert_flag(value)
        private$.suppress_warnings <- value
      }
    },

    #' @field f_b \[`matrix` | `numeric`\]\cr
    #' Fixed offset vector for `b`. Writable within the current row dimension.
    f_b = function(value) {
      if (missing(value)) {
        value <- private$.f_b
        rownames(value) <- unlist(self$var_names, use.names = FALSE)
        colnames(value) <- " "
        value
      } else {
        b_count <- self$b_count
        if (is.matrix(value)) {
          checkmate::assert_matrix(
            value, mode = "numeric", any.missing = FALSE
          )
          dimnames(value) <- NULL
        } else if (is.numeric(value)) {
          value <- matrix(as.numeric(value), ncol = 1)
        } else {
          stop(
            "Argument `f_b` must be a matrix or numeric vector.",
            call. = FALSE
          )
        }
        rows <- nrow(value)
        if (rows == 0) {
          value <- private$.resize_f(value, sum(b_count))
        }
        private$.assert_rows(value, sum(b_count), "f_b")
        if (ncol(value) != 1) {
          stop("Argument `f_b` must have one column.", call. = FALSE)
        }
        dimnames(value) <- NULL
        private$.f_b <- value
      }
    },

    #' @field H_b \[`matrix`\]\cr
    #' Transformation matrix for `b`. Writable within the current row dimension.
    H_b = function(value) {
      if (missing(value)) {
        value <- private$.H_b
        b_names <- unlist(self$var_names, use.names = FALSE)
        rownames(value) <- if (length(b_names) > 0) {
          paste0("b[", seq_along(b_names), "]:", b_names)
        } else {
          character()
        }
        colnames(value) <- if (ncol(value) > 0) {
          paste0("theta_", seq_len(ncol(value)))
        } else {
          character()
        }
        value
      } else {
        b_count <- self$b_count
        checkmate::assert_matrix(
          value, mode = "numeric", any.missing = FALSE
        )
        if (any(!is.finite(value))) {
          stop("Argument `H_b` must contain only finite values.", call. = FALSE)
        }
        private$.assert_rows(value, sum(b_count), "H_b")
        dimnames(value) <- NULL
        private$.H_b <- value
        private$.length_theta_b <- ncol(private$.H_b)
      }
    },

    #' @field f_O \[`matrix` | `numeric`\]\cr
    #' Fixed offset vector for `Omega`. Writable within the current row
    #' dimension.
    f_O = function(value) {
      if (missing(value)) {
        value <- private$.f_O
        var_names <- self$var_names$random
        rownames(value) <- private$.vech_labels(var_names)
        colnames(value) <- " "
        value
      } else {
        b_count <- self$b_count
        rows_current <- as.integer(b_count[2] * (b_count[2] + 1) / 2)
        if (is.matrix(value)) {
          checkmate::assert_matrix(
            value, mode = "numeric", any.missing = FALSE
          )
          dimnames(value) <- NULL
        } else if (is.numeric(value)) {
          value <- matrix(as.numeric(value), ncol = 1)
        } else {
          stop(
            "Argument `f_O` must be a matrix or numeric vector.",
            call. = FALSE
          )
        }
        rows <- nrow(value)
        if (rows == 0) {
          value <- private$.resize_f(value, rows_current)
        }
        private$.assert_rows(value, rows_current, "f_O")
        if (ncol(value) != 1) {
          stop("Argument `f_O` must have one column.", call. = FALSE)
        }
        dimnames(value) <- NULL
        private$.f_O <- value
      }
    },

    #' @field H_O \[`matrix`\]\cr
    #' Transformation matrix for `Omega`. Writable within the current row
    #' dimension.
    H_O = function(value) {
      if (missing(value)) {
        value <- private$.H_O
        var_names <- self$var_names$random
        rownames(value) <- private$.vech_labels(var_names)
        colnames(value) <- if (ncol(value) > 0) {
          paste0("theta_", self$length_theta_b + seq_len(ncol(value)))
        } else {
          character()
        }
        value
      } else {
        b_count <- self$b_count
        rows_current <- as.integer(b_count[2] * (b_count[2] + 1) / 2)
        checkmate::assert_matrix(
          value, mode = "numeric", any.missing = FALSE
        )
        if (any(!is.finite(value))) {
          stop("Argument `H_O` must contain only finite values.", call. = FALSE)
        }
        private$.assert_rows(value, rows_current, "H_O")
        dimnames(value) <- NULL
        private$.H_O <- value
        private$.length_theta_O <- ncol(private$.H_O)
      }
    },

    #' @field f_L \[`matrix` | `numeric`\]\cr
    #' Fixed offset vector for `Sigma`. Writable within the current row
    #' dimension.
    f_L = function(value) {
      if (missing(value)) {
        value <- private$.f_L
        alt_names <- self$alt_names
        rownames(value) <- private$.vech_labels(alt_names)
        colnames(value) <- " "
        value
      } else {
        if (is.matrix(value)) {
          checkmate::assert_matrix(
            value, mode = "numeric", any.missing = FALSE
          )
          dimnames(value) <- NULL
        } else if (is.numeric(value)) {
          value <- matrix(as.numeric(value), ncol = 1)
        } else {
          stop(
            "Argument `f_L` must be a matrix or numeric vector.",
            call. = FALSE
          )
        }
        alt <- self$alt
        rows_current <- alt * (alt + 1) / 2
        rows <- nrow(value)
        if (rows == 0) {
          value <- private$.resize_f(value, rows_current)
        }
        private$.assert_rows(value, rows_current, "f_L")
        if (ncol(value) != 1) {
          stop("Argument `f_L` must have one column.", call. = FALSE)
        }
        dimnames(value) <- NULL
        private$.f_L <- value
      }
    },

    #' @field H_L \[`matrix`\]\cr
    #' Transformation matrix for `Sigma`. Writable within the current row
    #' dimension.
    H_L = function(value) {
      if (missing(value)) {
        value <- private$.H_L
        alt_names <- self$alt_names
        rownames(value) <- private$.vech_labels(alt_names)
        colnames(value) <- if (ncol(value) > 0) {
          paste0(
            "theta_",
            self$length_theta_b + self$length_theta_O + seq_len(ncol(value))
          )
        } else {
          character()
        }
        value
      } else {
        checkmate::assert_matrix(
          value, mode = "numeric", any.missing = FALSE
        )
        if (any(!is.finite(value))) {
          stop("Argument `H_L` must contain only finite values.", call. = FALSE)
        }
        alt <- self$alt
        rows_current <- alt * (alt + 1) / 2
        private$.assert_rows(value, rows_current, "H_L")
        dimnames(value) <- NULL
        private$.H_L <- value
        private$.length_theta_L <- ncol(private$.H_L)
      }
    },

    #' @field alt (read only) \[`integer(1)`\]\cr
    #' Number of alternatives.
    alt = function(value) {
      if (!missing(value)) {
        stop("Field `alt` is read-only.", call. = FALSE)
      } else {
        private$.alt
      }
    },

    #' @field b_count (read only) \[`integer(2)`\]\cr
    #' Number of fixed and random coefficients.
    b_count = function(value) {
      if (!missing(value)) {
        stop("Field `b_count` is read-only.", call. = FALSE)
      }
      as.integer(lengths(self$var_names))
    },

    #' @field var_names (read only) \[`list()`\]\cr
    #' Names of fixed and random coefficients.
    var_names = function(value) {
      if (!missing(value)) {
        stop("Field `var_names` is read-only.", call. = FALSE)
      }
      alt_names <- private$.alt_names
      kept_alt_names <- alt_names[-private$.ref_alt_index]
      fixed <- character()
      random <- character()
      for (var in private$.covariate_types[[1]]) {
        if (var %in% private$.re) {
          random <- c(random, var)
        } else {
          fixed <- c(fixed, var)
        }
      }
      for (var in private$.covariate_types[[2]]) {
        names <- paste0(var, "_", kept_alt_names)
        if (var %in% private$.re) {
          random <- c(random, names)
        } else {
          fixed <- c(fixed, names)
        }
      }
      if (private$.ASC) {
        names <- paste0("ASC_", kept_alt_names)
        if ("ASC" %in% private$.re) {
          random <- c(random, names)
        } else {
          fixed <- c(fixed, names)
        }
      }
      for (var in private$.covariate_types[[3]]) {
        names <- paste0(var, "_", alt_names)
        if (var %in% private$.re) {
          random <- c(random, names)
        } else {
          fixed <- c(fixed, names)
        }
      }
      list(fixed = fixed, random = random)
    },

    #' @field length_theta_b (read only) \[`integer(1)`\]\cr
    #' Number of free parameters for `b`.
    length_theta_b = function(value) {
      if (missing(value)) {
        private$.length_theta_b
      } else {
        stop("Field `length_theta_b` is read-only.", call. = FALSE)
      }
    },

    #' @field length_theta_O (read only) \[`integer(1)`\]\cr
    #' Number of free parameters for `Omega`.
    length_theta_O = function(value) {
      if (missing(value)) {
        private$.length_theta_O
      } else {
        stop("Field `length_theta_O` is read-only.", call. = FALSE)
      }
    },

    #' @field length_theta_L (read only) \[`integer(1)`\]\cr
    #' Number of free parameters for `Sigma`.
    length_theta_L = function(value) {
      if (missing(value)) {
        private$.length_theta_L
      } else {
        stop("Field `length_theta_L` is read-only.", call. = FALSE)
      }
    },

    #' @field length_theta (read only) \[`integer(1)`\]\cr
    #' Total number of free model parameters.
    length_theta = function(value) {
      if (missing(value)) {
        private$.length_theta_b +
          private$.length_theta_O +
          private$.length_theta_L
      } else {
        stop("Field `length_theta` is read-only.", call. = FALSE)
      }
    },

    #' @field mod_list (read only) \[`list()`\]\cr
    #' Model components in the list format expected by [ll_cml]:
    #' - `length_theta_b`, `length_theta_O`, and `length_theta_L`, the number
    #'   of free parameters in the `b`, `Omega`, and `Sigma` parameter blocks.
    #' - `b_count`, the number of fixed and random coefficients.
    #' - `alt`, the number of choice alternatives.
    #' - `H_b` and `f_b`, the restriction matrix and offset vector for `b`.
    #' - `H_O` and `f_O`, the restriction matrix and offset vector for
    #'   `vech(chol(Omega))`.
    #' - `H_L` and `f_L`, the restriction matrix and offset vector for
    #'   `vech(chol(Sigma))`.
    mod_list = function(value) {
      if (!missing(value)) {
        stop("Field `mod_list` is read-only.", call. = FALSE)
      }
      list(
        length_theta_b = self$length_theta_b,
        length_theta_O = self$length_theta_O,
        length_theta_L = self$length_theta_L,
        b_count = self$b_count,
        alt = self$alt,
        H_b = self$H_b,
        f_b = as.numeric(self$f_b),
        H_O = self$H_O,
        f_O = as.numeric(self$f_O),
        H_L = self$H_L,
        f_L = as.numeric(self$f_L)
      )
    }
  )
)
