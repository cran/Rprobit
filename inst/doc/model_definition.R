## ----global-options, include = FALSE------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE
)

## ----setup, echo = FALSE------------------------------------------------------
library(Rprobit)

## ----basic-formula, eval = FALSE----------------------------------------------
# formula <- choice ~ V1 | V2 | V3

## ----car-formula, eval = FALSE------------------------------------------------
# formula <- choice ~  PS + Main | HI + HS + loc + CF + 0 | range

## ----random-effects, eval = FALSE---------------------------------------------
# re <- c("V1", "V3")

## ----create-basic-model, eval = TRUE------------------------------------------
formula <- choice ~ price | 1 | time
re <- "ASC"
m_basic <- mod$new(formula, re, alt = 3)

## ----print-basic-model, eval = TRUE-------------------------------------------
print(m_basic)

## ----configure-model, eval = TRUE---------------------------------------------
alt_names <- c("train", "bus", "car")

m <- mod$new(
  formula = formula, re = re, alt = 3,
  alt_names = alt_names,            # alternative labels
  ref_alt = "train",                # reference alternative (for type-2 effects)
  corr_re = FALSE,                  # structure for Omega
  error_struc = "diag"              # structure for Sigma
)

print(m)

## ----add-random-coefficient, eval = TRUE--------------------------------------
m$suppress_warnings <- TRUE
m$re <- c("ASC", "price")
m

## ----reset-random-coefficient, eval = TRUE------------------------------------
m$suppress_warnings <- FALSE
m$re <- c("ASC")
m

## ----restrict-coefficient, eval = TRUE----------------------------------------
m$H_b <- m$H_b[, -1]
m$f_b[1] <- 1
m

## ----build-parameter----------------------------------------------------------
set.seed(1)
theta <- runif(m$length_theta)
(parameter <- m$build_parameter(theta))

## ----build-theta--------------------------------------------------------------
theta2 <- m$build_theta(
  b = parameter$b,
  Omega = parameter$Omega,
  Sigma = parameter$Sigma
)
all.equal(theta, theta2)

## ----check-identifiable-model, eval = TRUE------------------------------------
m$H_b <- diag(6)
m$check_identifiability()

## ----check-unidentifiable-model, eval = TRUE, warning = TRUE------------------
m2 <- m$clone()
m2$H_L <- diag(6)[, -c(1, 2, 3)]
m2$f_L <- matrix(0, ncol = 1, nrow = 6)
m2$check_identifiability()

