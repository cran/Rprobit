## ----global-options, include = FALSE------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE
)

## ----setup, echo = FALSE------------------------------------------------------
library(Rprobit)

## ----setup-train-model--------------------------------------------------------
data("Train", package = "mlogit")

train_model <- setup(
  formula = choice ~ price + time + change | 0 | comfort,
  df = Train,
  column_decider_id = "id",
  column_occasion_id = "choiceid"
)
train_model$mod$H_b <- train_model$mod$H_b[, -1]
train_model$mod$f_b[1] <- -1
train_model$mod$H_L <- matrix(c(0, 0, 1), ncol = 1)
train_model$mod$f_L <- c(0, 0, 0)

print(train_model)

## ----inspect-train-model------------------------------------------------------
train_model$choice_data

## ----inspect-train-specification----------------------------------------------
train_model$mod

## ----fit-train-model----------------------------------------------------------
train_model$fit(print.level = 0, cml_type = "probit")

## ----inspect-fit-theta--------------------------------------------------------
train_model$theta_fit

## ----inspect-fit-ll-----------------------------------------------------------
train_model$ll_fit

## ----inspect-fit-vcov---------------------------------------------------------
train_model$theta_vcov

## ----inspect-fit-hessian------------------------------------------------------
train_model$hessian_fit

## ----inspect-fit-score-covariance---------------------------------------------
train_model$Vhat_fit

## ----summarize-fit------------------------------------------------------------
summary(train_model)

## ----refit-train-model--------------------------------------------------------
train_model$mod$H_b <- diag(5)
train_model$mod$f_b[1] <- 0
train_model$mod$H_L <- matrix(0, 3, 0)
train_model$mod$f_L[3] <- 1

train_model$fit(print.level = 0, cml_type = "probit")

## ----summarize-refit----------------------------------------------------------
summary(train_model)

## ----setup-mode-model---------------------------------------------------------
data("ModeCanada", package = "mlogit")

mode_model <- setup(
  formula = choice ~ cost + ivt | income | 0,
  df = ModeCanada,
  format = "long",
  column_decider_id = "case",
  column_alternative = "alt",
  ref_alt = "train"
)

## ----fit-mode-model-----------------------------------------------------------
mode_model$fit(cml_type = "FP")

## ----summarize-mode-model-----------------------------------------------------
summary(mode_model)

## ----simulate-choice-data-----------------------------------------------------
set.seed(1)
sim_model <- setup(
  formula = choice ~ quality + price | income | time,
  re = "quality",
  alt_names = c("bus", "car", "train"),
  ref_alt = "bus",
  corr_re = FALSE,
  error_struc = "diag",
  N = 300,
  Tp = rep(3, 300),
  sd = 2
)

print(sim_model)

## ----inspect-simulated-data---------------------------------------------------
head(sim_model$choice_data$df)

## ----inspect-simulated-parameters---------------------------------------------
sim_model$choice_data$theta_dgp

## ----inspect-simulated-counts-------------------------------------------------
sim_model$choice_data$choice_counts

## ----fit-simulated-model------------------------------------------------------
sim_model$fit(
  theta_init = NULL,
  print.level = 0,
  cml_type = "probit"
)

## ----summarize-simulated-fit--------------------------------------------------
summary(sim_model)

## ----compare-simulated-fit----------------------------------------------------
rbind(
  "fit" = sim_model$theta_fit,
  "dgp" = sim_model$choice_data$theta_dgp
) |> round(2)

