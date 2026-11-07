## ----global-options, include = FALSE------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE
)

## ----setup, echo = FALSE------------------------------------------------------
library(Rprobit)

## ----create-train-model-------------------------------------------------------
data("Train", package = "mlogit")

train_model <- setup(
  formula = choice ~ price + time + change | 0 | comfort,
  df = Train,
  re = c("time"),
  column_decider_id = "id",
  column_occasion_id = "choiceid"
)
train_model$mod$H_b <- train_model$mod$H_b[, -1]
train_model$mod$f_b[1] <- -0.001
train_model$mod$H_L <- matrix(c(0, 0, 1), ncol = 1)
train_model$mod$f_L <- c(0, 0, 0)

print(train_model)

## ----fit-probit-model---------------------------------------------------------
start_time <- Sys.time()
train_model$fit(cml_type = "probit")
end_time <- Sys.time()
elapsed_time <- end_time - start_time

## ----inspect-probit-time------------------------------------------------------
elapsed_time

## ----compare-fit-init---------------------------------------------------------
rbind(
  "init" = train_model$theta_init,
  "fit" = train_model$theta_fit
) |> round(2)

## ----inspect-probit-fit-------------------------------------------------------
summary(train_model)

## ----inspect-probit-coefficients----------------------------------------------
coef(train_model)

## ----inspect-probit-confidence------------------------------------------------
confint(train_model)

## ----inspect-probit-vcov------------------------------------------------------
vcov(train_model)

## ----fit-ap-model-------------------------------------------------------------
start_time <- Sys.time()
train_2 <- fit(train_model, theta_init = NULL, print.level = 0, cml_type = "AP")
end_time <- Sys.time()
elapsed_time2 <- end_time - start_time

## ----inspect-ap-time----------------------------------------------------------
elapsed_time2

## ----inspect-ap-fit-----------------------------------------------------------
summary(train_2)

## ----compare-fits-------------------------------------------------------------
c(logLik(train_model), logLik(train_2))

## ----compare-fit-errors-------------------------------------------------------
V_probit <- vcov(train_model)
V_ap <- vcov(train_2)

cbind(sqrt(diag(V_probit)), sqrt(diag(V_ap)))

## ----model-diagnostics--------------------------------------------------------
res <- residuals(train_model)
apply(res, 2, mean)

## ----plot-diagnostics, fig.align = "center", fig.width = 7, fig.height = 5----
probs <- predict(train_model)
plot(probs[, 1], res[, 1], xlab = "Choice probabilities", ylab = "Residuals")

## ----predict-new-data---------------------------------------------------------
prediction_df <- Train[1:6, ]
prediction_df[, "price_A"] <- 3000
prediction_df[, "time_A"] <- 120
prediction_df[, "choice"] <- NA
prediction_data <- choice_data$new(
  df = prediction_df,
  column_decider_id = "id"
)
new_probabilities <- predict(
  train_model,
  newdata = prediction_data,
  type = "probabilities"
)
new_probabilities

