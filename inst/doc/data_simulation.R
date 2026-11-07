## ----global-options, include = FALSE------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE
)

## ----setup, echo = FALSE------------------------------------------------------
library(Rprobit)

## ----mode-canada-data---------------------------------------------------------
library(mlogit)
data("ModeCanada", package = "mlogit")
head(ModeCanada)

## ----create-mode-data---------------------------------------------------------
d_mode <- choice_data$new(
  df = ModeCanada,
  format = "long",
  column_decider_id = "case",
  column_alternative = "alt"
)

print(d_mode)

## ----inspect-mode-data-frame--------------------------------------------------
head(d_mode$df, 8)

## ----inspect-mode-data--------------------------------------------------------
d_mode$choice_counts

## ----inspect-mode-frequencies-------------------------------------------------
d_mode$choice_frequencies

## ----train-data---------------------------------------------------------------
data("Train", package = "mlogit")

head(Train)

## ----create-train-data--------------------------------------------------------
d_train <- choice_data$new(
  df = Train,
  column_decider_id = "id",
  column_occasion_id = "choiceid"
)
print(d_train)

## ----inspect-train-data-------------------------------------------------------
head(d_train$df, 6)

## ----build-design-------------------------------------------------------------
m_mode <- mod$new(
  formula = choice ~ cost + ivt | income | 0,
  alt = length(d_mode$alt_names),
  alt_names = d_mode$alt_names,
  ref_alt = d_mode$alt_names[1]
)

design_mode <- d_mode$build_design(m_mode)
design_mode[[1]]$X[[1]]

## ----build-second-design------------------------------------------------------
m_mode2 <- mod$new(
  formula = choice ~ cost | income | ivt,
  alt = length(d_mode$alt_names),
  alt_names = d_mode$alt_names,
  ref_alt = d_mode$alt_names[1]
)

design_mode <- d_mode$build_design(m_mode2)
design_mode[[1]]$X[[1]]

## ----difference-design--------------------------------------------------------
design_mode <- d_mode$build_design(m_mode, difference = TRUE)
design_mode[[1]]$X[[1]]

## ----difference-second-design-------------------------------------------------
design_mode <- d_mode$build_design(m_mode2, difference = TRUE)
design_mode[[1]]$X[[1]]

## ----simulate-data------------------------------------------------------------
set.seed(1)
m_sim <- mod$new(
  formula = choice ~ quality + price | income | time,
  re = "quality",
  alt = 2,
  alt_names = c("bus", "car"),
  ref_alt = "bus"
)

d_sim <- choice_data$new()
theta <- rnorm(m_sim$length_theta)
theta[1L] <- -1 # negative price coefficient

simulated <- d_sim$simulate(
  mod = m_sim,
  theta = theta,
  N = 3,
  Tp = c(1, 2, 3),
  sd = 2 # standard deviation of centered normal distribution for covariates
)

head(simulated, 6)

## ----inspect-simulated-data---------------------------------------------------
head(d_sim$df, 6)

## ----inspect-simulated-counts-------------------------------------------------
d_sim$choice_counts

## ----redraw-choices-----------------------------------------------------------
set.seed(1)
modified <- simulated
modified$price_bus <- modified$price_bus - 5

redrawn <- d_sim$redraw_choices(
  df = modified,
  mod = m_sim,
  column_decider_id = "id"
)

head(redrawn)

## ----inspect-redrawn-counts---------------------------------------------------
d_sim$choice_counts

