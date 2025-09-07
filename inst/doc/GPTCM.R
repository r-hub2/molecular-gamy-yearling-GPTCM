## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
# rm(list = ls())
# 
# # simulate data
# set.seed(123)
# n <- 200 # subjects
# p <- 10 # variable selection predictors
# L <- 3 # cell types
# library(GPTCM)
# set.seed(123)
# dat <- simData(n, p, L)
# 
# # KM curve
# library(survival)
# library(survminer)
# fit.km <- survival::survfit(Surv(time, event) ~ 1, data = dat$survObj)
# ggsurv <- survminer::ggsurvplot(fit.km,
#   conf.int = TRUE,
#   xlab = "Follow-up time (year)",
#   ylab = "Survival probability (%)",
#   legend = "none",
#   risk.table = TRUE,
#   cumevents = TRUE,
#   palette = "jco",
#   risk.table.title = "Number of patients at risk",
#   tables.height = 0.1,
#   tables.theme = theme_cleantable(),
#   tables.y.text = FALSE,
#   ggtheme = theme_light()
# )
# ggsurv$plot <- ggsurv$plot +
#   theme(
#     axis.text = element_text(size = 15),
#     axis.title = element_text(size = 15, face = "bold")
#   )
# ggsurv

## -----------------------------------------------------------------------------
# ## run Bayesian GPTCM
# set.seed(123)
# fit <- GPTCM(dat, nIter = 1100, burnin = 100)
# 
# # draw time-dependent Brier scores
# plotBrier(dat,
#   datMCMC = fit,
#   time.star = 3,
#   xlab = "Evalutation time points",
#   ylab = "Prediction error"
# )

## -----------------------------------------------------------------------------
# # show cel-type-specific effects
# plotCoeff(dat, datMCMC = fit, estimator = "beta", bandwidth = 0.02)
# # show BVS
# plotCoeff(dat, datMCMC = fit, estimator = "gamma")

## -----------------------------------------------------------------------------
# # show cel-type-specific effects
# plotCoeff(dat, datMCMC = fit, estimator = "zeta", bandwidth = 0.01)
# # show BVS
# plotCoeff(dat, datMCMC = fit, estimator = "eta")

