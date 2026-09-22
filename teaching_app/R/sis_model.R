## Simple population-level SIS model for trachoma teaching purposes.
## No age structure / immunity classes here (see trachoma_first_principles.R
## for the full model) - this is deliberately the minimal SIS system used to
## introduce the concept in the teaching app, mirroring WODIN's "basic" apps.
##
## State:  S (susceptible), I (infected), with S + I = N held constant.
## Params: beta  - effective contact/transmission rate (yr^-1)
##         gamma - recovery (spontaneous clearance) rate (yr^-1)
##         N     - population size
##         mda_coverage - fraction of the population successfully treated
##                        at each MDA round (moves I -> S)
##         mda_interval - years between MDA rounds (Inf / NA = no MDA)
##         mda_start    - year of first MDA round

library(deSolve)

sis_derivs <- function(t, state, parms) {
  S <- state[["S"]]
  I <- state[["I"]]
  N <- parms$N
  foi <- parms$beta * I / N
  dS <- -foi * S + parms$gamma * I
  dI <-  foi * S - parms$gamma * I
  list(c(S = dS, I = dI))
}

mda_event <- function(t, state, parms) {
  I_treated <- state[["I"]] * (1 - parms$mda_coverage)
  c(S = parms$N - I_treated, I = I_treated)
}

#' Run the SIS trachoma model
#'
#' @param beta transmission rate (yr^-1)
#' @param gamma recovery rate (yr^-1)
#' @param N population size
#' @param I0_prop initial proportion infected
#' @param t_end simulation horizon (years)
#' @param dt output time step (years)
#' @param mda_coverage coverage per MDA round (0 disables MDA)
#' @param mda_interval years between MDA rounds
#' @param mda_start year of the first MDA round
#' @return data.frame with columns time, S, I, prevalence
run_sis_model <- function(beta = 4,
                           gamma = 2,
                           N = 1000,
                           I0_prop = 0.3,
                           t_end = 20,
                           dt = 0.05,
                           mda_coverage = 0,
                           mda_interval = 1,
                           mda_start = 2) {
  parms <- list(beta = beta, gamma = gamma, N = N,
                mda_coverage = mda_coverage)
  state0 <- c(S = N * (1 - I0_prop), I = N * I0_prop)
  times <- seq(0, t_end, by = dt)

  use_mda <- mda_coverage > 0 && is.finite(mda_interval) && mda_interval > 0
  mda_times <- if (use_mda) seq(mda_start, t_end, by = mda_interval) else numeric(0)
  mda_times <- mda_times[mda_times <= t_end]

  out <- if (length(mda_times) > 0) {
    ode(y = state0, times = times, func = sis_derivs, parms = parms,
        events = list(func = mda_event, time = mda_times))
  } else {
    ode(y = state0, times = times, func = sis_derivs, parms = parms)
  }

  df <- as.data.frame(out)
  df$prevalence <- df$I / N
  attr(df, "mda_times") <- mda_times
  attr(df, "R0") <- beta / gamma
  attr(df, "equilibrium_prevalence") <- max(0, 1 - gamma / beta)
  df
}
