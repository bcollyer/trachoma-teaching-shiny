## Population-level "SIS ladder" model for trachoma teaching purposes.
##
## No age structure (single well-mixed population), but uses a configurable
## number of infection-history "rungs" (n_rungs, default 5) with the same
## continuous nu0/nu1/nuExponent (recovery) and lExponent (infectivity
## decay) formulation as the full 300-history model in
## trachoma_first_principles_300_history.R - just without the age
## dimension. This replaces the earlier hand-picked 2-rung (gamma1/gamma2,
## sigma/rho) version with the same mechanism the full model actually uses.
##
## State:  S_1..S_n (susceptible, by number of past infections),
##         I_1..I_n (currently infected, by number of past infections),
##         TS, TT (irreversible disease burden, driven by repeat infection)
##
## A person's j-th infection recovers into S_{j+1} (one rung up the immunity
## ladder); recovery from the *last* rung n loops back to S_n (capped, same
## as the full model). There is no separate susceptibility reduction by
## rung - immunity here is entirely "encoded" via faster recovery (higher
## nu_j) and lower infectiousness (lower l_j) at higher rungs, matching the
## full model exactly.
##
## Params: beta       - transmission intensity scaling the force of infection
##         nu0        - recovery rate at the first infection (rung 1)
##         nu1        - asymptotic recovery rate at high rungs
##         nuExponent - speed of transition from nu0 toward nu1 across rungs
##         lExponent  - decline in relative infectiousness with rung
##         n_rungs    - number of infection-history rungs (>= 1)
##         r_TS       - rate of scarring accumulation, driven by repeat infection (rungs >= 2)
##         r_TT       - rate of progression TS -> TT
##         alpha      - rate of TT resolution, e.g. trichiasis surgery (not in
##                      the full model - added here as a teaching control)
##         mu         - background birth/death rate = 1 / life expectancy
##         mda_coverage/mda_interval/mda_start - MDA cures I_j -> S_j at the
##           *same* rung, for every rung, matching make_mda_events() in the
##           full model.

library(deSolve)

ladder_derivs <- function(t, state, parms) {
  n <- parms$n_rungs
  S <- state[1:n]
  I <- state[(n + 1):(2 * n)]
  TS <- state[[2 * n + 1]]
  TT <- state[[2 * n + 2]]
  N <- parms$N
  mu <- parms$mu

  j_idx <- 1:n
  inf_load <- exp(-parms$lExponent * (j_idx - 1))
  nu <- (parms$nu0 - parms$nu1) * exp(-parms$nuExponent * (j_idx - 1)) + parms$nu1

  lambda <- parms$beta * sum(inf_load * I) / N
  incidence <- S * lambda
  recovery <- I * nu

  dS <- -incidence
  if (n > 1) dS[2:n] <- dS[2:n] + recovery[1:(n - 1)]
  dS[n] <- dS[n] + recovery[n]          # recovery from the top rung is capped there
  dS[1] <- dS[1] + mu * N
  dS <- dS - mu * S

  dI <- incidence - recovery - mu * I

  repeat_I <- if (n > 1) sum(I[2:n]) else 0
  new_TS <- parms$r_TS * (repeat_I / N) * max(N - TS, 0)
  new_TT <- parms$r_TT * (TS / N) * max(N - TT, 0)
  dTS <- new_TS - mu * TS
  dTT <- new_TT - (parms$alpha + mu) * TT

  list(c(dS, dI, dTS, dTT))
}

mda_event <- function(t, state, parms) {
  n <- parms$n_rungs
  cvr <- parms$mda_coverage
  S <- state[1:n]
  I <- state[(n + 1):(2 * n)]
  cured <- I * cvr
  c(S + cured, I - cured, state[[2 * n + 1]], state[[2 * n + 2]])
}

#' Run the no-age SIS ladder trachoma model
#'
#' @param beta transmission intensity (yr^-1)
#' @param nu0 recovery rate at rung 1 (yr^-1)
#' @param nu1 asymptotic recovery rate at high rungs (yr^-1)
#' @param nuExponent speed of nu0 -> nu1 transition across rungs
#' @param lExponent decline in relative infectiousness across rungs
#' @param n_rungs number of infection-history rungs
#' @param r_TS rate of scarring accumulation from repeat infection (yr^-1)
#' @param r_TT rate of progression TS -> TT (yr^-1)
#' @param alpha rate of TT resolution (e.g. surgery) (yr^-1)
#' @param life_expectancy years (sets mu = 1/life_expectancy)
#' @param N population size
#' @param I0_prop initial proportion of N infected (all in I_1, fresh introduction)
#' @param t_end simulation horizon (years)
#' @param dt output time step (years)
#' @param mda_coverage coverage per MDA round (0 disables MDA)
#' @param mda_interval years between MDA rounds
#' @param mda_start year of the first MDA round
#' @return data.frame with time, S_1..S_n, I_1..I_n, TS, TT,
#'   prevalence_infection, prevalence_TT
run_ladder_model <- function(beta = 3,
                              nu0 = 0.885594,
                              nu1 = 4.208717,
                              nuExponent = 0.751195,
                              lExponent = 0.05277,
                              n_rungs = 5,
                              r_TS = 0.10,
                              r_TT = 0.02,
                              alpha = 0.04,
                              life_expectancy = 55,
                              N = 1000,
                              I0_prop = 0.05,
                              t_end = 80,
                              dt = 0.1,
                              mda_coverage = 0,
                              mda_interval = 1,
                              mda_start = 20) {
  n <- as.integer(n_rungs)
  mu <- 1 / life_expectancy
  parms <- list(beta = beta, nu0 = nu0, nu1 = nu1, nuExponent = nuExponent,
                lExponent = lExponent, n_rungs = n, r_TS = r_TS, r_TT = r_TT,
                alpha = alpha, mu = mu, N = N, mda_coverage = mda_coverage)

  S0 <- rep(0, n); S0[1] <- N * (1 - I0_prop)
  I0 <- rep(0, n); I0[1] <- N * I0_prop
  state0 <- setNames(c(S0, I0, 0, 0),
                      c(sprintf("S_%d", 1:n), sprintf("I_%d", 1:n), "TS", "TT"))
  times <- seq(0, t_end, by = dt)

  use_mda <- mda_coverage > 0 && is.finite(mda_interval) && mda_interval > 0
  mda_times <- if (use_mda) seq(mda_start, t_end, by = mda_interval) else numeric(0)
  mda_times <- mda_times[mda_times <= t_end]

  out <- if (length(mda_times) > 0) {
    ode(y = state0, times = times, func = ladder_derivs, parms = parms,
        events = list(func = mda_event, time = mda_times))
  } else {
    ode(y = state0, times = times, func = ladder_derivs, parms = parms)
  }

  df <- as.data.frame(out)
  s_cols <- sprintf("S_%d", 1:n)
  i_cols <- sprintf("I_%d", 1:n)
  df$S_total <- rowSums(df[, s_cols, drop = FALSE])
  df$I_total <- rowSums(df[, i_cols, drop = FALSE])
  df$prevalence_infection <- df$I_total / N
  df$prevalence_TT <- df$TT / N

  attr(df, "mda_times") <- mda_times
  attr(df, "n_rungs") <- n
  attr(df, "R0_invasion") <- beta / (nu0 + mu)
  df
}
