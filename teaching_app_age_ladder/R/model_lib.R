## Full model wrapper: adds the 18-group age structure from
## trachoma_first_principles_300_history.R on top of the infection-
## history "ladder" mechanics (nu0/nu1/nuExponent/lExponent, TS/TT) used in
## ../../teaching_app_ladder. Rather than reimplementing, this sources the
## real, working model file directly and just wraps it for Shiny use.
##
## The model file is vendored (copied) into this app's own R/ folder,
## rather than sourced via "../trachoma_first_principles_300_history.R",
## so this app is self-contained for shinyapps.io deployment (rsconnect
## only bundles files inside the app directory) - keep both copies in
## sync if the original in the SIR repo root changes.
##
## Reuses from trachoma_first_principles_300_history.R (unchanged):
##   make_parms_300(), run_model_300(), calc_R0_ngm(), age_widths, age_mids
## The demo/plotting code at the bottom of that file is guarded by
## `if (sys.nframe() == 0)`, so sourcing it here does not run it.

source(file.path("R", "trachoma_first_principles_300_history.R"))

age_band_labels <- c(
  "0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8", "8-9", "9-10",
  "10-15", "15-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-90"
)

#' Run the age-structured infection-history ladder model
#'
#' @param transRate transmission intensity (beta) (yr^-1)
#' @param epsilon age-assortative mixing (0 = random, 1 = within-age only)
#' @param nu0,nu1,nuExponent recovery-rate-by-rung parameters (see teaching_app_ladder)
#' @param lExponent infectiousness-decay-by-rung parameter
#' @param alpha birth rate (scales total population size)
#' @param k_TS,k_TT scarring / trichiasis accumulation rates
#' @param m_hist number of infection-history rungs
#' @param years simulation horizon (years)
#' @param dt_out output time step (years)
#' @param eq_window_years number of final years averaged to define "equilibrium"
#' @param mda_coverage,mda_start,mda_interval,mda_rounds,drug_efficacy MDA
#'   schedule (coverage 0 disables MDA; mda_rounds caps the number of
#'   rounds delivered, starting at mda_start, regardless of years)
#' @return list(age_df, ts_df, R0, mda_times)
run_age_ladder_model <- function(transRate = 3,
                                  epsilon = 0.5,
                                  nu0 = 0.885594,
                                  nu1 = 4.208717,
                                  nuExponent = 0.751195,
                                  lExponent = 0.05277,
                                  alpha = 0.04,
                                  k_TS = 0.08,
                                  k_TT = 0.02,
                                  m_hist = 5,
                                  years = 100,
                                  dt_out = 1,
                                  eq_window_years = 20,
                                  mda_coverage = 0,
                                  mda_start = 40,
                                  mda_interval = 1,
                                  mda_rounds = 10,
                                  drug_efficacy = 1.0) {
  parms <- make_parms_300(
    transRate = transRate, epsilon = epsilon, nu0 = nu0, nu1 = nu1,
    nuExponent = nuExponent, lExponent = lExponent, alpha = alpha,
    k_TS = k_TS, k_TT = k_TT
  )
  parms$m_hist <- m_hist

  mda_schedule <- NULL
  mda_times <- numeric(0)
  if (mda_coverage > 0 && is.finite(mda_interval) && mda_interval > 0 && mda_rounds > 0) {
    mda_times <- mda_start + (seq_len(mda_rounds) - 1) * mda_interval
    mda_times <- mda_times[mda_times <= years]
    if (length(mda_times) > 0) {
      mda_schedule <- data.frame(time = mda_times, coverage = mda_coverage)
    }
  }

  eq <- run_model_300(
    parms, years = years, dt_out = dt_out,
    eq_window_years = min(eq_window_years, years),
    mda_schedule = mda_schedule, drug_efficacy = drug_efficacy
  )

  ## Equilibrium (tail-window-averaged) age-aggregated prevalence
  N_age <- rowSums(eq$S + eq$I)
  okN <- N_age > 0
  prev_inf <- ifelse(okN, rowSums(eq$I) / N_age, NA_real_)
  prev_TS  <- ifelse(okN, eq$TS / N_age, NA_real_)
  prev_TT  <- ifelse(okN, eq$TT / N_age, NA_real_)

  age_df <- data.frame(
    age_mid = age_mids,
    age_label = factor(age_band_labels, levels = age_band_labels),
    prevalence_infection = prev_inf,
    prevalence_TS = prev_TS,
    prevalence_TT = prev_TT
  )

  ## Full time series, aggregated across all ages and rungs
  sol_df <- as.data.frame(eq$sol)
  names(sol_df)[1] <- "time"
  S_cols <- grep("^S_a", names(sol_df))
  I_cols <- grep("^I_a", names(sol_df))
  ts_df <- data.frame(
    time = sol_df$time,
    prevalence_infection = rowSums(sol_df[, I_cols, drop = FALSE]) /
      rowSums(sol_df[, c(S_cols, I_cols), drop = FALSE])
  )

  r0 <- tryCatch(calc_R0_ngm(parms), error = function(e) NA_real_)

  list(age_df = age_df, ts_df = ts_df, R0 = r0, mda_times = mda_times)
}
