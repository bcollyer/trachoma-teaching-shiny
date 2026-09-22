#!/usr/bin/env Rscript

if (!requireNamespace("deSolve", quietly = TRUE)) {
  stop("Package 'deSolve' is required. Install with install.packages('deSolve').")
}

library(deSolve)

state_name <- function(prefix, a, j = NULL) {
  if (is.null(j)) {
    sprintf("%s_a%02d", prefix, a)
  } else {
    sprintf("%s_a%02d_j%03d", prefix, a, j)
  }
}

make_initial_state <- function(n = 18, m = 100) {
  S <- matrix(0, nrow = n, ncol = m)
  I <- matrix(0, nrow = n, ncol = m)
  TS <- rep(0, n)
  TT <- rep(0, n)

  S[, 1] <- c(
    15.7468, 15.6310, 15.5161, 15.4020, 15.2888, 15.2559,
    15.2232, 15.1905, 15.1578, 15.1253, 75.3436, 74.7597,
    146.4487, 140.3511, 131.5182, 117.0465, 92.0111, 68.9837
  )

  S[1, 1] <- S[1, 1] - 1
  I[1, 1] <- 1

  c(
    setNames(
      as.vector(S),
      as.vector(outer(1:n, 1:m, Vectorize(function(a, j) state_name("S", a, j)))))
    ,
    setNames(
      as.vector(I),
      as.vector(outer(1:n, 1:m, Vectorize(function(a, j) state_name("I", a, j)))))
    ,
    setNames(TS, vapply(1:n, function(a) state_name("TS", a), character(1))),
    setNames(TT, vapply(1:n, function(a) state_name("TT", a), character(1)))
  )
}

trachoma_rhs_300 <- function(t, y, p) {
  n <- p$n_age
  m <- p$m_hist

  expected_len <- 2 * n * m + 2 * n # S + I + TS + TT
  if (length(y) != expected_len) {
    stop(sprintf(
      "State length mismatch at t=%.6f: got %d, expected %d (n=%d, m=%d).",
      t, length(y), expected_len, n, m
    ))
  }

  # Guard against non-finite state
  if (any(!is.finite(y))) y[!is.finite(y)] <- 0

  n_cells <- n * m
  S  <- matrix(y[1:n_cells], nrow = n, ncol = m, byrow = FALSE)
  I  <- matrix(y[(n_cells + 1):(2 * n_cells)], nrow = n, ncol = m, byrow = FALSE)
  TS <- y[(2 * n_cells + 1):(2 * n_cells + n)]
  TT <- y[(2 * n_cells + n + 1):(2 * n_cells + 2 * n)]

  # Keep states non-negative
  S  <- pmax(S, 0)
  I  <- pmax(I, 0)
  TS <- pmax(TS, 0)
  TT <- pmax(TT, 0)

  # Population by age
  N_age <- rowSums(S + I)
  N_tot <- sum(N_age)

  # Mixing matrix
  if (is.finite(N_tot) && N_tot > 0) {
    w_random <- matrix(rep(N_age / N_tot, each = n), nrow = n, ncol = n)
  } else {
    w_random <- matrix(0, nrow = n, ncol = n)
  }
  w <- (1 - p$epsilon) * w_random + p$epsilon * diag(n)

  # Infection-history-dependent infectivity
  j_idx <- 1:m
  inf_load <- exp(-p$lExponent * (j_idx - 1))
  weighted_I <- as.vector(I %*% inf_load)

  inf_frac <- rep(0, n)
  okN <- is.finite(N_age) & (N_age > 0)
  inf_frac[okN] <- weighted_I[okN] / N_age[okN]

  lambda <- as.vector(p$transRate * (w %*% inf_frac))
  lambda[!is.finite(lambda)] <- 0
  lambda <- pmax(lambda, 0)

  # Recovery by history stage
  nu <- (p$nu0 - p$nu1) * exp(-p$nuExponent * (j_idx - 1)) + p$nu1

  # Broadcast helpers
  lambda_mat <- matrix(lambda, nrow = n, ncol = m)
  nu_mat     <- matrix(nu, nrow = n, ncol = m, byrow = TRUE)
  mu_mat     <- matrix(rep(p$mu, m), nrow = n, ncol = m)
  aw_mat     <- matrix(rep(p$age_width, m), nrow = n, ncol = m)

  incidence <- S * lambda_mat
  recovery  <- I * nu_mat

  # Aging flows (vectorised)
  age_out_S <- S / aw_mat
  age_out_I <- I / aw_mat

  age_in_S <- matrix(0, nrow = n, ncol = m)
  age_in_I <- matrix(0, nrow = n, ncol = m)
  if (n > 1) {
    age_in_S[2:n, ] <- S[1:(n - 1), , drop = FALSE] / p$age_width[1:(n - 1)]
    age_in_I[2:n, ] <- I[1:(n - 1), , drop = FALSE] / p$age_width[1:(n - 1)]
  }

  # Base S/I dynamics
  dS <- age_in_S - age_out_S - incidence - mu_mat * S
  dI <- age_in_I - age_out_I + incidence - recovery - mu_mat * I

  # Births into S[1,1]
  dS[1, 1] <- dS[1, 1] + p$alpha * N_tot

  # History progression: I[a,j] recover to S[a,j+1], capped at j=m
  if (m > 1) {
    dS[, 2:m] <- dS[, 2:m] + recovery[, 1:(m - 1), drop = FALSE]
    dS[, m]   <- dS[, m] + recovery[, m]
  } else {
    dS[, 1] <- dS[, 1] + recovery[, 1]
  }

  # TS/TT burden dynamics
  repeat_I <- if (m > 1) rowSums(I[, 2:m, drop = FALSE]) else rep(0, n)
  repeat_prev <- rep(0, n)
  repeat_prev[okN] <- repeat_I[okN] / N_age[okN]

  ts_prev <- rep(0, n)
  ts_prev[okN] <- TS[okN] / N_age[okN]

  new_TS <- p$k_TS * repeat_prev * pmax(N_age - TS, 0)
  new_TT <- p$k_TT * ts_prev * pmax(N_age - TT, 0)

  age_out_TS <- TS / p$age_width
  age_out_TT <- TT / p$age_width

  age_in_TS <- rep(0, n)
  age_in_TT <- rep(0, n)
  if (n > 1) {
    age_in_TS[2:n] <- TS[1:(n - 1)] / p$age_width[1:(n - 1)]
    age_in_TT[2:n] <- TT[1:(n - 1)] / p$age_width[1:(n - 1)]
  }

  dTS <- age_in_TS - age_out_TS + new_TS - p$mu * TS
  dTT <- age_in_TT - age_out_TT + new_TT - p$mu * TT

  deriv <- c(as.vector(dS), as.vector(dI), dTS, dTT)
  deriv[!is.finite(deriv)] <- 0

  list(deriv)
}

normalise_mda_schedule <- function(mda_schedule) {
  if (is.null(mda_schedule)) return(NULL)

  if (is.data.frame(mda_schedule)) {
    if (!all(c("time", "coverage") %in% names(mda_schedule))) {
      stop("mda_schedule data.frame must contain columns named 'time' and 'coverage'.")
    }
    sch <- mda_schedule[, c("time", "coverage"), drop = FALSE]
  } else if (is.list(mda_schedule) && all(c("time", "coverage") %in% names(mda_schedule))) {
    sch <- data.frame(time = mda_schedule$time, coverage = mda_schedule$coverage)
  } else if (is.matrix(mda_schedule) && ncol(mda_schedule) == 2) {
    sch <- data.frame(time = mda_schedule[, 1], coverage = mda_schedule[, 2])
    names(sch) <- c("time", "coverage")
  } else {
    stop("mda_schedule must be NULL, a data.frame(time, coverage), a list(time=..., coverage=...), or a two-column matrix.")
  }

  if (nrow(sch) == 0) return(NULL)

  sch$time <- as.numeric(sch$time)
  sch$coverage <- as.numeric(sch$coverage)

  if (any(!is.finite(sch$time))) stop("mda_schedule$time contains non-finite values.")
  if (any(!is.finite(sch$coverage))) stop("mda_schedule$coverage contains non-finite values.")
  if (any(sch$coverage < 0 | sch$coverage > 1)) stop("mda_schedule$coverage must be in [0, 1].")

  sch <- sch[order(sch$time), , drop = FALSE]
  rownames(sch) <- NULL
  sch
}

make_mda_events <- function(mda_schedule, drug_efficacy, n, m) {
  sch <- normalise_mda_schedule(mda_schedule)
  if (is.null(sch)) return(NULL)

  eff <- as.numeric(drug_efficacy)
  if (!is.finite(eff)) stop("drug_efficacy must be finite.")
  if (eff < 0 || eff > 1) stop("drug_efficacy must be in [0, 1].")

  times <- sch$time
  coverages <- sch$coverage

  event_func <- function(t, y, parms) {
    idx <- which(abs(times - t) < 1e-10)
    if (length(idx) == 0) return(y)

    coverage_t <- sum(coverages[idx])
    eff_by_age <- rep(eff, n)
    if (n >= 1) eff_by_age[1] <- 0.25 * eff
    p_move_by_age <- pmax(pmin(coverage_t * eff_by_age, 1), 0)
    if (!any(p_move_by_age > 0)) return(y)

    n_cells <- n * m
    S <- matrix(y[1:n_cells], nrow = n, ncol = m, byrow = FALSE)
    I <- matrix(y[(n_cells + 1):(2 * n_cells)], nrow = n, ncol = m, byrow = FALSE)

    move_mat <- matrix(p_move_by_age, nrow = n, ncol = m)
    moved <- move_mat * I
    S <- S + moved
    I <- I - moved

    c(as.vector(S), as.vector(I), y[(2 * n_cells + 1):length(y)])
  }

  list(func = event_func, time = times)
}

run_trachoma_300 <- function(duration = 110, dt_out = 0.1, m_hist = 100,
                             mda_schedule = NULL, drug_efficacy = 1.0) {
  n <- 18

  params <- list(
      n_age = n,                 # Number of age groups.
      m_hist = m_hist,           # Number of infection-history stages.
      transRate = 1.0,           # Transmission intensity scaling force of infection.
      nu0 = 0.885594,            # Recovery rate at first infection-history stage.
      nu1 = 4.208717,            # Asymptotic recovery rate at high history stage.
      nuExponent = 0.751195,     # Speed of transition from nu0 toward nu1 across history.
      lExponent = 0.05277,       # Decline in infectivity with repeated infections.
      epsilon = 0.5,             # Assortative age mixing (0=random, 1=within-age only).
      alpha = 0.04,              # Birth rate into youngest, naive susceptible compartment.
      k_TS = 0.08,               # Rate of TS accumulation from repeat-infection burden.
      k_TT = 0.02,               # Rate of TT accumulation from TS burden.
    mu = c(0.09378, rep(0.00738, 4), rep(0.00215, 5), 0.00125,
        0.00176, 0.00292, 0.00458, 0.0069, 0.01266, 0.02645, 0.06349), # Age-specific mortality rates.
      age_width = c(rep(1.0, 10), rep(5.0, 2), rep(10.0, 5), 20.0) # Width of each age bin in years.
  )

  y0 <- make_initial_state(n = n, m = m_hist)
  times <- seq(0, duration, by = dt_out)
  mda_events <- make_mda_events(mda_schedule, drug_efficacy, n = n, m = m_hist)

  ode_args <- list(
    y = y0,
    times = times,
    func = trachoma_rhs_300,
    parms = params,
    method = "lsodes",
    rtol = 1e-6,
    atol = 1e-8,
    maxsteps = 2e5
  )

  if (!is.null(mda_events)) {
    ode_args$events <- mda_events
  }

  out <- do.call(deSolve::ode, ode_args)

  as.data.frame(out)
}

summarise_prevalence <- function(sim_df, n_age = 18, m_hist = 100) {
  out <- data.frame(time = sim_df$time)

  for (a in 1:n_age) {
    s_cols <- sprintf("S_a%02d_j%03d", a, 1:m_hist)
    i_cols <- sprintf("I_a%02d_j%03d", a, 1:m_hist)

    N_age <- rowSums(sim_df[, c(s_cols, i_cols), drop = FALSE])
    I_age <- rowSums(sim_df[, i_cols, drop = FALSE])

    out[[sprintf("inf_prev_a%02d", a)]] <- ifelse(N_age > 0, I_age / N_age, NA_real_)
    out[[sprintf("ts_prev_a%02d", a)]] <- ifelse(N_age > 0, sim_df[[sprintf("TS_a%02d", a)]] / N_age, NA_real_)
    out[[sprintf("tt_prev_a%02d", a)]] <- ifelse(N_age > 0, sim_df[[sprintf("TT_a%02d", a)]] / N_age, NA_real_)
  }

  out
}

if (sys.nframe() == 0) {
  sim <- run_trachoma_300(duration = 30, dt_out = 0.5)
  prev <- summarise_prevalence(sim)

  message("Simulation complete.")
  print(utils::head(prev[, 1:7]))
}

age_widths <- c(rep(1, 10), 5, 5, rep(10, 5), 20)
age_mids   <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5,
                7.5, 8.5, 9.5, 12, 17, 25, 35, 45, 55, 65, 80)
n_age <- length(age_widths)

make_parms_300 <- function(transRate = 1.8,
                           epsilon = 0.5,
                           nu0 = 0.885594,
                           nu1 = 4.208717,
                           nuExponent = 0.751195,
                           lExponent = 0.05277,
                           alpha = 0.04,
                           k_TS = 0.08,
                           k_TT = 0.02) {
  list(
    n_age = 18,                # Number of age groups.
    m_hist = 30,               # Number of infection-history stages.
    transRate = transRate,     # Transmission intensity scaling force of infection.
    nu0 = nu0,                 # Recovery rate at first infection-history stage.
    nu1 = nu1,                 # Asymptotic recovery rate at high history stage.
    nuExponent = nuExponent,   # Speed of transition from nu0 toward nu1 across history.
    lExponent = lExponent,     # Decline in infectivity with repeated infections.
    epsilon = epsilon,         # Assortative age mixing (0=random, 1=within-age only).
    alpha = alpha,             # Birth rate into youngest, naive susceptible compartment.
    k_TS = k_TS,               # Rate of TS accumulation from repeat-infection burden.
    k_TT = k_TT,               # Rate of TT accumulation from TS burden.
    mu = c(
      0.09378, rep(0.00738, 4), rep(0.00215, 5), 0.00125,
      0.00176, 0.00292, 0.00458, 0.00690, 0.01266, 0.02645, 0.06349
    ),                         # Age-specific mortality rates.
    age_width = age_widths     # Width of each age bin in years.
  )
}

run_model_300 <- function(parms, years = 500, dt_out = 0.25, eq_window_years = 20,
                          mda_schedule = NULL, drug_efficacy = 1.0) {
  # If mda_schedule is provided, events are applied at those exact model times.
  y0 <- make_initial_state(n = parms$n_age, m = parms$m_hist)
  mda_events <- make_mda_events(
    mda_schedule,
    drug_efficacy,
    n = parms$n_age,
    m = parms$m_hist
  )

  ode_args <- list(
    y = y0,
    times = seq(0, years, by = dt_out),
    func = trachoma_rhs_300,
    parms = parms,
    method = "lsoda",
    rtol = 1e-8,
    atol = 1e-10,
    hmax = 0.05,
    maxsteps = 1e6
  )

  if (!is.null(mda_events)) {
    ode_args$events <- mda_events
  }

  sol <- do.call(deSolve::ode, ode_args)

  sol_df <- as.data.frame(sol)
  tmax <- max(sol_df$time)
  eq_df <- sol_df[sol_df$time >= (tmax - eq_window_years), , drop = FALSE]
  ybar <- colMeans(eq_df[, -1, drop = FALSE])

  n <- parms$n_age
  m <- parms$m_hist
  n_cells <- n * m

  S  <- matrix(pmax(ybar[1:n_cells], 0), nrow = n, ncol = m, byrow = FALSE)
  I  <- matrix(pmax(ybar[(n_cells + 1):(2 * n_cells)], 0), nrow = n, ncol = m, byrow = FALSE)
  TS <- pmax(ybar[(2 * n_cells + 1):(2 * n_cells + n)], 0)
  TT <- pmax(ybar[(2 * n_cells + n + 1):(2 * n_cells + 2 * n)], 0)

  # Per-rung age vectors for direct access, e.g. S_rungs[["S_j003"]].
  S_rungs <- setNames(lapply(seq_len(m), function(j) S[, j]), sprintf("S_j%03d", seq_len(m)))
  I_rungs <- setNames(lapply(seq_len(m), function(j) I[, j]), sprintf("I_j%03d", seq_len(m)))

  list(
    S = S,
    I = I,
    S_rungs = S_rungs,
    I_rungs = I_rungs,
    S1 = S[, 1],
    I1 = I[, 1],
    S2 = rowSums(S[, 2:m, drop = FALSE]),
    I2 = rowSums(I[, 2:m, drop = FALSE]),
    TS = TS,
    TT = TT,
    sol = sol
  )
}

calc_R0_ngm <- function(parms, N_age_dfe = NULL) {
  n <- as.integer(parms$n_age)
  m <- as.integer(parms$m_hist)

  if (is.null(N_age_dfe)) {
    y0 <- make_initial_state(n = n, m = m)
    n_cells <- n * m
    S0 <- matrix(y0[1:n_cells], nrow = n, ncol = m, byrow = FALSE)
    I0 <- matrix(y0[(n_cells + 1):(2 * n_cells)], nrow = n, ncol = m, byrow = FALSE)
    N_age <- rowSums(S0 + I0)
  } else {
    N_age <- as.numeric(N_age_dfe)
  }

  if (length(N_age) != n) stop("N_age_dfe must have length n_age.")
  if (any(!is.finite(N_age)) || any(N_age <= 0)) stop("N_age_dfe must be finite and > 0 in all age groups.")

  N_tot <- sum(N_age)
  w_random <- matrix(rep(N_age / N_tot, each = n), nrow = n, ncol = n)
  w <- (1 - parms$epsilon) * w_random + parms$epsilon * diag(n)

  j_idx <- 1:m
  inf_load <- exp(-parms$lExponent * (j_idx - 1))
  nu <- (parms$nu0 - parms$nu1) * exp(-parms$nuExponent * (j_idx - 1)) + parms$nu1

  dim_I <- n * m
  F <- matrix(0, nrow = dim_I, ncol = dim_I)
  V <- matrix(0, nrow = dim_I, ncol = dim_I)

  idx <- function(a, j) as.integer((j - 1L) * n + a)

  # New infections only enter I[a,1] at the DFE.
  for (a in 1:n) {
    row <- idx(a, 1)
    for (b in 1:n) {
      for (j in 1:m) {
        col <- idx(b, j)
        F[row, col] <- parms$transRate * w[a, b] * inf_load[j]
      }
    }
  }

  # Multiply by susceptible mass in recipient age class and divide by source age size.
  for (a in 1:n) {
    row <- idx(a, 1)
    for (b in 1:n) {
      cols_b <- vapply(1:m, function(j) idx(b, j), integer(1))
      F[row, cols_b] <- F[row, cols_b] * (N_age[a] / N_age[b])
    }
  }

  # V captures non-infection transitions in the infected subsystem.
  for (j in 1:m) {
    for (a in 1:n) {
      k <- idx(a, j)
      V[k, k] <- 1 / parms$age_width[a] + nu[j] + parms$mu[a]
      if (a > 1) {
        V[k, idx(a - 1, j)] <- -1 / parms$age_width[a - 1]
      }
    }
  }

  K <- F %*% solve(V)
  eig <- eigen(K, only.values = TRUE)$values
  Re(eig[which.max(Re(eig))])
}

calc_Re_ngm <- function(parms, S_age, N_age = NULL) {
  n <- as.integer(parms$n_age)
  m <- as.integer(parms$m_hist)

  S_age <- as.numeric(S_age)
  if (is.null(N_age)) {
    N_age <- as.numeric(S_age)
  } else {
    N_age <- as.numeric(N_age)
  }

  if (length(S_age) != n) stop("S_age must have length n_age.")
  if (length(N_age) != n) stop("N_age must have length n_age.")
  if (any(!is.finite(S_age)) || any(S_age < 0)) stop("S_age must be finite and >= 0.")
  if (any(!is.finite(N_age)) || any(N_age <= 0)) stop("N_age must be finite and > 0.")
  if (any(S_age > N_age + 1e-12)) stop("S_age cannot exceed N_age in any age group.")

  N_tot <- sum(N_age)
  w_random <- matrix(rep(N_age / N_tot, each = n), nrow = n, ncol = n)
  w <- (1 - parms$epsilon) * w_random + parms$epsilon * diag(n)

  j_idx <- 1:m
  inf_load <- exp(-parms$lExponent * (j_idx - 1))
  nu <- (parms$nu0 - parms$nu1) * exp(-parms$nuExponent * (j_idx - 1)) + parms$nu1

  dim_I <- n * m
  F <- matrix(0, nrow = dim_I, ncol = dim_I)
  V <- matrix(0, nrow = dim_I, ncol = dim_I)

  idx <- function(a, j) as.integer((j - 1L) * n + a)

  for (a in 1:n) {
    row <- idx(a, 1)
    for (b in 1:n) {
      for (j in 1:m) {
        col <- idx(b, j)
        F[row, col] <- parms$transRate * w[a, b] * inf_load[j]
      }
    }
  }

  # At current state, susceptibility scales by S_age[a] rather than N_age[a].
  for (a in 1:n) {
    row <- idx(a, 1)
    for (b in 1:n) {
      cols_b <- vapply(1:m, function(j) idx(b, j), integer(1))
      F[row, cols_b] <- F[row, cols_b] * (S_age[a] / N_age[b])
    }
  }

  for (j in 1:m) {
    for (a in 1:n) {
      k <- idx(a, j)
      V[k, k] <- 1 / parms$age_width[a] + nu[j] + parms$mu[a]
      if (a > 1) {
        V[k, idx(a - 1, j)] <- -1 / parms$age_width[a - 1]
      }
    }
  }

  K <- F %*% solve(V)
  eig <- eigen(K, only.values = TRUE)$values
  Re(eig[which.max(Re(eig))])
}


# Demo/plotting only: skipped when this file is sourced as a library (e.g. by fitting scripts).
if (sys.nframe() == 0) {
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').")
}
library(ggplot2)

prev_data <- data.frame(
  age_mid = c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,
              13,18,25,35,45,55,65,80),
  pos     = c(3,7,4,4,4,5,5,5,5,5,14,6,11,10,1,2,2,1),
  n_test  = c(49,64,48,55,52,65,47,62,52,41,194,116,141,130,89,45,43,24)
)
prev_data$prev  <- prev_data$pos / prev_data$n_test
prev_data$lower <- qbeta(0.025, prev_data$pos + 0.5, prev_data$n_test - prev_data$pos + 0.5)
prev_data$upper <- qbeta(0.975, prev_data$pos + 0.5, prev_data$n_test - prev_data$pos + 0.5)

parms_fit <- make_parms_300(transRate = 4.8, epsilon = 0.5)
eq_fit <- run_model_300(parms_fit)

N_age_fit <- rowSums(eq_fit$S + eq_fit$I)
prev_inf  <- rowSums(eq_fit$I) / N_age_fit
prev_TS   <- eq_fit$TS / N_age_fit
prev_TT   <- eq_fit$TT / N_age_fit

model_df <- data.frame(
  age_mid  = age_mids,
  prev_inf = prev_inf,
  prev_TS  = prev_TS,
  prev_TT  = prev_TT
)

pct <- function(x) paste0(round(x * 100), "%")

## Plot A
beta_hat <- parms_fit$transRate
epsilon_hat <- parms_fit$epsilon

p_inf <- ggplot() +
  geom_ribbon(
    data = prev_data,
    aes(x = age_mid, ymin = lower, ymax = upper),
    fill = "#1565C0", alpha = 0.20
  ) +
  geom_point(
    data = prev_data,
    aes(x = age_mid, y = prev),
    colour = "#1565C0", size = 2.5, shape = 16
  ) +
  geom_line(
    data = model_df,
    aes(x = age_mid, y = prev_inf),
    colour = "#C62828", linewidth = 1.2
  ) +
  scale_x_continuous(breaks = c(0, 5, 10, 20, 40, 60, 80)) +
  scale_y_continuous(labels = pct, limits = c(0, NA)) +
  labs(
    x = "Age (years)",
    y = "Infection prevalence",
    title = sprintf(
      "Model fit to hypoendemic data   (\u03b2 = %.2f yr\u207b\u00b9, \u03b5 = %.2f)",
      beta_hat, epsilon_hat
    ),
    caption = "Blue points = data with 95% CI; Red line = 300-stage infection-history SIS model"
  ) +
  theme_minimal(base_size = 13)

## Plot B
seq_long <- rbind(
  data.frame(age_mid = age_mids, prev = prev_TS, sequela = "Trachomatous scarring (TS)"),
  data.frame(age_mid = age_mids, prev = prev_TT, sequela = "Trichiasis (TT)")
)

p_seq <- ggplot(seq_long, aes(x = age_mid, y = prev, colour = sequela)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(
    values = c(
      "Trachomatous scarring (TS)" = "#2E7D32",
      "Trichiasis (TT)" = "#E65100"
    ),
    name = NULL
  ) +
  scale_x_continuous(breaks = c(0, 5, 10, 20, 40, 60, 80)) +
  scale_y_continuous(labels = pct) +
  labs(
    x = "Age (years)",
    y = "Disease prevalence",
    title = "Age-profile of trachoma disease sequelae",
    caption = "Scarring accumulates with repeat infections (history stages j >= 2)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

## Plot C
selected_rungs <- c(1, 2, 3, 5, 10, 20, 30)
selected_rungs <- selected_rungs[selected_rungs >= 1 & selected_rungs <= parms_fit$m_hist]
include_other_rungs <- TRUE

comp_blocks <- list()

for (j in selected_rungs) {
  comp_blocks[[length(comp_blocks) + 1]] <- data.frame(
    age_mid = age_mids,
    proportion = eq_fit$S[, j] / N_age_fit,
    state = sprintf("S%02d susceptible", j)
  )

  comp_blocks[[length(comp_blocks) + 1]] <- data.frame(
    age_mid = age_mids,
    proportion = eq_fit$I[, j] / N_age_fit,
    state = sprintf("I%02d infected", j)
  )
}

if (include_other_rungs) {
  other_rungs <- setdiff(seq_len(parms_fit$m_hist), selected_rungs)

  if (length(other_rungs) > 0) {
    comp_blocks[[length(comp_blocks) + 1]] <- data.frame(
      age_mid = age_mids,
      proportion = rowSums(eq_fit$S[, other_rungs, drop = FALSE]) / N_age_fit,
      state = "S other rungs"
    )

    comp_blocks[[length(comp_blocks) + 1]] <- data.frame(
      age_mid = age_mids,
      proportion = rowSums(eq_fit$I[, other_rungs, drop = FALSE]) / N_age_fit,
      state = "I other rungs"
    )
  }
}

comp_long <- do.call(rbind, comp_blocks)
comp_long$state <- factor(comp_long$state, levels = unique(comp_long$state))

n_states <- length(levels(comp_long$state))
comp_cols <- setNames(grDevices::hcl.colors(n_states, palette = "Dynamic"), levels(comp_long$state))

p_comp <- ggplot(comp_long, aes(x = age_mid, y = proportion, fill = state)) +
  geom_area(position = "stack") +
  scale_fill_manual(values = comp_cols, name = NULL) +
  scale_x_continuous(breaks = c(0, 5, 10, 20, 40, 60, 80)) +
  scale_y_continuous(labels = pct) +
  labs(
    x = "Age (years)",
    y = "Proportion of age group",
    title = "Population composition by age (selected history rungs)",
    caption = "Selected S/I history rungs shown explicitly; remaining rungs optionally aggregated"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

## Plot D
beta_scenarios <- list(
  "Hypoendemic (\u03b2 = 1.8)" = 1.8,
  "Mesoendemic (\u03b2 = 5)" = 2.5,
  "Hyperendemic (\u03b2 = 15)" = 3.5,
  "Very high (\u03b2 = 30)" = 15
)

sens_list <- lapply(names(beta_scenarios), function(nm) {
  eq_k <- run_model_300(
    make_parms_300(
      transRate = beta_scenarios[[nm]],
      epsilon = parms_fit$epsilon
    ),
    years = 200,
    dt_out = 0.25,
    eq_window_years = 30
  )

  N_k <- pmax(rowSums(eq_k$S + eq_k$I), 1e-12)
  prev_k <- pmin(pmax(rowSums(eq_k$I) / N_k, 0), 1)

  data.frame(
    age_mid = age_mids,
    prev = prev_k,
    scenario = nm
  )
})

sens_df <- do.call(rbind, sens_list)
sens_df$scenario <- factor(sens_df$scenario, levels = names(beta_scenarios))

p_sens <- ggplot(sens_df, aes(x = age_mid, y = prev, colour = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_colour_brewer(palette = "OrRd", name = "Setting") +
  scale_x_continuous(breaks = c(0, 5, 10, 20, 40, 60, 80)) +
  scale_y_continuous(labels = pct) +
  labs(
    x = "Age (years)",
    y = "Infection prevalence",
    title = "Prevalence profiles across endemicity levels",
    caption = "Higher transmission shifts the peak to younger ages"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

print(p_inf)
print(p_seq)
print(p_comp)
print(p_sens)

## Plot E: epidemic trajectory
sol_df <- as.data.frame(eq_fit$sol)
names(sol_df)[1] <- "time"

S_cols <- grep("^S_a", names(sol_df))
I_cols <- grep("^I_a", names(sol_df))

traj_df <- data.frame(
  time = sol_df$time,
  prev = rowSums(sol_df[, I_cols, drop = FALSE]) /
    rowSums(sol_df[, c(S_cols, I_cols), drop = FALSE])
)

p_traj <- ggplot(traj_df, aes(x = time, y = prev)) +
  geom_line(colour = "#C62828", linewidth = 0.9) +
  scale_y_continuous(labels = pct) +
  xlim(199, 220) +
  labs(
    x = "Time (years)",
    y = "Overall infection prevalence",
    title = "Epidemic trajectory — convergence to endemic equilibrium",
    caption = sprintf("Seeded at t=0; equilibrium \u2248 %.1f%%",
                      100 * tail(traj_df$prev, 1))
  ) +
  theme_minimal(base_size = 13)

print(p_traj)

r0_fit <- calc_R0_ngm(parms_fit)
message(sprintf("R0 (NGM) for current fitted parameters: %.4f", r0_fit))

s_age_fit <- rowSums(eq_fit$S)
n_age_fit <- rowSums(eq_fit$S + eq_fit$I)
re_fit <- calc_Re_ngm(parms_fit, S_age = s_age_fit, N_age = n_age_fit)
message(sprintf("Re (NGM) for current fitted equilibrium state: %.4f", re_fit))
} # end sys.nframe() == 0 demo block