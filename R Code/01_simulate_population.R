# 01_simulate_population.R
# Simulate a finite population of 10,000 Deal or No Deal contestants.
# Each game follows a stylised region-specific board and opening schedule,
# and all winnings are converted to USD.
#
# Strata  : region (UK, US, AU)
# Primary : mean winnings per contestant (USD)
# Second. : proportion of contestants who won less than the opening
#           expected value of their region's board

set.seed(20261019)

N <- 10000

# ---- Region settings ------------------------------------------------------

# Prize boards in local currency.
# UK: Channel 4 version, 22 boxes, 1p to GBP 250,000
# US: NBC version, 26 cases, $0.01 to $1,000,000
# AU: assumed high-stakes board; not claimed to reproduce a verified series
boards <- list(
  UK = c(0.01, 0.10, 0.50, 1, 5, 10, 50, 100, 250, 500, 750,
         1000, 3000, 5000, 10000, 15000, 20000, 35000, 50000,
         75000, 100000, 250000),
  US = c(0.01, 1, 5, 10, 25, 50, 75, 100, 200, 300, 400, 500, 750,
         1000, 5000, 10000, 25000, 50000, 75000, 100000, 200000,
         300000, 400000, 500000, 750000, 1000000),
  AU = c(0.05, 0.25, 0.50, 0.75, 1, 5, 10, 25, 50, 75, 100, 250, 500,
         1000, 2500, 5000, 7500, 10000, 25000, 50000, 75000, 100000,
         250000, 500000, 1000000, 2000000)
)

# Boxes opened before each banker offer. Every schedule leaves the
# contestant's box and one other box unopened after the final round.
rounds <- list(
  UK = c(5, 3, 3, 3, 3, 3),
  US = c(6, 5, 4, 3, 2, 1, 1, 1, 1),
  AU = c(6, 5, 4, 3, 2, 1, 1, 1, 1)
)

# Fixed modelling conversion factors to USD; not dated market quotations
fx_to_usd <- c(UK = 1.27, US = 1.00, AU = 0.66)

# Stratum population sizes (assumed; must sum to N)
N_h <- c(UK = 4000, US = 4000, AU = 2000)
stopifnot(sum(N_h) == N)

# ---- Game mechanics -------------------------------------------------------

# Banker offer as a fraction of the expected value of the unopened boxes.
# Low early on, rising towards the expected value in the last round,
# with a little random variation from game to game.
banker_fraction <- function(r, R) {
  base <- 0.15 + 0.85 * ((r - 1) / (R - 1))^1.5
  min(base * exp(rnorm(1, 0, 0.10)), 1.10)
}

# Certainty equivalent of a gamble paying each value in x with equal
# probability, under CRRA utility with risk aversion rho and wealth w.
cert_equiv <- function(x, rho, w) {
  if (abs(rho - 1) < 1e-6) {
    exp(mean(log(w + x))) - w
  } else {
    mean((w + x)^(1 - rho))^(1 / (1 - rho)) - w
  }
}

# Play one game. The contestant accepts the first offer that is at least
# the certainty equivalent of the unopened boxes (a myopic rule: future
# offers are not taken into account).
play_game <- function(board, opens, rho, w) {
  boxes  <- sample(board)      # prizes randomly assigned to boxes
  own    <- boxes[1]           # contestant's box
  others <- boxes[-1]          # opened in this (random) order
  R <- length(opens)

  for (r in seq_len(R)) {
    others    <- others[-seq_len(opens[r])]
    remaining <- c(own, others)
    offer     <- banker_fraction(r, R) * mean(remaining)
    if (offer >= cert_equiv(remaining, rho, w)) {
      return(c(winnings = offer, deal = 1, round_ended = r))
    }
  }
  c(winnings = own, deal = 0, round_ended = R)
}

# ---- Contestants ----------------------------------------------------------

region <- rep(names(N_h), times = N_h)

age    <- round(pmin(pmax(rnorm(N, 40, 12), 18), 75))
gender <- sample(c("Female", "Male"), N, replace = TRUE)

# Relative risk aversion, rising slightly with age
risk_aversion <- pmax(rnorm(N, 0.4 + 0.01 * (age - 40), 0.3), 0)

# Baseline wealth in USD, converted to local currency for the game
wealth_usd <- rlnorm(N, meanlog = log(20000), sdlog = 0.8)

games <- t(vapply(seq_len(N), function(i) {
  reg <- region[i]
  play_game(boards[[reg]], rounds[[reg]], risk_aversion[i],
            wealth_usd[i] / fx_to_usd[[reg]])
}, numeric(3)))

opening_ev_usd <- sapply(boards, mean) * fx_to_usd

population <- data.frame(
  id             = sprintf("C%05d", seq_len(N)),
  region         = region,
  age            = age,
  gender         = gender,
  risk_aversion  = round(risk_aversion, 3),
  wealth_usd     = round(wealth_usd, 2),
  deal           = games[, "deal"],
  round_ended    = games[, "round_ended"],
  winnings_usd   = round(games[, "winnings"] * fx_to_usd[region], 2),
  opening_ev_usd = round(opening_ev_usd[region], 2)
)
population$below_ev <- as.integer(population$winnings_usd <
                                    population$opening_ev_usd)

# ---- True finite-population values ----------------------------------------

true_mean <- mean(population$winnings_usd)
true_prop <- mean(population$below_ev)

strata_summary <- do.call(rbind, lapply(split(population, population$region),
  function(d) data.frame(
    region    = d$region[1],
    N_h       = nrow(d),
    W_h       = nrow(d) / N,
    mean_h    = mean(d$winnings_usd),
    sd_h      = sd(d$winnings_usd),
    prop_h    = mean(d$below_ev),
    deal_rate = mean(d$deal)
  )))

print(strata_summary)
cat("True mean winnings (USD):", round(true_mean, 2), "\n")
cat("True proportion below opening EV:", round(true_prop, 4), "\n")

dir.create("data", showWarnings = FALSE)
write.csv(population, file.path("data", "population.csv"), row.names = FALSE)

stopifnot(nrow(population) == 10000, !anyDuplicated(population$id),
          all(is.finite(population$winnings_usd)),
          all(population$winnings_usd > 0), all(population$below_ev %in% 0:1))
