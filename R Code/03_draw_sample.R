# 03_draw_sample.R
# Section C: draw one stratified random sample (Neyman allocation by region)
# from the finite population, simulate non-response, and attach inclusion
# probabilities and sampling weights.
#
# Run 01_simulate_population.R first; this script reads data/population.csv
# (the full finite population, already saved as CSV by that script).

set.seed(20261020)

population <- read.csv(file.path("data", "population.csv"))

N         <- nrow(population)
z         <- qnorm(0.975)
rel_e     <- 0.05
resp_rate <- 0.90   # anticipated response rate, as in 02_sample_size.R

# ---- C1: final sample size (same calculation as Section B) ----------------

strata <- do.call(rbind, lapply(split(population, population$region),
  function(d) data.frame(region = d$region[1],
                         N_h    = nrow(d),
                         S_h    = sd(d$winnings_usd))))
strata$stratum_id <- seq_len(nrow(strata))

e_mean <- rel_e * mean(population$winnings_usd)
V      <- (e_mean / z)^2
n      <- sum(strata$N_h * strata$S_h)^2 /
          (N^2 * V + sum(strata$N_h * strata$S_h^2))

# Neyman allocation, then inflate for non-response
strata$n_h_resp <- ceiling(n * strata$N_h * strata$S_h / sum(strata$N_h * strata$S_h))
strata$n_h      <- ceiling(strata$n_h_resp / resp_rate)
stopifnot(all(strata$n_h <= strata$N_h))

print(strata)
cat("Respondents needed:", sum(strata$n_h_resp),
    "| Units selected:", sum(strata$n_h), "\n")

# ---- C3: select the sample ------------------------------------------------
# Within each stratum, draw n_h contestants by simple random sampling
# without replacement, independently across strata.

selected <- unlist(lapply(seq_len(nrow(strata)), function(i) {
  rows <- which(population$region == strata$region[i])
  rows[sample.int(length(rows), strata$n_h[i])]
}))

smp <- population[selected, ]
smp <- merge(smp, strata[, c("region", "stratum_id", "N_h", "n_h")],
             by = "region", sort = FALSE)

# ---- C4 and C5: inclusion probabilities and weights -----------------------
# Every contestant in stratum h has inclusion probability n_h / N_h > 0,
# so the design weight is N_h / n_h.

smp$pi_i     <- smp$n_h / smp$N_h
smp$w_design <- 1 / smp$pi_i
stopifnot(all(smp$pi_i > 0))

# ---- Non-response ---------------------------------------------------------
# Each selected contestant responds independently with probability
# resp_rate (assumed missing completely at random within strata). The
# weights of respondents are scaled up so they again add to N_h.

smp$responded <- rbinom(nrow(smp), 1, resp_rate)
r_h <- tapply(smp$responded, smp$region, sum)
stopifnot(all(r_h >= 2)) # variance estimation needs at least two respondents per stratum
smp$r_h     <- as.integer(r_h[smp$region])
smp$w_final <- ifelse(smp$responded == 1, smp$N_h / smp$r_h, 0)

# Survey variables are only observed for respondents; frame variables
# (region, age, gender, opening EV) are known for everyone selected.
survey_vars <- c("winnings_usd", "below_ev", "deal", "round_ended",
                 "risk_aversion", "wealth_usd")
smp[smp$responded == 0, survey_vars] <- NA

# Checks: design weights add to N_h in each stratum, final weights to N
stopifnot(all(abs(tapply(smp$w_design, smp$region, sum) -
                  strata$N_h[match(names(r_h), strata$region)]) < 1e-6))
stopifnot(abs(sum(smp$w_final) - N) < 1e-6)

# ---- C6: save the realised sample -----------------------------------------

smp <- smp[order(smp$stratum_id, smp$id), c(
  "id", "stratum_id", "region", "N_h", "n_h", "pi_i", "w_design",
  "responded", "r_h", "w_final",
  "winnings_usd", "below_ev", "deal", "round_ended",
  "opening_ev_usd", "age", "gender", "risk_aversion", "wealth_usd")]

print(aggregate(cbind(selected = 1, responded) ~ region, data = smp, FUN = sum))

write.csv(smp, file.path("data", "sample.csv"), row.names = FALSE)
