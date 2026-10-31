# 04_estimation.R
# Section D: design-based estimation with the survey package.
# Stratified sample (Neyman allocation by region) against an SRSWOR
# benchmark of the same size.
#
# Run 01 to 03 first; this script reads data/population.csv and
# data/sample.csv. D7 (empirical against estimated SE) needs repeated
# samples and is done in the Section E simulation script.

library(survey)

population <- read.csv(file.path("data", "population.csv"))
smp        <- read.csv(file.path("data", "sample.csv"))

N         <- nrow(population)
rel_e     <- 0.05
resp_rate <- 0.90

true_mean <- mean(population$winnings_usd)
true_prop <- mean(population$below_ev)

dir.create("output", showWarnings = FALSE)

# Estimates are made from respondents only, using the
# non-response-adjusted weights N_h / r_h.
resp <- smp[smp$responded == 1, ]

# ---- D1: survey design objects --------------------------------------------

# Conditional on each realised r_h, independent equal-probability response
# leaves an SRSWOR of r_h in region h. Thus the respondent-level fpc is valid
# under the simulated response mechanism; this does not cover informative response.
# Stratified: strata = region, weights = N_h / r_h, fpc = N_h
des_str <- svydesign(ids = ~1, strata = ~stratum_id, weights = ~w_final,
                     fpc = ~N_h, data = resp)

# SRSWOR benchmark: select the same number of contestants as the
# stratified sample, apply the same response rate, weights N / r
set.seed(20261021)
srs_sel <- population[sample.int(N, nrow(smp)), ]
srs_sel$responded <- rbinom(nrow(srs_sel), 1, resp_rate)
srs <- srs_sel[srs_sel$responded == 1, ]
srs$N_pop   <- N
srs$w_final <- N / nrow(srs)

des_srs <- svydesign(ids = ~1, weights = ~w_final, fpc = ~N_pop, data = srs)

# Store selection probabilities and weights for every benchmark unit too.
srs_sel$pi_i <- nrow(srs_sel) / N
srs_sel$w_design <- 1 / srs_sel$pi_i
srs_sel$w_final <- ifelse(srs_sel$responded == 1, N / nrow(srs), 0)
srs_sel[srs_sel$responded == 0, c("winnings_usd", "below_ev", "deal",
  "round_ended", "risk_aversion", "wealth_usd")] <- NA
write.csv(srs_sel, file.path("data", "sample_srswor.csv"), row.names = FALSE)

# ---- D2 and D3: estimates -------------------------------------------------

summarise_est <- function(design, var, true_value, design_name) {
  est <- svymean(as.formula(paste0("~", var)), design, deff = TRUE)
  ci  <- confint(est)
  se  <- as.numeric(SE(est))
  data.frame(
    design          = design_name,
    parameter       = var,
    n_respondents   = nrow(design$variables),
    estimate        = as.numeric(coef(est)),
    variance        = se^2,
    se              = se,
    ci_lower        = ci[1, 1],
    ci_upper        = ci[1, 2],
    cv              = se / as.numeric(coef(est)),
    deff            = as.numeric(deff(est)),
    half_width      = (ci[1, 2] - ci[1, 1]) / 2,
    target_hw       = rel_e * true_value,
    true_value      = true_value,
    ci_covers_truth = ci[1, 1] <= true_value & true_value <= ci[1, 2]
  )
}

results <- rbind(
  summarise_est(des_str, "winnings_usd", true_mean, "Stratified (Neyman)"),
  summarise_est(des_srs, "winnings_usd", true_mean, "SRSWOR"),
  summarise_est(des_str, "below_ev",     true_prop, "Stratified (Neyman)"),
  summarise_est(des_srs, "below_ev",     true_prop, "SRSWOR")
)
results$meets_5pct <- results$half_width <= results$target_hw
print(results, digits = 4)

# Estimates within each region
print(svyby(~winnings_usd + below_ev, ~region, des_str, svymean))

# ---- D4: inclusion probabilities and weights ------------------------------
# pi_i = n_h / N_h is the chance a contestant in region h is selected;
# w_design = 1 / pi_i is how many population members each selected
# contestant represents. After non-response, w_final = N_h / r_h spreads
# each region's total over its respondents only.

weights_table <- unique(smp[, c("region", "N_h", "n_h", "pi_i", "w_design",
                                "r_h", "w_final")])
weights_table <- weights_table[weights_table$w_final > 0, ]
print(weights_table)

# ---- D5: naive against design-adjusted against true -----------------------
# The naive mean ignores the weights. The US is sampled at a higher rate
# than its population share, so the naive mean is pulled towards the US.

comparison <- data.frame(
  parameter       = c("Mean winnings (USD)", "Proportion below opening EV"),
  naive           = c(mean(resp$winnings_usd), mean(resp$below_ev)),
  design_adjusted = results$estimate[results$design == "Stratified (Neyman)"],
  true_value      = c(true_mean, true_prop)
)
comparison$naive_error    <- comparison$naive - comparison$true_value
comparison$adjusted_error <- comparison$design_adjusted - comparison$true_value
print(comparison, digits = 4)

# ---- D6 and D8: precision of the two designs ------------------------------

precision <- data.frame(
  parameter = c("Mean winnings", "Proportion below opening EV"),
  se_strat  = results$se[results$design == "Stratified (Neyman)"],
  se_srswor = results$se[results$design == "SRSWOR"]
)
precision$variance_ratio <- precision$se_strat^2 / precision$se_srswor^2
print(precision, digits = 4)

write.csv(results,    file.path("output", "D_estimates.csv"),      row.names = FALSE)
write.csv(comparison, file.path("output", "D_naive_vs_design.csv"), row.names = FALSE)

write.csv(weights_table, "output/D_weights.csv", row.names = FALSE)
write.csv(precision, "output/D_precision.csv", row.names = FALSE)
