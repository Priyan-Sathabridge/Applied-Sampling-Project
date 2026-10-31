# 05_simulation.R
# Section E (and D7): repeat the whole sampling and estimation procedure
# R times for the stratified design and the SRSWOR benchmark.
#
# Each repetition: select the sample, simulate non-response, re-weight,
# estimate the mean and the proportion with their SEs and 95% CIs.
# The estimators are written out in closed form for speed; the first
# realised sample is checked against survey for both outcomes and designs.
#
# Run 01 to 03 first; this script reads data/population.csv and
# data/sample.csv (for the stratum sample sizes).

library(survey)
library(ggplot2)

set.seed(20261022)

R         <- 1000
resp_rate <- 0.90
rel_e     <- 0.05
z         <- qnorm(0.975)

population <- read.csv(file.path("data", "population.csv"))
smp        <- read.csv(file.path("data", "sample.csv"))

N         <- nrow(population)
true_mean <- mean(population$winnings_usd)
true_prop <- mean(population$below_ev)

# Stratum sizes and sample sizes, as used in Section C
alloc <- unique(smp[, c("region", "N_h", "n_h")])
n_total <- sum(alloc$n_h)
rows_by_region <- split(seq_len(N), population$region)

dir.create("output",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ---- Estimators -----------------------------------------------------------

# Stratified mean and variance from respondents:
# ybar_st = sum W_h ybar_h,  v = sum W_h^2 (1 - r_h / N_h) s_h^2 / r_h
est_strat <- function(y, region, N_h) {
  W_h  <- N_h / sum(N_h)
  ybar <- tapply(y, region, mean)[names(N_h)]
  s2   <- tapply(y, region, var)[names(N_h)]
  r_h  <- tapply(y, region, length)[names(N_h)]
  stopifnot(all(is.finite(r_h)), all(r_h >= 2))
  c(est = sum(W_h * ybar),
    se  = sqrt(sum(W_h^2 * (1 - r_h / N_h) * s2 / r_h)))
}

# SRSWOR mean and variance from respondents: v = (1 - r / N) s^2 / r
est_srs <- function(y, N) {
  r <- length(y)
  stopifnot(r >= 2)
  c(est = mean(y), se = sqrt((1 - r / N) * var(y) / r))
}

N_h <- setNames(alloc$N_h, alloc$region)

# ---- One repetition -------------------------------------------------------

one_rep <- function() {
  # Stratified: SRSWOR of n_h within each region, then non-response
  sel <- unlist(lapply(alloc$region, function(h) {
    rows <- rows_by_region[[h]]
    rows[sample.int(length(rows), alloc$n_h[alloc$region == h])]
  }))
  sel  <- sel[rbinom(length(sel), 1, resp_rate) == 1]
  d    <- population[sel, ]
  st_m <- est_strat(d$winnings_usd, d$region, N_h)
  st_p <- est_strat(d$below_ev,     d$region, N_h)

  # SRSWOR benchmark of the same selected size, same response rate
  sel2 <- sample.int(N, n_total)
  sel2 <- sel2[rbinom(length(sel2), 1, resp_rate) == 1]
  d2   <- population[sel2, ]
  sr_m <- est_srs(d2$winnings_usd, N)
  sr_p <- est_srs(d2$below_ev,     N)

  rbind(
    data.frame(design = "Stratified (Neyman)", parameter = "Mean winnings",
               est = st_m[["est"]], se = st_m[["se"]], n_resp = nrow(d)),
    data.frame(design = "SRSWOR", parameter = "Mean winnings",
               est = sr_m[["est"]], se = sr_m[["se"]], n_resp = nrow(d2)),
    data.frame(design = "Stratified (Neyman)", parameter = "Proportion below opening EV",
               est = st_p[["est"]], se = st_p[["se"]], n_resp = nrow(d)),
    data.frame(design = "SRSWOR", parameter = "Proportion below opening EV",
               est = sr_p[["est"]], se = sr_p[["se"]], n_resp = nrow(d2))
  )
}

# ---- Check the closed-form estimator against the survey package -----------

resp <- smp[smp$responded == 1, ]
des  <- svydesign(ids = ~1, strata = ~stratum_id, weights = ~w_final,
                  fpc = ~N_h, data = resp)
chk_survey <- svymean(~winnings_usd, des)
chk_manual <- est_strat(resp$winnings_usd, resp$region, N_h)
stopifnot(abs(coef(chk_survey) - chk_manual[["est"]]) < 1e-6,
          abs(SE(chk_survey)   - chk_manual[["se"]])  < 1e-6)

# Check both outcomes and the benchmark without drawing any random numbers.
for (variable in c("winnings_usd", "below_ev")) {
  a <- svymean(reformulate(variable), des)
  b <- est_strat(resp[[variable]], resp$region, N_h)
  stopifnot(abs(coef(a) - b["est"]) < 1e-6,
            abs(SE(a) - b["se"]) < 1e-6)
  ds <- svydesign(ids = ~1, weights = ~I(N / nrow(resp)),
                  fpc = ~I(rep(N, nrow(resp))), data = resp)
  a <- svymean(reformulate(variable), ds)
  b <- est_srs(resp[[variable]], N)
  stopifnot(abs(coef(a) - b["est"]) < 1e-6,
            abs(SE(a) - b["se"]) < 1e-6)
}

# ---- E1 and E2: repeat R times --------------------------------------------

reps <- do.call(rbind, lapply(seq_len(R), function(i) cbind(rep = i, one_rep())))
reps$lower <- reps$est - z * reps$se
reps$upper <- reps$est + z * reps$se
reps$true  <- ifelse(reps$parameter == "Mean winnings", true_mean, true_prop)

write.csv(reps, file.path("output", "E_replicates.csv"), row.names = FALSE)

# ---- E3 to E5 and D7: performance measures --------------------------------

perf <- do.call(rbind, lapply(split(reps, list(reps$design, reps$parameter)),
  function(d) {
    true     <- d$true[1]
    coverage <- mean(d$lower <= true & true <= d$upper)
    data.frame(
      design         = d$design[1],
      parameter      = d$parameter[1],
      true_value     = true,
      mean_estimate  = mean(d$est),
      bias           = mean(d$est) - true,
      mc_se_bias     = sd(d$est) / sqrt(nrow(d)),
      avg_estimated_variance = mean(d$se^2),
      relative_bias  = (mean(d$est) - true) / true,
      empirical_se   = sd(d$est),
      avg_estimated_se = mean(d$se),
      rmse           = sqrt(mean((d$est - true)^2)),
      coverage       = coverage,
      mc_se_coverage = sqrt(coverage * (1 - coverage) / nrow(d)),
      avg_ci_width   = mean(d$upper - d$lower),
      pct_meeting_5pct = mean(z * d$se <= rel_e * true)
    )
  }))
rownames(perf) <- NULL

# Empirical design effect: variance of the stratified estimates over
# variance of the SRSWOR estimates
perf$empirical_deff <- NA
for (p in unique(perf$parameter)) {
  v_st  <- perf$empirical_se[perf$parameter == p & perf$design == "Stratified (Neyman)"]^2
  v_srs <- perf$empirical_se[perf$parameter == p & perf$design == "SRSWOR"]^2
  perf$empirical_deff[perf$parameter == p & perf$design == "Stratified (Neyman)"] <- v_st / v_srs
  perf$empirical_deff[perf$parameter == p & perf$design == "SRSWOR"] <- 1
}

print(perf, digits = 4)
write.csv(perf, file.path("output", "E_performance.csv"), row.names = FALSE)

# ---- Figures --------------------------------------------------------------

p_dist <- ggplot(reps, aes(est, fill = design)) +
  geom_density(alpha = 0.4) +
  geom_vline(aes(xintercept = true), linetype = "dashed") +
  facet_wrap(~parameter, scales = "free") +
  labs(x = "Estimate", y = "Density", fill = "Design",
       title = paste("Sampling distribution of the estimates over", R, "repetitions")) +
  theme_minimal()
ggsave(file.path("figures", "E_sampling_distributions.png"), p_dist, width = 9, height = 4.5)

# First 100 confidence intervals for the mean under each design
ci_plot <- reps[reps$parameter == "Mean winnings" & reps$rep <= 100, ]
ci_plot$covers <- ci_plot$lower <= ci_plot$true & ci_plot$true <= ci_plot$upper
p_ci <- ggplot(ci_plot, aes(rep, est, ymin = lower, ymax = upper, colour = covers)) +
  geom_pointrange(size = 0.15) +
  geom_hline(yintercept = true_mean, linetype = "dashed") +
  facet_wrap(~design) +
  scale_colour_manual(values = c(`TRUE` = "grey40", `FALSE` = "red")) +
  labs(x = "Repetition", y = "Mean winnings (USD)", colour = "Covers true mean",
       title = "95% confidence intervals for mean winnings, first 100 repetitions") +
  theme_minimal()
ggsave(file.path("figures", "E_confidence_intervals.png"), p_ci, width = 9, height = 4.5)
