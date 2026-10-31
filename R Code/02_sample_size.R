# 02_sample_size.R
# Section B: sample-size determination.
# Primary design : stratified random sampling by region, Neyman allocation
# Benchmark      : SRSWOR
# Primary        : mean winnings per contestant (USD)
# Secondary      : proportion who won less than the opening expected value
#
# Run 01_simulate_population.R first; this script reads data/population.csv.

library(ggplot2)

population <- read.csv(file.path("data", "population.csv"))

N         <- nrow(population)
z         <- qnorm(0.975)   # 95% confidence
rel_e     <- 0.05           # half-width = 5% of the baseline value
resp_rate <- 0.90           # anticipated response rate (assumed)

dir.create("figures", showWarnings = FALSE)

# ---- Prior information ----------------------------------------------------
# The baseline values and stratum variances are treated as known prior
# information, as the project brief allows.

strata <- do.call(rbind, lapply(split(population, population$region),
  function(d) data.frame(
    region = d$region[1],
    N_h    = nrow(d),
    S_h    = sd(d$winnings_usd),
    P_h    = mean(d$below_ev)
  )))
strata$W_h  <- strata$N_h / N
strata$Sp_h <- sqrt(strata$N_h / (strata$N_h - 1) * strata$P_h * (1 - strata$P_h))

mu0 <- mean(population$winnings_usd)
P0  <- mean(population$below_ev)
S_y <- sd(population$winnings_usd)
S_p <- sqrt(N / (N - 1) * P0 * (1 - P0))

# ---- Sample-size formulas -------------------------------------------------

# SRSWOR: n0 = (z S / e)^2, then the finite population correction
# n = n0 / (1 + n0 / N)
n_srs <- function(e, S, N) {
  n0 <- (z * S / e)^2
  n0 / (1 + n0 / N)
}

# Stratified with Neyman allocation (n_h proportional to N_h S_h), with fpc:
# n = (sum N_h S_h)^2 / (N^2 V + sum N_h S_h^2),  where V = (e / z)^2
# For tight precision, unconstrained Neyman allocation can exceed N_h.
# Set such strata to census, then optimise allocation over the remaining
# strata. Census strata contribute zero variance; population weights remain
# N_h / N for the original population. Returns a continuous lower bound.
capped_neyman <- function(e, N_h, S_h) {
  stopifnot(e > 0, all(N_h > 1), all(S_h > 0))
  V <- (e / z)^2
  N <- sum(N_h)
  nh <- numeric(length(N_h))
  free <- rep(TRUE, length(N_h))
  repeat {
    c_mult <- sum(N_h[free] * S_h[free]) /
      (N^2 * V + sum(N_h[free] * S_h[free]^2))
    nh[free] <- c_mult * N_h[free] * S_h[free]
    hit <- free & nh >= N_h
    if (!any(hit)) break
    nh[hit] <- N_h[hit]
    free[hit] <- FALSE
    if (!any(free)) break
  }
  nh
}
n_neyman <- function(e, N_h, S_h) sum(capped_neyman(e, N_h, S_h))

neyman_alloc <- function(n, N_h, S_h) n * N_h * S_h / sum(N_h * S_h)

# Variance of the estimated mean under each design
var_strat <- function(n_h, N_h, S_h) {
  W_h <- N_h / sum(N_h)
  sum(W_h^2 * (1 - n_h / N_h) * S_h^2 / n_h)
}
var_srs <- function(n, S, N) (1 - n / N) * S^2 / n

# ---- B3 to B5: minimum sample size and allocation -------------------------

e_mean <- rel_e * mu0
e_prop <- rel_e * P0

n_req <- data.frame(
  parameter = c("Mean winnings", "Proportion below opening EV"),
  baseline  = c(mu0, P0),
  half_width = c(e_mean, e_prop),
  n_srswor  = ceiling(c(n_srs(e_mean, S_y, N), n_srs(e_prop, S_p, N))),
  n_strat   = ceiling(c(n_neyman(e_mean, strata$N_h, strata$S_h),
                        n_neyman(e_prop, strata$N_h, strata$Sp_h)))
)
print(n_req)

# The allocation is driven by the primary parameter (the mean)
strata$n_h       <- ceiling(neyman_alloc(n_neyman(e_mean, strata$N_h, strata$S_h),
                                         strata$N_h, strata$S_h))
strata$n_h_final <- ceiling(strata$n_h / resp_rate)   # allow for non-response
stopifnot(all(strata$n_h_final <= strata$N_h))

n_total       <- sum(strata$n_h)
n_total_final <- sum(strata$n_h_final)

# Design effect of the stratified design relative to SRSWOR of the same size
deff_mean <- var_strat(strata$n_h, strata$N_h, strata$S_h) /
             var_srs(n_total, S_y, N)

# Check that the same sample also meets the 5% bound for the proportion
hw_prop <- z * sqrt(var_strat(strata$n_h, strata$N_h, strata$Sp_h))
deff_prop <- var_strat(strata$n_h, strata$N_h, strata$Sp_h) /
             var_srs(n_total, S_p, N)

print(strata)
cat("Total sample before non-response:", n_total, "\n")
cat("Total sample to select (response rate", resp_rate, "):", n_total_final, "\n")
cat("Design effect, mean:", round(deff_mean, 3), "\n")
cat("Half-width for the proportion with this sample:", round(hw_prop, 4),
    "(target", round(e_prop, 4), ")\n")
cat("Design effect, proportion:", round(deff_prop, 3), "\n")

# ---- B6: sample size against margin of error ------------------------------
# Margin of error expressed as a percentage of the baseline mean.

rel_grid <- seq(0.01, 0.10, by = 0.0025)

curve_b6 <- rbind(
  data.frame(rel_e = rel_grid, design = "Stratified (Neyman)",
             n = sapply(rel_grid * mu0, n_neyman, N_h = strata$N_h, S_h = strata$S_h)),
  data.frame(rel_e = rel_grid, design = "SRSWOR",
             n = sapply(rel_grid * mu0, n_srs, S = S_y, N = N))
)

p_b6 <- ggplot(curve_b6, aes(100 * rel_e, n, colour = design)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 5, linetype = "dashed") +
  labs(x = "Half-width of 95% CI (% of baseline mean)",
       y = "Required sample size",
       colour = "Design",
       title = "Required sample size against margin of error, mean winnings") +
  theme_minimal()
ggsave(file.path("figures", "B6_n_vs_margin.png"), p_b6, width = 7, height = 4.5)

# ---- B7 and B8: sensitivity to the assumed mean ---------------------------
# Held fixed: N, stratum sizes N_h, stratum SDs S_h, capacity-constrained Neyman allocation rule,
# 95% confidence. Only the assumed mean changes, so the absolute margin
# (a percentage of the assumed mean) shrinks as the assumed mean falls.

k_vals <- c(0.25, 0.50, 0.75, 1.00)

curve_b7 <- do.call(rbind, lapply(k_vals, function(k)
  data.frame(rel_e = rel_grid,
             assumed = paste0(100 * k, "% of baseline"),
             n = sapply(rel_grid * k * mu0, n_neyman,
                        N_h = strata$N_h, S_h = strata$S_h))))
curve_b7$assumed <- factor(curve_b7$assumed,
                           levels = paste0(100 * k_vals, "% of baseline"))

p_b7 <- ggplot(curve_b7, aes(100 * rel_e, n, colour = assumed)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 5, linetype = "dashed") +
  geom_hline(yintercept = N, linetype = "dotted") +
  labs(x = "Half-width of 95% CI (% of assumed mean)",
       y = "Required sample size (stratified, Neyman)",
       colour = "Assumed mean",
       title = "Sensitivity of required sample size to the assumed mean") +
  theme_minimal()
ggsave(file.path("figures", "B7_sensitivity_mean.png"), p_b7, width = 7, height = 4.5)

table_b7 <- data.frame(
  assumed_mean = k_vals * mu0,
  pct_baseline = 100 * k_vals,
  n_strat = ceiling(sapply(rel_e * k_vals * mu0, n_neyman,
                           N_h = strata$N_h, S_h = strata$S_h)),
  n_srswor = ceiling(sapply(rel_e * k_vals * mu0, n_srs, S = S_y, N = N))
)
print(table_b7)

# ---- B9: proportion, P from 0.10 to 0.90 ----------------------------------
# Held fixed: N, the absolute half-width e_prop (5% of the baseline P),
# 95% confidence, and the Neyman design effect at the baseline. The
# variance P(1 - P) is largest at P = 0.5, so n peaks there.

deff_p_fixed <- sum(strata$W_h * strata$Sp_h)^2 / S_p^2   # Neyman, ignoring fpc

P_grid <- seq(0.10, 0.90, by = 0.01)
S_grid <- sqrt(N / (N - 1) * P_grid * (1 - P_grid))

curve_b9 <- rbind(
  data.frame(P = P_grid, design = "SRSWOR",
             n = n_srs(e_prop, S_grid, N)),
  data.frame(P = P_grid, design = "Stratified (Neyman)",
             n = n_srs(e_prop, sqrt(deff_p_fixed) * S_grid, N))
)

p_b9 <- ggplot(curve_b9, aes(P, n, colour = design)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = P0, linetype = "dotted") +
  labs(x = "Assumed proportion P",
       y = "Required sample size",
       colour = "Design",
       title = paste0("Required sample size against P (half-width fixed at ",
                      round(e_prop, 4), ")")) +
  theme_minimal()
ggsave(file.path("figures", "B9_n_vs_P.png"), p_b9, width = 7, height = 4.5)

table_b9 <- data.frame(P = c(0.10, 0.25, 0.50, P0, 0.90))
table_b9$PQ       <- table_b9$P * (1 - table_b9$P)
table_b9$n_srswor <- ceiling(n_srs(e_prop, sqrt(N / (N - 1) * table_b9$PQ), N))
print(table_b9)

# Save planning outputs so the report uses the executed calculations.
dir.create("output", showWarnings = FALSE)
write.csv(n_req, "output/B_sample_size_summary.csv", row.names = FALSE)
write.csv(strata, "output/B_allocation.csv", row.names = FALSE)
write.csv(table_b7, "output/B_sensitivity.csv", row.names = FALSE)
write.csv(table_b9, "output/B_proportion_sensitivity.csv", row.names = FALSE)

# Verify feasibility and target variance throughout BOTH sensitivity grids.
for (k in k_vals) for (margin in rel_grid) {
  a <- capped_neyman(margin * k * mu0, strata$N_h, strata$S_h)
  stopifnot(all(a > 0), all(a <= strata$N_h),
    abs(var_strat(a, strata$N_h, strata$S_h) -
        (margin * k * mu0 / z)^2) < 1e-5)
}
