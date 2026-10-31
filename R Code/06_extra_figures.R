# 06_extra_figures.R
# Four additional figures for the report:
#   A_winnings_by_region.png      - distribution of winnings by region (A5, A6)
#   A_banker_offer_curve.png      - banker offer as % of expected value (A5)
#   B_allocation_comparison.png   - n needed under each allocation (B5)
#   E_half_width_distribution.png - realised CI half-widths vs the 5% target (E, D8)
#
# Run 01 to 05 first; this script reads data/population.csv and
# output/E_replicates.csv. It draws no random numbers.

library(ggplot2)

population <- read.csv(file.path("data", "population.csv"))
reps       <- read.csv(file.path("output", "E_replicates.csv"))

N     <- nrow(population)
z     <- qnorm(0.975)
rel_e <- 0.05

dir.create("figures", showWarnings = FALSE)

# ---- 1. Winnings by region ------------------------------------------------
# Log scale, because winnings run from cents to over a million dollars.
# Diamonds mark each region's mean; crosses mark the opening expected value.

region_means <- aggregate(winnings_usd ~ region, population, mean)
opening_ev   <- unique(population[, c("region", "opening_ev_usd")])

p1 <- ggplot(population, aes(region, winnings_usd)) +
  geom_violin(fill = "grey85", colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.4, outlier.alpha = 0.3) +
  geom_point(data = region_means, shape = 23, size = 3, fill = "steelblue") +
  geom_point(data = opening_ev, aes(y = opening_ev_usd), shape = 4, size = 3,
             stroke = 1.2, colour = "firebrick") +
  # amounts under $1 keep their cents, so the bottom labels are not all "$0"
  scale_y_log10(labels = function(x) ifelse(x < 1, paste0("$", x), scales::dollar(x)),
                breaks = 10^(-2:6)) +
  labs(x = "Region (stratum)", y = "Winnings (USD, log scale)",
       title = "Distribution of winnings by region",
       subtitle = "Blue diamond = mean winnings; red cross = opening expected value of the board") +
  theme_minimal()
ggsave(file.path("figures", "A_winnings_by_region.png"), p1, width = 7, height = 4.5)

# ---- 2. Banker offer curve ------------------------------------------------
# Same rule as banker_fraction() in 01_simulate_population.R:
# base fraction 0.15 + 0.85 ((r - 1) / (R - 1))^1.5, multiplied by
# exp(N(0, 0.10^2)) noise and capped at 1.10. The band shows the middle 95%
# of that noise.

offer_curve <- function(show, R) {
  r    <- seq_len(R)
  base <- 0.15 + 0.85 * ((r - 1) / (R - 1))^1.5
  data.frame(show  = show, round = r, base = base,
             lower = base * exp(-z * 0.10),
             upper = pmin(base * exp(z * 0.10), 1.10))
}
offers <- rbind(offer_curve("UK (6 rounds)", 6),
                offer_curve("US and Australia (9 rounds)", 9))

p2 <- ggplot(offers, aes(round, base, colour = show, fill = show)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.15)) +
  scale_x_continuous(breaks = 1:9) +
  labs(x = "Round", y = "Offer as % of expected value of unopened boxes",
       colour = "Show", fill = "Show",
       title = "Banker offers rise towards the expected value as the game goes on",
       subtitle = "Line = median offer; band = middle 95% of offers; dashed line = expected value") +
  theme_minimal()
ggsave(file.path("figures", "A_banker_offer_curve.png"), p2, width = 7, height = 4.5)

# ---- 3. Required sample size under each allocation ------------------------
# For allocation shares a_h (n_h = n a_h), the variance of the stratified
# mean with fpc is  V = (1/n) sum W_h^2 S_h^2 / a_h - (1/N) sum W_h S_h^2,
# so the n that gives V = (e / z)^2 is
#   n = sum(W_h^2 S_h^2 / a_h) / ((e / z)^2 + sum(W_h S_h^2) / N).

strata <- do.call(rbind, lapply(split(population, population$region),
  function(d) data.frame(
    region = d$region[1],
    N_h    = nrow(d),
    S_h    = sd(d$winnings_usd),
    Sp_h   = sqrt(nrow(d) / (nrow(d) - 1) * mean(d$below_ev) * (1 - mean(d$below_ev)))
  )))
strata$W_h <- strata$N_h / N

n_alloc <- function(e, S_h, a_h) {
  sum(strata$W_h^2 * S_h^2 / a_h) / ((e / z)^2 + sum(strata$W_h * S_h^2) / N)
}
n_srs <- function(e, S) {
  n0 <- (z * S / e)^2
  n0 / (1 + n0 / N)
}

# Allocation shares. The design uses Neyman allocation for the mean, so the
# proportion is shown under that same allocation.
shares <- list(
  "Proportional"           = strata$W_h,
  "Equal"                  = rep(1 / nrow(strata), nrow(strata)),
  "Neyman for mean (used)" = strata$N_h * strata$S_h / sum(strata$N_h * strata$S_h)
)

e_mean <- rel_e * mean(population$winnings_usd)
e_prop <- rel_e * mean(population$below_ev)
S_p    <- sqrt(N / (N - 1) * mean(population$below_ev) * (1 - mean(population$below_ev)))

alloc <- rbind(
  data.frame(parameter = "Mean winnings", design = "SRSWOR",
             n = n_srs(e_mean, sd(population$winnings_usd))),
  data.frame(parameter = "Mean winnings", design = names(shares),
             n = sapply(shares, n_alloc, e = e_mean, S_h = strata$S_h)),
  data.frame(parameter = "Proportion below opening EV", design = "SRSWOR",
             n = n_srs(e_prop, S_p)),
  data.frame(parameter = "Proportion below opening EV", design = names(shares),
             n = sapply(shares, n_alloc, e = e_prop, S_h = strata$Sp_h))
)
alloc$n      <- ceiling(alloc$n)
alloc$design <- factor(alloc$design, levels = c("SRSWOR", names(shares)))
alloc$used   <- alloc$design == "Neyman for mean (used)"
print(alloc)

p3 <- ggplot(alloc, aes(design, n, fill = used)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4, size = 3.5) +
  facet_wrap(~parameter, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "steelblue", `FALSE` = "grey65"), guide = "none") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Respondents needed for a 5% half-width",
       title = "Neyman allocation needs the fewest respondents for the mean",
       subtitle = "The same allocation is the least efficient choice for the proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path("figures", "B_allocation_comparison.png"), p3, width = 9, height = 4.5)

# ---- 4. Realised half-widths against the 5% target ------------------------
# Half-width of each repetition's 95% CI as a percentage of the true mean.

hw <- reps[reps$parameter == "Mean winnings", ]
hw$half_width_pct <- 100 * z * hw$se / hw$true

met <- aggregate(half_width_pct ~ design, hw, function(x) mean(x <= 100 * rel_e))
met$label <- paste0(round(100 * met$half_width_pct), "% of samples met the target")

p4 <- ggplot(hw, aes(half_width_pct, fill = design)) +
  geom_histogram(binwidth = 0.05, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 100 * rel_e, linetype = "dashed") +
  geom_text(data = met, aes(x = 100 * rel_e, y = Inf, label = label, colour = design),
            hjust = -0.05, vjust = c(2, 4), size = 3.5, inherit.aes = FALSE,
            show.legend = FALSE) +
  labs(x = "Half-width of 95% CI (% of true mean)", y = "Number of repetitions",
       fill = "Design",
       title = "Realised precision of the mean over 1,000 repetitions",
       subtitle = "Dashed line = 5% target") +
  theme_minimal()
ggsave(file.path("figures", "E_half_width_distribution.png"), p4, width = 7, height = 4.5)
