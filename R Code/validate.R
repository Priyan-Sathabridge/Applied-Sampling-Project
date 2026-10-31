# Submission integrity checks, independent of the report's narrative.
pop <- read.csv("data/population.csv")
s <- read.csv("data/sample.csv")
r <- s[s$responded == 1, ]
x <- read.csv("output/E_replicates.csv")
d <- read.csv("output/D_estimates.csv")
stopifnot(nrow(pop) == 10000, !anyDuplicated(pop$id), !anyDuplicated(s$id),
          all(s$id %in% pop$id), all(s$pi_i > 0 & s$pi_i <= 1),
          all(abs(s$w_design * s$pi_i - 1) < 1e-10),
          all(is.na(s$winnings_usd[s$responded == 0])),
          abs(sum(r$w_final) - nrow(pop)) < 1e-8,
          nrow(x) == 4000, all(table(x$design, x$parameter) == 1000),
          all(is.finite(x$se)), all(x$se > 0),
          all(x$lower < x$est & x$est < x$upper),
          all(abs(d$se^2 - d$variance) < 1e-6))
for (h in unique(s$region)) {
  stopifnot(abs(sum(r$w_final[r$region == h]) - sum(pop$region == h)) < 1e-8)
}
message("PASS: population, samples, probabilities, weights, intervals and 1,000 repetitions per design.")
