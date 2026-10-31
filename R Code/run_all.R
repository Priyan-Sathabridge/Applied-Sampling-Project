# Run from the project directory: Rscript run_all.R
# No automatic installation: fail with a clear message if packages are absent.
required <- c("survey", "ggplot2", "scales", "rmarkdown", "knitr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these R packages first: ", paste(missing, collapse = ", "))
if (!file.exists("01_simulate_population.R")) stop("Run from the submission project directory.")
scripts <- c("01_simulate_population.R", "02_sample_size.R", "03_draw_sample.R",
             "04_estimation.R", "05_simulation.R", "06_extra_figures.R")
for (script in scripts) {
  message("Running ", script)
  # Isolated environments prevent accidental dependence on previous scripts.
  sys.source(script, envir = new.env(parent = globalenv()))
}
source("validate.R", local = new.env(parent = globalenv()))
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
if (!"--no-render" %in% commandArgs(trailingOnly = TRUE))
  rmarkdown::render("report.Rmd", output_format = "pdf_document", quiet = TRUE)
