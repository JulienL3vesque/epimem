# Regenerate R `mem` memgoodness reference numbers for tests/test_equivalence.py.
#
#   Rscript tests/reference/generate_goodness_reference.R
#
# Writes goodness_reference.csv next to this script (default grid + a coarse grid).
suppressMessages(library(mem))
data(flucyl)

args   <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
out    <- if (length(script)) dirname(normalizePath(script)) else "."

names15 <- c("weeks", "non_missing_weeks", "true_positives", "false_positives",
             "true_negatives", "false_negatives", "sensitivity", "specificity",
             "positive_predictive_value", "negative_predictive_value",
             "positive_likelihood_ratio", "negative_likelihood_ratio",
             "percent_agreement", "matthews_corrcoef", "youden_index")

g_default <- memgoodness(flucyl, i.graph = FALSE)
g_coarse  <- memgoodness(flucyl, i.detection.values = seq(1, 5, 0.5), i.graph = FALSE)

ref <- rbind(
  data.frame(grid = "default", name = names15, value = as.numeric(g_default$results)),
  data.frame(grid = "coarse",  name = names15, value = as.numeric(g_coarse$results))
)
write.csv(ref, file.path(out, "goodness_reference.csv"), row.names = FALSE)
cat("wrote goodness_reference.csv to", out, "\n")
print(round(g_default$results, 6))
