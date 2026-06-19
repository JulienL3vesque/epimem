# Regenerate R `mem` memstability reference for tests/test_equivalence.py.
#
#   Rscript tests/reference/generate_stability_reference.R
#
# Writes stability_reference.csv (one row per expanding-window fit) next to this script.
suppressMessages(library(mem))
data(flucyl)

args   <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
out    <- if (length(script)) dirname(normalizePath(script)) else "."

res <- memstability(flucyl)
sd <- res$stability.data
sd$count <- as.integer(rownames(sd))
sd <- sd[, c("count", setdiff(colnames(sd), "count"))]
write.csv(sd, file.path(out, "stability_reference.csv"), row.names = FALSE)

cat("wrote stability_reference.csv to", out, "\n")
options(digits = 10)
print(res$stability.data)
