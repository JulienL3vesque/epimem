# Regenerate R `mem` memevolution reference for tests/test_equivalence.py.
#
#   Rscript tests/reference/generate_evolution_reference.R
#
# Writes evolution_reference.csv (one row per season-by-season fit, for both validation methods)
# next to this script. evolution.data already carries a "number" column (= n.seasons); we prepend
# the method name and a 0-based row index so the test can line rows up with mem_evolution.
suppressMessages(library(mem))
data(flucyl)

args   <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
out    <- if (length(script)) dirname(normalizePath(script)) else "."

evo <- function(method) {
  ed <- memevolution(flucyl, i.evolution.method = method)$evolution.data
  data.frame(method = method, row = 0:(nrow(ed) - 1), ed, row.names = NULL)
}

both <- rbind(evo("sequential"), evo("cross"))
write.csv(both, file.path(out, "evolution_reference.csv"), row.names = FALSE)

cat("wrote evolution_reference.csv to", out, "\n")
options(digits = 10)
print(both[, c("method", "row", "number", "epidemic", "medium", "high", "veryhigh")])
