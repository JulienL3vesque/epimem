# Regenerate R `mem` roc.analysis reference for tests/test_equivalence.py.
#
#   Rscript tests/reference/generate_roc_reference.R
#
# Writes roc_reference.csv (the sweep table) + roc_optimum.csv (winner per criterion)
# next to this script. Uses a coarse grid - the full seq(1,5,0.1) sweep is slow.
suppressMessages(library(mem))
data(flucyl)

args   <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
out    <- if (length(script)) dirname(normalizePath(script)) else "."

res <- roc.analysis(flucyl, i.param.values = seq(2, 3, 0.5), i.min.seasons = 6)

write.csv(res$roc.data, file.path(out, "roc_reference.csv"), row.names = FALSE)
opt <- res$optimum
write.csv(data.frame(criterion = names(opt), value = as.numeric(opt[1, ])),
          file.path(out, "roc_optimum.csv"), row.names = FALSE)

cat("wrote roc_reference.csv + roc_optimum.csv to", out, "\n")
print(res$optimum)
print(res$roc.data)
