# Regenerate R `mem` reference for the companion functions (memintensity, memtrend).
#
#   Rscript tests/reference/generate_companions_reference.R
#
# Writes intensity_reference.csv + trend_reference.csv next to this script.
suppressMessages(library(mem))
data(flucyl)

args   <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
out    <- if (length(script)) dirname(normalizePath(script)) else "."

epi <- memmodel(flucyl)

# memintensity: four labelled cut-points.
ii <- memintensity(epi)
write.csv(data.frame(label = colnames(ii$intensity.thresholds),
                     value = as.numeric(ii$intensity.thresholds)),
          file.path(out, "intensity_reference.csv"), row.names = FALSE)

# memtrend: ascending / descending change thresholds.
tr <- memtrend(epi)
write.csv(data.frame(name = c("ascending", "descending"),
                     value = as.numeric(tr$trend.thresholds)),
          file.path(out, "trend_reference.csv"), row.names = FALSE)

cat("wrote intensity_reference.csv + trend_reference.csv to", out, "\n")
print(ii$intensity.thresholds)
print(tr$trend.thresholds)

# Fidelity check: does memmodel store the raw data or a gap-filled copy? (matters for memtrend)
cat("memmodel$data identical to raw flucyl:", identical(epi$data, as.matrix(flucyl)), "\n")
