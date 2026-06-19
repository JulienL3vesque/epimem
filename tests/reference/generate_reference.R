# Regenerate the R `mem` reference numbers that tests/test_equivalence.py checks epimem against.
#
#   Rscript tests/reference/generate_reference.R
#
# Needs the mem package:  install.packages("mem")
# Writes flucyl.csv, flucyl_gappy.csv and reference_thresholds.csv next to this script.
suppressMessages(library(mem))
data(flucyl)

args   <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
out    <- if (length(script)) dirname(normalizePath(script)) else "."

named         <- c("epidemic_onset", "post_epidemic", "medium", "high", "very_high")
thresholds_of <- function(m) c(m$epidemic.thresholds, m$intensity.thresholds)

# Clean reference.
write.csv(flucyl, file.path(out, "flucyl.csv"), row.names = FALSE)
clean <- thresholds_of(memmodel(flucyl))

# Gappy reference: poke interior holes (incl. mid-season) so fill.missing actually runs.
g <- flucyl
g[16, 2] <- NA; g[17, 2] <- NA; g[15, 5] <- NA; g[18, 7] <- NA; g[19, 7] <- NA
write.csv(g, file.path(out, "flucyl_gappy.csv"), row.names = FALSE)
gappy <- thresholds_of(memmodel(g))

ref <- data.frame(
  case  = rep(c("clean", "gappy"), each = 5),
  name  = rep(named, 2),
  value = c(clean, gappy)
)
write.csv(ref, file.path(out, "reference_thresholds.csv"), row.names = FALSE)
cat("wrote reference files to", out, "\n")
print(ref)
