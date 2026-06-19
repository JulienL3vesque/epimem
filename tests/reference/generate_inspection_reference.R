# Regenerate R `mem` optimum.by.inspection reference for tests/test_equivalence.py.
#
#   Rscript tests/reference/generate_inspection_reference.R
#
# optimum.by.inspection is interactive in R (the analyst clicks each season's epidemic
# start/end). We reproduce its exact internal computation non-interactively with a fixed,
# reproducible inspection timing per season (each season's own 2.8 criterion optimum) and
# commit those timings, so the Python port can be checked against the same input.
suppressMessages(library(mem))
data(flucyl)
pdf(NULL)                       # swallow any graphics device

args   <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
out    <- if (length(script)) dirname(normalizePath(script)) else "."

i.data <- flucyl
pv <- seq(2.0, 3.0, 0.1)
anios <- ncol(i.data); n <- length(pv)
resultados.i <- array(dim = c(anios, 15, n))
timings <- array(dim = c(anios, 2))

for (i in 1:anios) {
  cur <- i.data[i]
  mt <- memtiming(cur, i.n.values = 5, i.method = 2, i.param = 2.8)
  t1 <- as.numeric(mt$optimum.map[4:5])          # inspection timing stand-in (start, end)
  timings[i, ] <- t1
  curva.map <- mem:::calcular.map(as.vector(as.matrix(cur)))
  for (j in 1:n) {
    t2 <- mem:::calcular.optimo(curva.map, 2, pv[j])$resultados[4:5]
    rj <- mem:::calcular.indicadores.2.timings(cur, t1, t2,
            i.timing.labels = c("insp", "x"), i.graph.file = FALSE)$indicadores
    resultados.i[i, , j] <- as.numeric(rj)
  }
}

resultado <- data.frame(apply(resultados.i, c(3, 2), sum, na.rm = TRUE))
resultado[7]  <- resultado[3] / (resultado[3] + resultado[6])
resultado[8]  <- resultado[5] / (resultado[5] + resultado[4])
resultado[9]  <- resultado[3] / (resultado[3] + resultado[4])
resultado[10] <- resultado[5] / (resultado[5] + resultado[6])
resultado[11] <- resultado[7] / (1 - resultado[8])
resultado[12] <- (1 - resultado[7]) / resultado[8]
resultado[13] <- (resultado[3] + resultado[5]) / (resultado[3] + resultado[4] + resultado[5] + resultado[6])
resultado[14] <- (resultado[3] * resultado[5] - resultado[4] * resultado[6]) /
  sqrt((resultado[3] + resultado[4]) * (resultado[3] + resultado[6]) *
       (resultado[5] + resultado[4]) * (resultado[5] + resultado[6]))
resultado[15] <- resultado[7] + resultado[8] - 1
resultado[sapply(resultado, is.nan)] <- NA

cols <- c("value", "weeks", "non_missing_weeks", "true_positives", "false_positives",
          "true_negatives", "false_negatives", "sensitivity", "specificity",
          "positive_predictive_value", "negative_predictive_value",
          "positive_likelihood_ratio", "negative_likelihood_ratio",
          "percent_agreement", "matthews_corrcoef", "youden_index")
resultados <- data.frame(value = pv, resultado)
colnames(resultados) <- cols
write.csv(resultados, file.path(out, "inspection_reference.csv"), row.names = FALSE)
write.csv(data.frame(season = 1:anios, start = timings[, 1], end = timings[, 2]),
          file.path(out, "inspection_timings.csv"), row.names = FALSE)

s <- resultados$sensitivity; sp <- resultados$specificity
r1 <- rank(-s, na.last = TRUE) + rank(-sp, na.last = TRUE)
r2 <- rank(-s * sp, na.last = TRUE)
r3 <- rank(-resultados$positive_likelihood_ratio, na.last = TRUE)
r4 <- rank(-resultados$negative_likelihood_ratio, na.last = TRUE)
qf <- abs(s - sp); qe <- 2 - s - sp; qs <- (1 - s)^2 + (1 - sp)^2
r5 <- rank(qf) + rank(qe) + rank(qs)
r6 <- rank(-resultados$percent_agreement, na.last = TRUE)
r7 <- rank(-resultados$matthews_corrcoef, na.last = TRUE)
r8 <- rank(-resultados$youden_index, na.last = TRUE)
opt <- data.frame(
  criterion = c("pos_likelihood", "neg_likelihood", "additive", "multiplicative",
                "mixed", "percent", "matthews", "youden"),
  value = c(pv[which.min(r3)], pv[which.min(r4)], pv[which.min(r1)], pv[which.min(r2)],
            pv[which.min(r5)], pv[which.min(r6)], pv[which.min(r7)], pv[which.min(r8)]))
write.csv(opt, file.path(out, "inspection_optimum.csv"), row.names = FALSE)

cat("wrote inspection_reference.csv + inspection_timings.csv + inspection_optimum.csv to", out, "\n")
print(opt)
print(round(resultados[, c("value", "sensitivity", "specificity", "matthews_corrcoef", "youden_index")], 6))
