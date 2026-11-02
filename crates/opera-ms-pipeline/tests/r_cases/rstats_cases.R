# Golden cases for opera-ms-pipeline::rstats (Q17): R 3.5.1 used as a black
# box. Writes, one case per line:
#   modes <TAB> values (comma-separated) <TAB> R's printed modes (space-separated)
#         <TAB> R's density at each mode / its maximum (space-separated)
#   qnbinom <TAB> p <TAB> size <TAB> mu <TAB> R's result
# Usage: Rscript rstats_cases.R <seed> <number of density cases> <number of qnbinom cases> > cases.tsv
args <- commandArgs(TRUE)
set.seed(as.integer(args[1]))
n_density <- as.integer(args[2])
n_q <- as.integer(args[3])
modes <- function(values) {
  # compute_mode in cluster_strain.pl
  span <- 11
  dens <- density(values)
  series <- dens$y
  z <- embed(series, span)
  s <- span %/% 2
  ind <- apply(z, 1, which.max)
  v <- ind == (1 + s)
  result <- c(rep(FALSE, s), v)
  result <- result[1:(length(result) - s)]
  out <- capture.output(print(dens$x[result]))
  tokens <- unlist(strsplit(paste(out, collapse = " "), "[[:space:]]+"))
  tokens <- tokens[tokens != "" & !grepl("[", tokens, fixed = TRUE)]
  rel <- dens$y[result] / max(dens$y)
  c(paste(tokens, collapse = " "), paste(sprintf("%.3e", rel), collapse = " "))
}
for (i in seq_len(n_density)) {
  k <- sample(1:4, 1)
  values <- c()
  for (j in seq_len(k)) {
    mu <- exp(runif(1, log(2), log(500)))
    size <- exp(runif(1, log(0.5), log(500)))
    values <- c(values, rnbinom(sample(c(50, 500, 3000, 15000), 1), size = size, mu = mu))
  }
  if (runif(1) < 0.3) values <- c(values, rnbinom(sample(1:200, 1), size = 5, mu = exp(runif(1, log(1000), log(50000)))))
  values <- sample(values)
  m <- modes(values)
  cat("modes\t", paste(values, collapse = ","), "\t", m[1], "\t", m[2], "\n", sep = "")
}
for (i in seq_len(n_q)) {
  p <- sample(c(0.02, 0.98), 1)
  size <- signif(exp(runif(1, log(0.1), log(2000))), 15)
  mu <- signif(exp(runif(1, log(0.01), log(50000))), 15)
  cat("qnbinom\t", sprintf("%.17g", p), "\t", sprintf("%.17g", size), "\t", sprintf("%.17g", mu), "\t", qnbinom(p, size = size, mu = mu), "\n", sep = "")
}
