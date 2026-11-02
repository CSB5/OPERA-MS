# Simulated species for opera-ms-pipeline::cluster_strain (Q17): does the
# strain assignment change when the modes come from R (with its FFT noise
# modes) instead of the Rust density? One case per block:
#   case <TAB> dispersion
#   contig <TAB> name <TAB> length <TAB> windows (space-separated)   (repeated)
#   rmodes <TAB> tokens as cluster_strain.pl splits R's output
# Usage: Rscript strain_cases.R <seed> <number of cases> [outlier factor, default 200] > cases.tsv
args <- commandArgs(TRUE)
set.seed(as.integer(args[1]))
n_cases <- as.integer(args[2])
window_size <- 340
outlier <- if (length(args) > 2) as.numeric(args[3]) else 200
for (case in seq_len(n_cases)) {
  dispersion <- signif(exp(runif(1, log(5), log(200))), 8)
  n_strains <- sample(1:3, 1)
  base <- exp(runif(1, log(5), log(300)))
  covs <- base * c(1, runif(2, 0.15, 0.8))[1:n_strains]
  names <- c(); lengths <- c(); windows <- list()
  k <- 0
  for (s in seq_len(n_strains)) {
    n_contigs <- sample(20:300, 1)
    for (i in seq_len(n_contigs)) {
      k <- k + 1
      nw <- sample(2:300, 1)
      size <- exp(runif(1, log(3), log(300)))
      mu <- covs[s]
      if (runif(1) < 0.05) mu <- mu * sample(2:20, 1)  # repeat contig
      w <- rnbinom(nw, size = size, mu = mu)
      if (runif(1) < 0.02) w[sample(nw, 1)] <- rnbinom(1, size = 2, mu = mu * outlier)  # outlier window
      names[k] <- sprintf("c%05d", k); lengths[k] <- nw * window_size + sample(0:339, 1); windows[[k]] <- w
    }
  }
  # cluster_strain.pl: windows by decreasing int(mean), ties by name, mean <= 100000
  means <- sapply(windows, function(w) trunc(mean(w)))
  o <- order(-means, names, method = "radix")
  o <- o[means[o] <= 100000]
  values <- unlist(windows[o])
  dens <- density(values)
  z <- embed(dens$y, 11); s <- 5
  ind <- apply(z, 1, which.max)
  result <- c(rep(FALSE, s), ind == (1 + s))
  result <- result[1:(length(result) - s)]
  out <- capture.output(cat(capture.output(print(length(values))), "\n"), print(dens$x[result]))
  tokens <- unlist(strsplit(paste(out, collapse = " "), "[[:space:]]+"))
  tokens <- tokens[tokens != ""]
  tokens <- tokens[-1]  # the Perl drops the first token ("[1]")
  cat("case\t", dispersion, "\n", sep = "")
  for (i in seq_along(names)) cat("contig\t", names[i], "\t", lengths[i], "\t", paste(windows[[i]], collapse = " "), "\n", sep = "")
  cat("rmodes\t", paste(tokens, collapse = " "), "\n", sep = "")
}
