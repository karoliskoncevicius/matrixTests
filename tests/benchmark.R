
dims <- c(100, 1000, 10000, 100000, 1000000, 10000000, 100000000)
to <- te <- tw <- tp <- oe <- bt <- numeric(length(dims))

for(i in 1:length(dims)) {
  d <- dims[i]
  X <- matrix(rnorm(d), ncol=10)
  Y <- matrix(rnorm(d), ncol=10)
  g <- rep(c("A","B"), each=5)
  to[i] <- system.time(ttest_onegroup(X))[3]
  te[i] <- system.time(ttest_equalvar(X, Y))[3]
  tw[i] <- system.time(ttest_welch(X, Y))[3]
  tp[i] <- system.time(ttest_paired(X, Y))[3]
  oe[i] <- system.time(oneway_equalvar(X, g))[3]
  bt[i] <- system.time(bartlett(X, g))[3]
  print(i)
}

