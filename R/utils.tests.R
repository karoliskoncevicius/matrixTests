# Obtain p-values and confidence intervals for t-tests
do_ttest <- function(mx, mu, stder, alt, df, conf) {
  res <- matrix(numeric(), nrow=length(mx), ncol=4)
  colnames(res) <- c("t", "p", "cl", "ch")

  df[df<=0] <- NA

  res[,1] <- (mx-mu)/stder

  inds <- alt=="less"
  if(any(inds)) {
    res[inds,2] <- stats::pt(res[inds,1], df[inds])
    res[inds,3] <- rep(-Inf, sum(inds))
    res[inds,4] <- res[inds,1] + stats::qt(conf[inds], df[inds])
  }

  inds <- alt=="greater"
  if(any(inds)) {
    res[inds,2] <- stats::pt(res[inds,1], df[inds], lower.tail=FALSE)
    res[inds,3] <- res[inds,1] - stats::qt(conf[inds], df[inds])
    res[inds,4] <- rep(Inf, sum(inds))
  }

  inds <- alt=="two.sided"
  if(any(inds)) {
    res[inds,2] <- 2 * stats::pt(-abs(res[inds,1]), df[inds])
    intrange    <- stats::qt(1 - (1-conf[inds])/2, df[inds])
    res[inds,3] <- res[inds,1] - intrange
    res[inds,4] <- res[inds,1] + intrange
  }

  res[,3:4] <- mu + res[,3:4] * stder

  res
}

# Obtain p-values and confidence intervals for Pearson correlation test
do_pearson <- function(r, df, alt, conf) {
  res <- matrix(numeric(), nrow=length(r), ncol=4)
  colnames(res) <- c("stat", "p", "cl", "ch")

  df[df<=0] <- NA

  res[,1] <- sqrt(df)*r / sqrt(1 - r^2)
  z <- atanh(r)
  sigma <- 1/sqrt(df-1)

  inds <- alt=="less"
  if(any(inds)) {
    res[inds,2] <- stats::pt(res[inds,1], df[inds])
    res[inds,3] <- rep(-Inf, sum(inds))
    res[inds,4] <- z[inds] + sigma[inds] * stats::qnorm(conf[inds])
  }

  inds <- alt=="greater"
  if(any(inds)) {
    res[inds,2] <- stats::pt(res[inds,1], df[inds], lower.tail=FALSE)
    res[inds,3] <- z[inds] - sigma[inds] * stats::qnorm(conf[inds])
    res[inds,4] <- rep(Inf, sum(inds))
  }

  inds <- alt=="two.sided"
  if(any(inds)) {
    res[inds,2] <- 2 * pmin(stats::pt(res[inds,1], df[inds]), stats::pt(res[inds,1], df[inds], lower.tail=FALSE))
    res[inds,3] <- z[inds] + -sigma[inds] * stats::qnorm((1 + conf[inds])/2)
    res[inds,4] <- z[inds] + sigma[inds] * stats::qnorm((1 + conf[inds])/2)
  }

  res[,3] <- tanh(res[,3])
  res[,4] <- tanh(res[,4])

  res
}
