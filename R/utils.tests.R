# Function used to obtain p-values and confidence intervals for t-tests
do_ttest <- function(mx, mu, stder, alt, df, conf) {
  res <- matrix(numeric(), nrow=length(mx), ncol=4)
  colnames(res) <- c("t", "p", "cl", "ch")

  res[,1] <- (mx-mu)/stder

  inds <- alt=="less"
  if(any(inds)) {
    res[inds,2]  <- pt(res[inds,1], df[inds])
    res[inds,3]  <- rep(-Inf, sum(inds))
    res[inds,4]  <- res[inds,1] + qt(conf[inds], df[inds])
  }

  inds <- alt=="greater"
  if(any(inds)) {
    res[inds,2]  <- pt(res[inds,1], df[inds], lower.tail=FALSE)
    res[inds,3]  <- res[inds,1] - qt(conf[inds], df[inds])
    res[inds,4]  <- rep(Inf, sum(inds))
  }

  inds <- alt=="two.sided"
  if(any(inds)) {
    res[inds,2]  <- 2 * pt(-abs(res[inds,1]), df[inds])
    intrange     <- qt(1 - (1-conf[inds])/2, df[inds])
    res[inds,3]  <- res[inds,1] - intrange
    res[inds,4]  <- res[inds,1] + intrange
  }

  res[,3:4] <- mu + res[,3:4] * stder

  res
}

