# Function used to obtain p-values and confidence intervals for t-tests
do_ttest <- function(mx, mu, stder, alt, df, conf) {
  pvals <- numeric(length(mx))
  cints <- matrix(nrow=length(mx), ncol=2)

  t <- (mx - mu)/stder

  inds <- alt=="less"
  pvals[inds]  <- pt(t[inds], df[inds])
  cints[inds,] <- cbind(-Inf, t[inds] + qt(conf[inds], df[inds]))

  inds <- alt=="greater"
  pvals[inds]  <- pt(t[inds], df[inds], lower.tail=FALSE)
  cints[inds,] <- cbind(t[inds] - qt(conf[inds], df[inds]), Inf)

  inds <- alt=="two.sided"
  pvals[inds]  <- 2 * pt(-abs(t[inds]), df[inds])
  intrange     <- qt(1 - (1-conf[inds])/2, df[inds])
  cints[inds,] <- t[inds] + cbind(-intrange, intrange)

  cints <- mu + cints * stder

  list(t.statistic=t, p.value=pvals, ci.low=cints[,1], ci.high=cints[,2])
}

