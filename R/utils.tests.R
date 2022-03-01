# Obtain p-values and confidence intervals for t-tests
do_ttest <- function(mx, mu, stder, alt, df, conf) {
  res <- matrix(numeric(), nrow=length(mx), ncol=4)
  colnames(res) <- c("t", "p", "cl", "ch")

  df[df<=0] <- NA

  res[,1] <- (mx-mu)/stder

  inds <- alt=="less"
  if(any(inds)) {
    res[inds,2] <- stats::pt(res[inds,1], df[inds])
    res[inds,3] <- rep.int(-Inf, sum(inds))
    res[inds,4] <- res[inds,1] + stats::qt(conf[inds], df[inds])
  }

  inds <- alt=="greater"
  if(any(inds)) {
    res[inds,2] <- stats::pt(res[inds,1], df[inds], lower.tail=FALSE)
    res[inds,3] <- res[inds,1] - stats::qt(conf[inds], df[inds])
    res[inds,4] <- rep.int(Inf, sum(inds))
  }

  inds <- alt=="two.sided"
  if(any(inds)) {
    res[inds,2] <- 2 * stats::pt(-abs(res[inds,1]), df[inds])
    intrange    <- stats::qt(1 - (1-conf[inds])*0.5, df[inds])
    res[inds,3] <- res[inds,1] - intrange
    res[inds,4] <- res[inds,1] + intrange
  }

  res[,3:4] <- mu + res[,3:4] * stder

  res
}

# Obtain p-values and confidence intervals for f-tests
do_ftest <- function(est, rat, alt, df1, df2, conf) {
  res <- matrix(numeric(), nrow=length(est), ncol=4)
  colnames(res) <- c("f", "p", "cl", "ch")

  df1[df1 <= 0] <- NA
  df2[df2 <= 0] <- NA

  res[,1] <- est/rat

  inds <- alt=="less"
  if(any(inds)) {
    res[inds,2] <- stats::pf(res[inds,1], df1[inds], df2[inds])
    res[inds,3] <- 0
    res[inds,4] <- est[inds] / stats::qf(1 - conf[inds], df1[inds], df2[inds])
  }

  inds <- alt=="greater"
  if(any(inds)) {
    res[inds,2] <- 1 - stats::pf(res[inds,1], df1[inds], df2[inds])
    res[inds,3] <- est[inds] / stats::qf(conf[inds], df1[inds], df2[inds])
    res[inds,4] <- Inf
  }

  inds <- alt=="two.sided"
  if(any(inds)) {
    pval <- stats::pf(res[inds,1], df1[inds], df2[inds])
    beta <- (1 - conf[inds]) * 0.5
    res[inds,2] <- 2 * pmin(pval, 1 - pval)
    res[inds,3] <- est[inds] / stats::qf(1-beta, df1[inds], df2[inds])
    res[inds,4] <- est[inds] / stats::qf(beta, df1[inds], df2[inds])
  }

  res
}


do_wilcox_1_exact <- function(stat, n, alt) {
  res <- rep.int(NA_real_, length(stat))

  case <- stat > (n * (n+1)*0.25)

  # if n == 0 then we leave p-value as NA, hence: n!=0
  inds <- alt=="two.sided" & n!=0 & case
  if(any(inds)) {
    res[inds] <- stats::psignrank(stat[inds]-1, n[inds], lower.tail=FALSE)
    res[inds] <- pmin(2*res[inds], 1)
  }

  inds <- alt=="two.sided" & n!=0 & !case
  if(any(inds)) {
    res[inds] <- stats::psignrank(stat[inds], n[inds])
    res[inds] <- pmin(2*res[inds], 1)
  }

  inds <- alt=="less" & n!=0
  if(any(inds)) {
    res[inds] <- stats::psignrank(stat[inds], n[inds])
  }

  inds <- alt=="greater" & n!=0
  if(any(inds)) {
    res[inds] <- stats::psignrank(stat[inds]-1, n[inds], lower.tail=FALSE)
  }

  res
}

do_wilcox_1_approx <- function(stat, n, alt, nties, correct) {
  res <- rep.int(NA_real_, length(stat))

  z <- stat - n * (n+1)*0.25
  correction <- rep.int(0, length(stat))
  correction[correct & alt=="two.sided"] <- sign(z[correct & alt=="two.sided"]) * 0.5
  correction[correct & alt=="greater"]   <- 0.5
  correction[correct & alt=="less"   ]   <- -0.5
  z <- z - correction

  sigma <- sqrt(n * (n+1) * (2*n + 1)/24 - rowSums(nties^3 - nties, na.rm=TRUE)/48)
  z <- z/sigma


  inds <- alt=="two.sided"
  if(any(inds)) {
    res[inds] <- 2 * pmin(stats::pnorm(z[inds]), stats::pnorm(z[inds], lower.tail=FALSE))
  }

  inds <- alt=="greater"
  if(any(inds)) {
    res[inds] <- stats::pnorm(z[inds], lower.tail=FALSE)
  }

  inds <- alt=="less"
  if(any(inds)) {
    res[inds] <- stats::pnorm(z[inds])
  }

  res
}

do_wilcox_2_exact <- function(stat, nx, ny, alt) {
  res <- rep.int(NA_real_, length(stat))

  case <- stat > (nx*ny*0.5)


  # if nx == 0 or ny == 0 then we leave p-value as NA, hence: nx!=0 & ny!=0
  inds <- alt=="two.sided" & nx!=0 & ny!=0 & case
  if(any(inds)) {
    res[inds] <- stats::pwilcox(stat[inds]-1, nx[inds], ny[inds], lower.tail=FALSE)
    res[inds] <- pmin(2*res[inds], 1)
  }

  inds <- alt=="two.sided" & nx!=0 & ny!=0 & !case
  if(any(inds)) {
    res[inds] <- stats::pwilcox(stat[inds], nx[inds], ny[inds])
    res[inds] <- pmin(2*res[inds], 1)
  }

  inds <- alt=="greater" & nx!=0 & ny!=0
  if(any(inds)) {
    res[inds] <- stats::pwilcox(stat[inds]-1, nx[inds], ny[inds], lower.tail=FALSE)
  }

  inds <- alt=="less" & nx!=0 & ny!=0
  if(any(inds)) {
    res[inds] <- stats::pwilcox(stat[inds], nx[inds], ny[inds])
  }

  res
}

do_wilcox_2_approx <- function(stat, nx, ny, alt, nties, correct) {
  res <- rep.int(NA_real_, length(stat))

  z <- stat - nx*ny*0.5
  correction <- rep.int(0, length(stat))
  correction[correct & alt=="two.sided"] <- sign(z[correct & alt=="two.sided"]) * 0.5
  correction[correct & alt=="greater"]   <- 0.5
  correction[correct & alt=="less"   ]   <- -0.5
  z <- z - correction

  sigma <- sqrt((nx*ny/12) * ((nx+ny+1) - rowSums(nties^3 - nties, na.rm=TRUE) / ((nx+ny) * (nx+ny-1))))
  z <- z/sigma


  inds <- alt=="two.sided"
  if(any(inds)) {
    res[inds] <- 2 * pmin(stats::pnorm(z[inds]), stats::pnorm(z[inds], lower.tail=FALSE))
  }

  inds <- alt=="greater"
  if(any(inds)) {
    res[inds] <- stats::pnorm(z[inds], lower.tail=FALSE)
  }

  inds <- alt=="less"
  if(any(inds)) {
    res[inds] <- stats::pnorm(z[inds])
  }

  res
}

# Obtain p-values and confidence intervals for Pearson correlation test
do_pearson <- function(r, df, alt, conf) {
  res <- matrix(numeric(), nrow=length(r), ncol=4)
  colnames(res) <- c("stat", "p", "cl", "ch")

  df[df<=0] <- NA

  res[,1] <- sqrt(df)*r / sqrt(1 - r*r)
  z <- atanh(r)
  sigma <- 1/sqrt(df-1)

  inds <- alt=="less"
  if(any(inds)) {
    res[inds,2] <- stats::pt(res[inds,1], df[inds])
    res[inds,3] <- rep.int(-Inf, sum(inds))
    res[inds,4] <- z[inds] + sigma[inds] * stats::qnorm(conf[inds])
  }

  inds <- alt=="greater"
  if(any(inds)) {
    res[inds,2] <- stats::pt(res[inds,1], df[inds], lower.tail=FALSE)
    res[inds,3] <- z[inds] - sigma[inds] * stats::qnorm(conf[inds])
    res[inds,4] <- rep.int(Inf, sum(inds))
  }

  inds <- alt=="two.sided"
  if(any(inds)) {
    res[inds,2] <- 2 * pmin(stats::pt(res[inds,1], df[inds]), stats::pt(res[inds,1], df[inds], lower.tail=FALSE))
    res[inds,3] <- z[inds] + -sigma[inds] * stats::qnorm((1 + conf[inds])*0.5)
    res[inds,4] <- z[inds] + sigma[inds] * stats::qnorm((1 + conf[inds])*0.5)
  }

  res[,3] <- tanh(res[,3])
  res[,4] <- tanh(res[,4])

  res
}


# Obtain statistics for linear regression
do_regression <- function(Y, X) {
  betas <- matrix(NA, nrow=ncol(X), ncol=nrow(Y))
  dfmod <- numeric(nrow(Y))
  dfres <- numeric(nrow(Y))
  sstot <- numeric(nrow(Y))
  ssres <- numeric(nrow(Y))
  rsqs  <- numeric(nrow(Y))
  fs    <- numeric(nrow(Y))

  nainds <- is.na(Y)
  groups <- rowSums(nainds)
  if(any(groups != 0)) {
    categ <- lapply(split(nainds[groups!=0,,drop=FALSE], row(nainds[groups!=0,,drop=FALSE])), which)
    groups[groups!=0] <- match(categ, unique(categ))
  }
  for(g in unique(groups)) {
    rowinds <- groups == g
    colinds <- !nainds[match(g, groups),]
    y   <- Y[rowinds,colinds,drop=FALSE]
    res <- stats::.lm.fit(X[colinds,,drop=FALSE], t(y))

    betas[,rowinds] <- res$coefficients
    betas[,rowinds][abs(betas[,rowinds]) < .Machine$double.eps] <- 0

    dfres[rowinds] <- sum(colinds) - res$rank
    dfmod[rowinds] <- res$rank - 1

    sstot[rowinds] <- rowSums((y - rowMeans(y, na.rm=TRUE))^2, na.rm=TRUE)
    sstot[rowinds][sstot[rowinds] < .Machine$double.eps] <- 0
    ssres[rowinds] <- colSums(res$residuals^2, na.rm=TRUE)
    ssres[rowinds][ssres[rowinds] < .Machine$double.eps] <- 0
    isequal <- abs(ssres[rowinds]-sstot[rowinds]) < .Machine$double.eps^0.5
    ssres[rowinds][isequal] <- sstot[rowinds][isequal]
    rsqs[rowinds] <- 1 - (ssres[rowinds]/sstot[rowinds])

    ssmod <- sstot[rowinds] - ssres[rowinds]
    msres <- ssres[rowinds] / dfres[rowinds]
    msmod <- ssmod / dfmod[rowinds]
    fs[rowinds] <- msmod / msres
  }

  ps <- stats::pf(fs, dfmod, dfres, lower.tail=FALSE)

  list(betas=betas, stats=data.frame(dfmod=dfmod, dfres=dfres, sstot=sstot, ssres=ssres, rsq=rsqs, f=fs, p=ps))
}

