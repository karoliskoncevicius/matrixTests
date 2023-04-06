library(matrixTests)

#--- functions -----------------------------------------------------------------

base_t_onesample <- function(mat, null=0, alternative="two.sided", conf=0.95) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)
  if(length(alternative)==1) alternative <- rep(alternative, nrow(mat))
  if(length(null)==1) null <- rep(null, nrow(mat))
  if(length(conf)==1) conf <- rep(conf, nrow(mat))

  mx <- vx <- nx <- tst <- p <- cl <- ch <- se <- df <- m0 <- cnf <- numeric(nrow(mat))
  al <- character(nrow(mat))
  for(i in 1:nrow(mat)) {
    vec <- na.omit(mat[i,])

    res <- t.test(vec, alternative=alternative[i], mu=null[i], conf.level=conf[i])

    # if p-value is NA turn stderr and df to NA as well
    if(is.na(res$p.value)) {
      res$stderr <- NA
      res$parameter <- NA
    }

    vx[i]  <- var(vec)
    nx[i]  <- length(vec)
    mx[i]  <- res$estimate
    tst[i] <- res$statistic
    p[i]   <- res$p.value
    cl[i]  <- res$conf.int[1]
    ch[i]  <- res$conf.int[2]
    se[i]  <- res$stderr
    df[i]  <- res$parameter
    m0[i]  <- res$null.value
    al[i]  <- res$alternative
    cnf[i] <- attr(res$conf.int, "conf.level")
  }

  data.frame(obs=nx, mean=mx, var=vx, stderr=se, df=df, statistic=tst,
             pvalue=p, conf.low=cl, conf.high=ch, mean.null=m0, alternative=al,
             conf.level=cnf, stringsAsFactors=FALSE
             )
}


#--- montecarlo ----------------------------------------------------------------

# 5 observations
x <- matrix(rnorm(5000), ncol=5)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
cfs  <- sample(seq(0, 1, length.out=nrow(x)))
res1 <- base_t_onesample(x, mus, alts, cfs)
res2 <- row_t_onesample(x, mus, alts, cfs)
stopifnot(all.equal(res1, res2))

# 20 observations
x <- matrix(rnorm(20000), ncol=20)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
cfs  <- sample(seq(0, 1, length.out=nrow(x)))
res1 <- base_t_onesample(x, mus, alts, cfs)
res2 <- row_t_onesample(x, mus, alts, cfs)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
res1 <- base_t_onesample(x)
res2 <- row_t_onesample(x)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0)
res1 <- base_t_onesample(x)
res2 <- row_t_onesample(x)
stopifnot(all.equal(res1, res2))

# large sample
x <- rnorm(10^6)
res1 <- base_t_onesample(x)
res2 <- row_t_onesample(x)
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# two numbers
x <- matrix(rnorm(6), ncol=2)
alt <- c("two.sided", "greater", "less")
res1 <- base_t_onesample(x, alternative=alt)
res2 <- row_t_onesample(x, alternative=alt)
stopifnot(all.equal(res1, res2))

# three numbers with NAs
x <- matrix(c(rnorm(6),NA,NA,NA), ncol=3)
alt <- c("two.sided", "greater", "less")
res1 <- base_t_onesample(x, alternative=alt)
res2 <- row_t_onesample(x, alternative=alt)
stopifnot(all.equal(res1, res2))


#--- parameter edge cases ------------------------------------------------------

# various corner cases with NAs
mus <- c(-Inf, -1, 0, 1, Inf)
alt <- c("l", "t", "g")
cfs <- c(0, 0.5, 1)
pars <- expand.grid(mus, alt, cfs, stringsAsFactors=FALSE)
x <- matrix(rnorm(10*nrow(pars)), ncol=10)
res1 <- base_t_onesample(x, pars[,1], pars[,2], pars[,3])
res2 <- row_t_onesample(x, pars[,1], pars[,2], pars[,3])
stopifnot(all.equal(res1, res2))

# NAs in confidence intervals
x <- matrix(rnorm(40), ncol=10)
cnf <- c(0.95, NA, 0.5, NA)
res1 <- base_t_onesample(x, conf=ifelse(is.na(cnf), 0.95, cnf))
res1[is.na(cnf), c("conf.level", "conf.low", "conf.high")] <- NA
res2 <- row_t_onesample(x, conf.level=cnf)
stopifnot(all.equal(res1, res2))

# null exactly equal to the mean
res1 <- base_t_onesample(c(1,2,3), null=2)
res2 <- row_t_onesample(c(1,2,3), null=2)
stopifnot(all.equal(res2$pvalue, 1))
stopifnot(all.equal(res1, res2))

