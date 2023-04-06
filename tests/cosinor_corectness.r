library(matrixTests)

#--- functions -----------------------------------------------------------------

cosinor_cosinor <- function(mat, time, per) {
  stopifnot(ncol(mat) == length(time))
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)
  n <- m <- amp <- acr <- dfm <- dfr <- rsq <- f <- p <- numeric(nrow(mat))

  for(i in 1:nrow(mat)) {
    dat <- data.frame(x=mat[i,], t=time)
    dat <- dat[!is.na(dat$x) & !is.na(dat$t),]
    res <- cosinor::cosinor.lm(x ~ time(t), data=dat, period=per)

    n[i]   <- nrow(dat)
    m[i]   <- res$coefficients[1]
    amp[i] <- res$coefficients[2]
    acr[i] <- abs(cosinor2::correct.acrophase(res) / pi * (per/2))
    dfm[i] <- cosinor2::cosinor.detect(res)[2]
    dfr[i] <- cosinor2::cosinor.detect(res)[3]
    rsq[i] <- cosinor2::cosinor.PR(res)[[2]]
    f[i]   <- cosinor2::cosinor.detect(res)[1]
    p[i]   <- cosinor2::cosinor.detect(res)[4]
  }

  data.frame(obs=n, mesor=m, amplitude=amp, acrophase=acr, rsquared=rsq,
             df.model=dfm, df.residual=dfr, statistic=f,
             pvalue=p, period=per
             )
}


#--- montecarlo ----------------------------------------------------------------

# equally spaced 10 time points, period = 10
X <- matrix(rnorm(10000), ncol=10)
t <- -5:4
res1 <- cosinor_cosinor(X, t, 10)
res2 <- row_cosinor(X, t, 10)
stopifnot(all.equal(res1, res2))

# equally spaced 5 time points, period = 2.2
X <- matrix(rnorm(5000), ncol=5)
t <- 1:5
res1 <- cosinor_cosinor(X, t, 2.2)
res2 <- row_cosinor(X, t, 2.2)
stopifnot(all.equal(res1, res2))

# random 10 time points, random period
X <- matrix(rnorm(10000), ncol=10)
t   <- runif(10, 1, 10)
per <- runif(1, 1, 10)
res1 <- cosinor_cosinor(X, t, per)
res2 <- row_cosinor(X, t, per)
stopifnot(all.equal(res1, res2))

# random 50 time points, random period
X <- matrix(rnorm(50000), ncol=50)
t   <- runif(50, 1, 5)
per <- runif(1, 1, 10)
res1 <- cosinor_cosinor(X, t, per)
res2 <- row_cosinor(X, t, per)
stopifnot(all.equal(res1, res2))

# small period
X <- matrix(rnorm(10000), ncol=10)
t <- c(1:9, 9.3)
res1 <- cosinor_cosinor(X, t, 0.301)
res2 <- row_cosinor(X, t, 0.301)
stopifnot(all.equal(res1, res2))

# large period
X <- matrix(rnorm(10000), ncol=10)
t <- 1:10
res1 <- cosinor_cosinor(X, t, 300)
res2 <- row_cosinor(X, t, 300)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# NOTE: extreme number test is skipped because cosinor_cosinor implementation is not robust against them

# large sample
x <- rnorm(10^6)
t <- runif(10^6, 0, 24)
res1 <- cosinor_cosinor(x, t, 24)
res2 <- row_cosinor(x, t, 24)
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# NOTE: matrixTests can work with 4 observations but cosinor_cosinor requires 5

# five observations
x <- rnorm(5)
t <- runif(5, 1, 24)
res1 <- cosinor_cosinor(x, t, 24)
res2 <- row_cosinor(x, t, 24)
stopifnot(all.equal(res1, res2))

# three unique time points each repeated 3 times
x <- rnorm(9)
t <- rep(runif(3, 1, 24), each=3)
res1 <- cosinor_cosinor(x, t, 24)
res2 <- row_cosinor(x, t, 24)
stopifnot(all.equal(res1, res2))

# three unique time points, one repeated 2 times and others 1 time
x <- rnorm(5)
t <- rep(runif(3, 1, 24), c(3,1,1))
res1 <- cosinor_cosinor(x, t, 24)
res2 <- row_cosinor(x, t, 24)
stopifnot(all.equal(res1, res2))

# 3 unique time points, all having the same values repeated two times
# NOTE: this produces a perfect fit, so cosinor_cosinor values are adjusted to match
x <- rep(rnorm(3), each=2)
t <- rep(runif(3, 1, 24), each=2)
res1 <- cosinor_cosinor(x, t, 24)
res1$statistic <- Inf; res1$pvalue <- 0
res2 <- suppressWarnings(row_cosinor(x, t, 24))
stopifnot(all.equal(res1, res2))


#--- constant values -----------------------------------------------------------

# all values are constant except one
x <- c(1,1,1,1,rnorm(1))
t <- runif(5, 1, 24)
res1 <- cosinor_cosinor(x, t, 24)
res2 <- row_cosinor(x, t, 24)
stopifnot(all.equal(res1, res2))

# two different values
x <- c(1,1,1,2,2,2)
t <- runif(6, 1, 24)
res1 <- cosinor_cosinor(x, t, 24)
res2 <- row_cosinor(x, t, 24)
stopifnot(all.equal(res1, res2))


#--- edge period values --------------------------------------------------------

# smallest viable period - 2x + epsilon compared to the gap
x <- rnorm(24)
t <- 1:24
p <- 2.0001
res1 <- cosinor_cosinor(x, t, p)
res2 <- row_cosinor(x, t, p)
stopifnot(all.equal(res1, res2))

# smallest viable period when only one small gap exists
x <- rnorm(24)
t <- c(1:23, 23.5)
p <- 1.0001
res1 <- cosinor_cosinor(x, t, p)
res2 <- row_cosinor(x, t, p)
stopifnot(all.equal(res1, res2))


#--- edge acrophase values -----------------------------------------------------

x <- sinpi(2*1:24/24)
res <- suppressWarnings(row_cosinor(x, 1:24, 24))
stopifnot(all.equal(res$acrophase, 6))

x <- sinpi(2*13:36/24)
res <- suppressWarnings(row_cosinor(x, 1:24, 24))
stopifnot(all.equal(res$acrophase, 18))

x <- cospi(2*1:24/24)
res <- suppressWarnings(row_cosinor(x, 1:24, 24))
stopifnot(all.equal(res$acrophase, 0))

x <- cospi(2*13:36/24)
res <- suppressWarnings(row_cosinor(x, 1:24, 24))
stopifnot(all.equal(res$acrophase, 12))

