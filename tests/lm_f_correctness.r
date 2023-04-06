library(matrixTests)

#--- functions -----------------------------------------------------------------

base_lm <- function(mat, mod1, mod0) {
  stopifnot(ncol(mat) == nrow(mod1))
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)
  n <- rsq1 <- rsq0 <- dfm <- dfr <- f <- p <- numeric(nrow(mat))
  betas <- matrix(nrow=nrow(mat), ncol=ncol(mod1))
  colnames(betas) <- paste0("beta.", 1:ncol(betas)-1)

  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- mat[i,!bad]
    lm0 <- lm(vec ~ 0 + mod0[!bad,])
    lm1 <- lm(vec ~ 0 + mod1[!bad,])
    res <- anova(lm0, lm1)

    n[i]      <- length(vec)
    betas[i,] <- coefficients(lm1)
    rsq0[i]   <- 1 - (res$RSS[1] / sum((X[i,] - mean(X[i,]))^2))
    rsq1[i]   <- 1 - (res$RSS[2] / sum((X[i,] - mean(X[i,]))^2))
    dfm[i]    <- res$Df[2]
    dfr[i]    <- res$Res.Df[2]
    f[i]      <- res$F[2]
    p[i]      <- res[2,6]
  }

  data.frame(obs=n, betas, rsquared.model=rsq1, rsquared.null=rsq0,
             df.model=dfm, df.residual=dfr, statistic=f, pvalue=p
             )
}


#--- montecarlo ----------------------------------------------------------------

# random 100 rows with various covariates
X  <- matrix(rnorm(1000), ncol=10)
m  <- model.matrix(~ rnorm(10) + sample(1:10) + runif(10, -1, 1) + sample(c(0,1), 10, replace=TRUE))
m0 <- m[,c(1,sample(2:ncol(m), 3))]
res1 <- base_lm(X, m, m0)
res2 <- row_lm_f(X, m, m0 )
stopifnot(all.equal(res1, res2))



#--- extreme numbers -----------------------------------------------------------

# NOTE: extreme number test is skipped because cosinor_cosinor implementation is not robust against them

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

