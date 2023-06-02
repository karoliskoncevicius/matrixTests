library(matrixTests)

#--- functions -----------------------------------------------------------------

moments_jarquebera <- function(mat) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  n <- skew <- kurt <- st <- df <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    vec <- na.omit(mat[i,])

    res <- moments::jarque.test(x=as.numeric(vec))

    # if p-value is NA turn df to NA as well
    df[i] <- ifelse(is.na(res$p.value), NA, 2)

    n[i]    <- length(vec)
    skew[i] <- moments::skewness(vec)
    kurt[i] <- moments::kurtosis(vec)
    st[i]   <- res$statistic
    df[i]   <- df
    p[i]    <- res$p.value
  }

  data.frame(obs=n, skewness=skew, kurtosis=kurt, df=df, statistic=st, pvalue=p)
}


#--- montecarlo ----------------------------------------------------------------

# 5 observations
x <- matrix(rnorm(5000), ncol=5)
res1 <- moments_jarquebera(x)
res2 <- row_jarquebera(x)
stopifnot(all.equal(res1, res2))

# 100 observations
x <- matrix(rnorm(100000), ncol=100)
res1 <- moments_jarquebera(x)
res2 <- row_jarquebera(x)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
       100000000000003, 100000000000002, 100000000000003, 100000000000000
       )
res1 <- moments_jarquebera(x)
res2 <- row_jarquebera(x)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
       1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
       )
res1 <- moments_jarquebera(x)
res2 <- row_jarquebera(x)
stopifnot(all.equal(res1, res2))

# NOTE: turned-off because of precission errors on architectures without long doubles
# large sample
# x <- rnorm(10^6)
# res1 <- moments_jarquebera(x)
# res2 <- row_jarquebera(x)
# stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# two numbers
x <- rnorm(2)
res1 <- moments_jarquebera(x)
res2 <- row_jarquebera(x)
stopifnot(all.equal(res1, res2))

# three numbers with NAs
x <- c(NA, rnorm(2))
res1 <- moments_jarquebera(x)
res2 <- row_jarquebera(x)
stopifnot(all.equal(res1, res2))

# four numbers
x <- rnorm(4)
res1 <- moments_jarquebera(x)
res2 <- row_jarquebera(x)
stopifnot(all.equal(res1, res2))
