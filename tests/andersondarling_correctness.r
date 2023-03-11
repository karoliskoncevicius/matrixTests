library(matrixTests)

#--- functions -----------------------------------------------------------------

nortest_andersondarling <- function(mat) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  n <- st <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    vec <- na.omit(mat[i,])

    res <- nortest::ad.test(x=as.numeric(vec))

    n[i]    <- length(vec)
    st[i]   <- res$statistic
    p[i]    <- res$p.value
  }

  data.frame(obs=n, statistic=st, pvalue=p)
}


#--- montecarlo ----------------------------------------------------------------

# 10 observations
x <- matrix(rnorm(10000), ncol=10)
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))

# 100 observations
x <- matrix(rnorm(100000), ncol=100)
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))

# 20 observations - highly non normal
x <- cbind(matrix(runif(10000), ncol=10), matrix(runif(10000, 10, 11), ncol=10))
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
       100000000000003, 100000000000002, 100000000000003, 100000000000000
       )
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
       1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
       )
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))

# extremely non normal
x <- 1:10000
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))

# large sample
x <- rnorm(10^6)
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# eight numbers
x <- rnorm(8)
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))

# eight numbers with NAs
x <- c(NA, rnorm(8))
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))

# nine numbers
x <- rnorm(9)
res1 <- nortest_andersondarling(x)
res2 <- row_andersondarling(x)
stopifnot(all.equal(res1, res2))

