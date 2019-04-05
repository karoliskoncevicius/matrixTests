context("correctness of jarquebera")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

moments_jarquebera <- function(mat) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  n <- skew <- kurt <- st <- df <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    vec <- na.omit(mat[i,])
    res <- moments::jarque.test(x=as.numeric(vec))

    n[i]    <- length(vec)
    skew[i] <- moments::skewness(vec)
    kurt[i] <- moments::kurtosis(vec)
    st[i]   <- res$statistic
    df[i]   <- 2
    p[i]    <- res$p.value
  }

  data.frame(obs=n, skewness=skew, kurtosis=kurt, df=df, statistic=st, pvalue=p)
}

################################################################################
########################### TEST ON A RANDOM SAMPLE ############################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(10000), ncol=10)
  X[sample(length(X), nrow(X))] <- NA

  t1 <- moments_jarquebera(X)
  t2 <- row_jarquebera(X)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  vals <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
  t1 <- moments_jarquebera(vals)
  t2 <- row_jarquebera(vals)
  expect_equal(t1, t2)

  # small numbers
  vals <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0)
  t1 <- moments_jarquebera(vals)
  t2 <- row_jarquebera(vals)
  expect_equal(t1, t2)
})


test_that("minumum allowed sample sizes give equal results", {
  # two numbers
  x <- matrix(rnorm(6), ncol=2)
  t1 <- moments_jarquebera(x)
  t2 <- row_jarquebera(x)
  expect_equal(t1, t2)

  # three numbers with NAs
  x <- matrix(rnorm(9), ncol=3); x[,3] <- NA
  t1 <- moments_jarquebera(x)
  t2 <- row_jarquebera(x)
  expect_equal(t1, t2)
})


################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warnign when row has less than 2 available observations", {
  wrn <- 'row_jarquebera: 1 of the rows had less than 2 total observations\\.\nFirst occurrence at row 1'
  nacolumns <- c("skewness", "kurtosis", "statistic", "pvalue")

  # single observation
  expect_warning(res <- row_jarquebera(1), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs, 1)
  expect_equal(res$df, 2)

  # single observation with some NA values
  expect_warning(res <- row_jarquebera(c(0,NA,NA)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs, 1)
  expect_equal(res$df, 2)

  # zero observations
  expect_warning(res <- row_jarquebera(NA_integer_), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs, 0)
  expect_equal(res$df, 2)
})


test_that("warning when a row has all constant values", {
  wrn <- 'row_jarquebera: 1 of the rows were essentially constant\\.\nFirst occurrence at row 1'
  nacolumns <- c("skewness", "kurtosis", "statistic", "pvalue")

  # two equal observations
  expect_warning(res <- row_jarquebera(c(1,1)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs, 2)
  expect_equal(res$df, 2)

  # three observations with some NA values
  expect_warning(res <- row_jarquebera(c(0,0,0,NA,NA)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs, 3)
  expect_equal(res$df, 2)
})

