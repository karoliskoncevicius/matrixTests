context("correctness of cosinor")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_cosinor <- function(mat, time, per) {
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

  data.frame(obs=n, mesor=m, amplitude=amp, acrophase=acr,
             df.model=dfm, df.residual=dfr, rsquared=rsq, statistic=f,
             pvalue=p, period=per
             )
}

################################################################################
########################### TEST ON A RANDOM SAMPLE ############################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  # equally spaced time points
  set.seed(14)
  X <- matrix(rnorm(10000), ncol=10)
  X[sample(length(X), nrow(X))] <- NA
  T <- 1:10

  t1 <- base_cosinor(X, T, 10)
  t2 <- row_cosinor(X, T, 10)

  expect_equal(t1, t2)

  # random time points
  set.seed(14)
  X <- matrix(rnorm(10000), ncol=10)
  X[sample(length(X), nrow(X))] <- NA
  per <- runif(1, 0 ,10)
  T   <- runif(10, 0 ,10)

  t1 <- base_cosinor(X, T, per)
  t2 <- row_cosinor(X, T, per)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("corner acrophases give correct results", {
  wave <- sinpi(2*1:24/24)
  expect_equal(suppressWarnings(row_cosinor(wave, 1:24, 24)$acrophase), 6)
  wave <- sinpi(2*13:36/24)
  expect_equal(suppressWarnings(row_cosinor(wave, 1:24, 24)$acrophase), 18)
  wave <- cospi(2*1:24/24)
  expect_equal(suppressWarnings(row_cosinor(wave, 1:24, 24)$acrophase), 0)
  wave <- cospi(2*13:36/24)
  expect_equal(suppressWarnings(row_cosinor(wave, 1:24, 24)$acrophase), 12)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning is shown when columns are removed because of NA times", {
  wrn <- '3 columns dropped due to missing time information'

  # 3 NAs
  expect_warning(res <- row_cosinor(1:11, c(1,2,3,4,NA,NA,NA,4,3,2,1)), wrn, all=TRUE)
  expect_equal(res$obs, 8)
})

test_that("warning when a infinite values are removed", {
  wrn <- 'row_cosinor: 1 of the rows had infinite observations that were removed\\.\nFirst occurrence at row 1'

  # -Inf and Inf among observations
  expect_warning(res <- row_cosinor(c(1,1,Inf,2,3,4), c(1,1,2,2,3,3)), wrn, all=TRUE)
  expect_equal(res$obs, 5)
})

test_that("warning when rows have exactly 3 complete observations", {
  wrn <- 'row_cosinor: 1 of the rows had exactly 3 complete observations: no p-values produced\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # standard case
  expect_warning(res <- row_cosinor(c(1,1,2), c(1,2,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_true(all(!is.na(res[,!colnames(res) %in% nacolumns])))

  # with NAs present
  expect_warning(res <- row_cosinor(c(1,1,4,NA), c(1,2,3,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_true(all(!is.na(res[,!colnames(res) %in% nacolumns])))
})

test_that("warning when rows have less than 3 complete observations", {
  wrn <- 'row_cosinor: 1 of the rows had less than 3 complete observations: no p-values produced, amplitude and acrophase will be unreliable\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # 1 observations
  expect_warning(res <- row_cosinor(1, 1), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # 2 observations
  expect_warning(res <- row_cosinor(c(1,2), c(1,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # with NAs present
  expect_warning(res <- row_cosinor(c(2,1,NA,NA), c(1,2,1,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
})

test_that("when less than 3 unique timepoints are available", {
  wrn <- 'row_cosinor: 1 of the rows had less than 3 unique timepoints within the specified period: amplitude and acrophase will be unreliable\\.\nFirst occurrence at row 1'

  # two distinct points, duplicated multiple times
  expect_warning(res <- row_cosinor(rnorm(4), c(1,2,1,2)), wrn, all=TRUE)

  # four points, but period is such that they are not unique
  expect_warning(res <- row_cosinor(rnorm(4), c(1,2,3,4), 2), wrn, all=TRUE)

  # with NAs present
  expect_warning(res <- row_cosinor(c(1,2,3,4,NA,NA), c(1,2,1,2,3,4)), wrn, all=TRUE)
})

test_that("when one of the rows has constant values", {
  wrn <- 'row_cosinor: 1 of the rows were essentially constant\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # only x
  expect_warning(res <- row_cosinor(c(1,1,1,1), c(1,2,3,4)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # with NAs present
  expect_warning(res <- row_cosinor(c(1,1,1,1,4), c(2,2,2,2,NA)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
})

test_that("warning when the fit is perfect", {
  wrn <- 'row_cosinor: 1 of the rows had essentially perfect fit\\.\nFirst occurrence at row 1'

  # perfect sine
  expect_warning(res <- row_cosinor(sin(2*pi*1:24/24), 1:24), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic))

  # only three distinct points present
  expect_warning(res <- row_cosinor(c(1:3,3), c(1:3,3), 24), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic))
})

