context("mat_ttest_single")


test_that("t.test(x) gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(1000), ncol=10)
  X[sample(length(X), 100)] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  mus  <- seq(-1, 1, length.out=nrow(X))
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- mat_ttest_single(X, alts, mus, cfs)
  t2 <- Map(t.test, split(X, row(X)), alternative=alts, mu=mus, conf.level=cfs)

  expect_equal(t1$x.mean, unname(sapply(t2, "[[", "estimate")))
  expect_equal(t1$x.var, apply(X, 1, var, na.rm=TRUE))
  expect_equal(t1$x.obs, rowSums(!is.na(X)))
  expect_equal(t1$t.statistic, unname(sapply(t2, "[[", "statistic")))
  expect_equal(t1$p.value, unname(sapply(t2, "[[", "p.value")))
  expect_equal(t1$ci.low, unname(sapply(t2, "[[", "conf.int")[1,]))
  expect_equal(t1$ci.high, unname(sapply(t2, "[[", "conf.int")[2,]))
  expect_equal(t1$stderr, unname((sapply(t2, "[[", "estimate")-mus) / sapply(t2, "[[", "statistic")))
  expect_equal(t1$df, unname(sapply(t2, "[[", "parameter")))
  expect_equal(t1$null.mean, unname(sapply(t2, "[[", "null.value")))
  expect_equal(t1$conf.level, unname(sapply(t2, function(x) attributes(x$conf.int)[[1]])))
  expect_equal(t1$alternative, unname(sapply(t2, "[[", "alternative")))
})

test_that("unsolvable situations produce warnings", {
  expect_warning(mat_ttest_single(NA_integer_), "1 of the rows had less than 2 'x' observations")
  expect_warning(mat_ttest_single(matrix(c(1,2), ncol=1)), "2 of the rows had less than 2 'x' observations")
  expect_warning(mat_ttest_single(c(1,1,1,1,1,1)), "1 of the rows were essentially constant")
})

test_that("invalid x values produce errors", {
  expect_error(mat_ttest_single(matrix(LETTERS[1:10], ncol=2)), "'x' must be a numeric matrix or a vector")
  expect_error(mat_ttest_single(list(1:10)), "'x' must be a numeric matrix or a vector")
  expect_error(mat_ttest_single(NA), "'x' must be a numeric matrix or a vector")
  expect_error(mat_ttest_single(NULL), "'x' must be a numeric matrix or a vector")
  expect_error(mat_ttest_single(), 'argument "x" is missing, with no default')
})

test_that("invalid 'alternative' values produce errors", {
  expect_error(mat_ttest_single(c(1,2), alternative=NULL), "'alternative' must be a character vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_single(c(1,2), alternative=NA), "all 'alternative' values must be in: two.sided, less, greater")
  expect_error(mat_ttest_single(c(1,2), alternative="w"), "all 'alternative' values must be in: two.sided, less, greater")
  expect_error(mat_ttest_single(c(1,2), alternative=c("t", "g")), "'alternative' must be a character vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_single(c(1,2), alternative=0), "all 'alternative' values must be in: two.sided, less, greater")
})

test_that("invalid 'mu' values produce errors", {
  expect_error(mat_ttest_single(c(1,2), mu=NULL), "'mu' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_single(c(1,2), mu=NA), "'mu' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_single(c(1,2), mu="a"), "'mu' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_single(c(1,2), mu=c(NA, 1)), "'mu' must be a numeric vector with length 1 or nrow\\(x\\)")
})

test_that("invalid 'conf.level' values produce errors", {
  expect_error(mat_ttest_single(c(1,2), conf.level=NULL), "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_single(c(1,2), conf.level=NA), "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_single(c(1,2), conf.level=NA_integer_), "all 'conf.level' values must be between: 0 and 1")
  expect_error(mat_ttest_single(c(1,2), conf.level="a"), "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_single(c(1,2), conf.level=c(NA, 1)), "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_single(c(1,2), conf.level=1.1), "all 'conf.level' values must be between: 0 and 1")
})



