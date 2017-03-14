context("mat_ttest_equalvar")


test_that("t.test(x) gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(1000), ncol=10)
  Y <- matrix(rnorm(900), ncol=9)
  X[sample(length(X), 100)] <- NA
  Y[sample(length(Y), 100)] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  mus  <- seq(-1, 1, length.out=nrow(X))
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- mat_ttest_equalvar(X, Y, alts, mus, cfs)
  t2 <- Map(t.test, split(X, row(X)), split(Y, row(Y)), var.equal=TRUE,
            alternative=alts, mu=mus, conf.level=cfs
            )

  expect_equal(t1$x.mean, unname(sapply(t2, "[[", "estimate")[1,]))
  expect_equal(t1$y.mean, unname(sapply(t2, "[[", "estimate")[2,]))
  expect_equal(t1$diff.mean, -unname(sapply(lapply(t2, getElement, "estimate"), diff)))
  expect_equal(t1$x.var, apply(X, 1, var, na.rm=TRUE))
  expect_equal(t1$y.var, apply(Y, 1, var, na.rm=TRUE))
  expect_equal(t1$pool.var, (t1$x.var*(t1$x.obs-1) + t1$y.var*(t1$y.obs-1)) / (t1$tot.obs-2))
  expect_equal(t1$x.obs, rowSums(!is.na(X)))
  expect_equal(t1$y.obs, rowSums(!is.na(Y)))
  expect_equal(t1$tot.obs, rowSums(!is.na(X))+rowSums(!is.na(Y)))
  expect_equal(t1$t.statistic, unname(sapply(t2, "[[", "statistic")))
  expect_equal(t1$p.value, unname(sapply(t2, "[[", "p.value")))
  expect_equal(t1$ci.low, unname(sapply(t2, "[[", "conf.int")[1,]))
  expect_equal(t1$ci.high, unname(sapply(t2, "[[", "conf.int")[2,]))
  expect_equal(t1$stderr, unname((-sapply(lapply(t2, getElement, "estimate"), diff) - mus) / sapply(t2, "[[", "statistic")))
  expect_equal(t1$df, unname(sapply(t2, "[[", "parameter")))
  expect_equal(t1$null.mean, unname(sapply(t2, "[[", "null.value")))
  expect_equal(t1$conf.level, unname(sapply(t2, function(x) attributes(x$conf.int)[[1]])))
  expect_equal(t1$alternative, unname(sapply(t2, "[[", "alternative")))
})

test_that("unsolvable situations produce warnings", {
  expect_warning(mat_ttest_equalvar(c(1,1,1), c(2,2,2)), "1 of the rows were essentially constant")
  expect_warning(mat_ttest_equalvar(matrix(c(1,2), ncol=1), matrix(c(1,2), ncol=1)), "2 of the rows had less than 3 total observations")
  expect_warning(mat_ttest_equalvar(NA_integer_, c(1,2,3)), "1 of the rows had zero 'x' observations")
  expect_warning(mat_ttest_equalvar(c(1,2,3), NA_integer_), "1 of the rows had zero 'y' observations")
})

test_that("invalid x and y values produce errors", {
  mat <- matrix(1:10, ncol=2)
  expect_error(mat_ttest_equalvar(matrix(LETTERS[1:10], ncol=2), mat), "'x' must be a numeric matrix or a vector")
  expect_error(mat_ttest_equalvar(mat, matrix(LETTERS[1:10], ncol=2)), "'y' must be a numeric matrix or a vector")
  expect_error(mat_ttest_equalvar(list(1:10), mat), "'x' must be a numeric matrix or a vector")
  expect_error(mat_ttest_equalvar(mat, list(1:10)), "'y' must be a numeric matrix or a vector")
  expect_error(mat_ttest_equalvar(NA, mat), "'x' must be a numeric matrix or a vector")
  expect_error(mat_ttest_equalvar(mat, NA), "'y' must be a numeric matrix or a vector")
  expect_error(mat_ttest_equalvar(NULL, mat), "'x' must be a numeric matrix or a vector")
  expect_error(mat_ttest_equalvar(mat, NULL), "'y' must be a numeric matrix or a vector")
  expect_error(mat_ttest_equalvar(), 'argument "x" is missing, with no default')
  expect_error(mat_ttest_equalvar(mat), 'argument "y" is missing, with no default')
  expect_error(mat_ttest_equalvar(mat, mat[1,,drop=FALSE]), "'x' and 'y' must have the same number of rows")
})

test_that("invalid 'alternative' values produce errors", {
  expect_error(mat_ttest_equalvar(1, 1, alternative=NULL), "'alternative' must be a character vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_equalvar(1, 1, alternative=NA), "all 'alternative' values must be in: two.sided, less, greater")
  expect_error(mat_ttest_equalvar(1, 1, alternative="w"), "all 'alternative' values must be in: two.sided, less, greater")
  expect_error(mat_ttest_equalvar(1, 1, alternative=c("t", "g")), "'alternative' must be a character vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_equalvar(1, 1, alternative=0), "all 'alternative' values must be in: two.sided, less, greater")
})

test_that("invalid 'mu' values produce errors", {
  expect_error(mat_ttest_equalvar(1, 1, mu=NULL), "'mu' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_equalvar(1, 1, mu=NA), "'mu' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_equalvar(1, 1, mu="a"), "'mu' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_equalvar(1, 1, mu=c(NA, 1)), "'mu' must be a numeric vector with length 1 or nrow\\(x\\)")
})

test_that("invalid 'conf.level' values produce errors", {
  expect_error(mat_ttest_equalvar(1, 1, conf.level=NULL), "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_equalvar(1, 1, conf.level=NA), "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_equalvar(1, 1, conf.level=NA_integer_), "all 'conf.level' values must be between: 0 and 1")
  expect_error(mat_ttest_equalvar(1, 1, conf.level="a"), "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_equalvar(1, 1, conf.level=c(NA, 1)), "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)")
  expect_error(mat_ttest_equalvar(1, 1, conf.level=1.1), "all 'conf.level' values must be between: 0 and 1")
})



