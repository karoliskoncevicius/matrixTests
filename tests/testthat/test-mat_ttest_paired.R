context("mat_ttest_paired")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

ttest_paired <- function(mat1, mat2, alt="two.sided", mu=0, conf=0.95) {
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  res <- Map(t.test, split(mat1, row(mat1)), split(mat2, row(mat2)),
             alternative=alt, mu=mu, conf=conf, paired=TRUE
             )

  # for pooled var
  data.frame(x.mean=rowMeans(mat1, na.rm=TRUE),
             y.mean=rowMeans(mat2, na.rm=TRUE),
             diff.mean=unname(sapply(res, "[[", "estimate")),
             x.var=apply(mat1, 1, var, na.rm=TRUE),
             y.var=apply(mat2, 1, var, na.rm=TRUE),
             diff.var=apply(mat1-mat2, 1, var, na.rm=TRUE),
             x.obs=rowSums(!is.na(mat1)),
             y.obs=rowSums(!is.na(mat2)),
             pair.obs=rowSums(!is.na(mat1-mat2)),
             t.statistic=unname(sapply(res, "[[", "statistic")),
             p.value=unname(sapply(res, "[[", "p.value")),
             ci.low=unname(sapply(res, "[[", "conf.int")[1,]),
             ci.high=unname(sapply(res, "[[", "conf.int")[2,]),
             stderr=unname((sapply(res, "[[", "estimate") - mu) / sapply(res, "[[", "statistic")),
             df=unname(sapply(res, "[[", "parameter")),
             null.mean=unname(sapply(res, "[[", "null.value")),
             conf.level=unname(sapply(res, function(x) attributes(x$conf.int)[[1]])),
             alternative=unname(sapply(res, "[[", "alternative")),
             stringsAsFactors=FALSE
             )
}

################################################################################
######################### TEST CONSISTENCY WITH T.TEST #########################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(1000), ncol=10)
  Y <- matrix(rnorm(1000), ncol=10)
  X[sample(length(X), 100)] <- NA
  Y[sample(length(Y), 100)] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  mus  <- seq(-1, 1, length.out=nrow(X))
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- ttest_paired(X, Y, alts, mus, cfs)
  t2 <- mat_ttest_paired(X, Y, alts, mus, cfs)

  expect_equal(t1, t2)
})

test_that("edge cases give equal results", {
  vals1 <- c(1.00000000000001, 1.00000000000002)
  vals2 <- c(1.00000000000001, 1.00000000000003)
  expect_equal(ttest_paired(vals1, vals2), mat_ttest_paired(vals1, vals2))
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("unsolvable situations produce warnings", {
  er <- "1 of the rows were essentially constant"
  expect_warning(mat_ttest_paired(c(1,1,1), c(2,2,2)), er)
  er <- "2 of the rows had less than 2 paired observations"
  expect_warning(mat_ttest_paired(matrix(1:2, ncol=1), matrix(1:2, ncol=1)), er)
  er1 <- "1 of the rows had less than 2 paired observations"
  er2 <- "1 of the rows had none 'y' observations"
  expect_warning(mat_ttest_paired(1:2, c(NA_integer_, NA_integer_)), er1)
  expect_warning(mat_ttest_paired(1:2, c(NA_integer_, NA_integer_)), er2)
})

################################################################################
############################### TEST PARAMETERS ################################
################################################################################

################################### X and Y ####################################

test_that("x cannot be missing", {
  er <- 'argument "x" is missing, with no default'
  expect_error(mat_ttest_paired(), er)
  expect_error(mat_ttest_paired(y=1:2), er)
})

test_that("y cannot be missing", {
  er <- 'argument "y" is missing, with no default'
  expect_error(mat_ttest_paired(1:2), er)
})

test_that("invalid x values produce errors", {
  mat <- matrix(1:10, ncol=2)
  er  <- "'x' must be a numeric matrix or a vector"
  expect_error(mat_ttest_paired(iris, mat), er)
  expect_error(mat_ttest_paired(matrix(LETTERS[1:10], ncol=2), mat), er)
  expect_error(mat_ttest_paired(list(1:10), mat), er)
  expect_error(mat_ttest_paired(NA, mat), er)
  expect_error(mat_ttest_paired(NULL, mat), er)
  expect_error(mat_ttest_paired(matrix(nrow=0, ncol=5), mat), er)
  expect_error(mat_ttest_paired(matrix(nrow=5, ncol=0), mat), er)
})

test_that("invalid y values produce errors", {
  mat <- matrix(1:10, ncol=2)
  er  <- "'y' must be a numeric matrix or a vector"
  expect_error(mat_ttest_paired(mat, iris), er)
  expect_error(mat_ttest_paired(mat, matrix(LETTERS[1:10], ncol=2)), er)
  expect_error(mat_ttest_paired(mat, list(1:10)), er)
  expect_error(mat_ttest_paired(mat, NA), er)
  expect_error(mat_ttest_paired(mat, NULL), er)
  expect_error(mat_ttest_paired(mat, matrix(nrow=0, ncol=5)), er)
  expect_error(mat_ttest_paired(mat, matrix(nrow=5, ncol=0)), er)
})

test_that("x and y have the same number of rows and columns", {
  mat <- matrix(1:20, nrow=5)
  er <- "'x' and 'y' must have the same number of rows"
  expect_error(mat_ttest_paired(mat, mat[1:2,]), er)
  er <- "'x' and 'y' must have the same number of columns"
  expect_error(mat_ttest_paired(mat, mat[,1:2]), er)
})

test_that("data.frame values pass when they are all numeric", {
  expect_error(mat_ttest_paired(iris[,-5], rev(iris[,-5])), NA)
})

################################# alternative ##################################

test_that("invalid 'alternative' length", {
  er <- "'alternative' must be a character vector with length 1 or nrow\\(x\\)"
  expect_error(mat_ttest_paired(1:2, 1:2, alternative=NULL), er)
  expect_error(mat_ttest_paired(1:2, 1:2, alternative=c("t", "g")), er)
})

test_that("invalid 'alternative' values", {
  er <- "all 'alternative' values must be in: two.sided, less, greater"
  expect_error(mat_ttest_paired(1:2, 1:2, alternative=NA), er)
  expect_error(mat_ttest_paired(1:2, 1:2, alternative="w"), er)
  expect_error(mat_ttest_paired(1:2, 1:2, alternative=0), er)
})

###################################### mu ######################################

test_that("invalid 'mu' type or length", {
  er <- "'mu' must be a numeric vector with length 1 or nrow\\(x\\)"
  expect_error(mat_ttest_paired(1:2, 1:2, mu=NULL), er)
  expect_error(mat_ttest_paired(1:2, 1:2, mu=NA), er)
  expect_error(mat_ttest_paired(1:2, 1:2, mu="a"), er)
  expect_error(mat_ttest_paired(1:2, 1:2, mu=c(NA, 1)), er)
})

test_that("edge cases for mu pass without errors", {
  expect_error(mat_ttest_paired(1:2, 2:1, mu=Inf), NA)
  expect_error(mat_ttest_paired(1:2, 2:1, mu=-Inf), NA)
  expect_error(mat_ttest_paired(1:2, 2:1, mu=NA_integer_), NA)
})

################################## conf.level ##################################

test_that("invalid 'conf.level' type or length", {
  er <- "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)"
  expect_error(mat_ttest_paired(1:2, 1:2, conf.level=NULL), er)
  expect_error(mat_ttest_paired(1:2, 1:2, conf.level=NA), er)
  expect_error(mat_ttest_paired(1:2, 1:2, conf.level="a"), er)
  expect_error(mat_ttest_paired(1:2, 1:2, conf.level=c(NA, 1)), er)
})

test_that("invalid 'conf.level' values", {
  er <- "all 'conf.level' values must be between: 0 and 1"
  expect_error(mat_ttest_paired(1:2, 1:2, conf.level=NA_integer_), er)
  expect_error(mat_ttest_paired(1:2, 1:2, conf.level=1.1), er)
})

