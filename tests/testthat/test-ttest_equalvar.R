context("ttest_equalvar")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_ttest_equalvar <- function(mat1, mat2, alt="two.sided", mu=0, conf=0.95) {
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  res <- Map(t.test, split(mat1, row(mat1)), split(mat2, row(mat2)),
             alternative=alt, mu=mu, conf=conf, var.equal=TRUE
             )

  # for pooled var
  xvar <- apply(mat1, 1, var, na.rm=TRUE)
  yvar <- apply(mat2, 1, var, na.rm=TRUE)
  xobs <- rowSums(!is.na(mat1))
  yobs <- rowSums(!is.na(mat2))
  tobs <- xobs+yobs
  pvar <- rep(0, nrow(mat1))
  pvar <- ifelse(xobs > 1, pvar + (xobs-1) * xvar, pvar)
  pvar <- ifelse(yobs > 1, pvar + (yobs-1) * yvar, pvar)
  pvar <- pvar/(xobs+yobs-2)

  data.frame(x.mean=unname(sapply(res, "[[", "estimate")[1,]),
             y.mean=unname(sapply(res, "[[", "estimate")[2,]),
             diff.mean=-unname(sapply(lapply(res, getElement, "estimate"), diff)),
             x.var=xvar,
             y.var=yvar,
             pool.var=pvar,
             x.obs=rowSums(!is.na(mat1)),
             y.obs=rowSums(!is.na(mat2)),
             tot.obs=rowSums(!is.na(mat1))+rowSums(!is.na(mat2)),
             t.statistic=unname(sapply(res, "[[", "statistic")),
             p.value=unname(sapply(res, "[[", "p.value")),
             ci.low=unname(sapply(res, "[[", "conf.int")[1,]),
             ci.high=unname(sapply(res, "[[", "conf.int")[2,]),
             stderr=unname((-sapply(lapply(res, getElement, "estimate"), diff) - mu) / sapply(res, "[[", "statistic")),
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

  t1 <- base_ttest_equalvar(X, Y, alts, mus, cfs)
  t2 <- ttest_equalvar(X, Y, alts, mus, cfs)

  expect_equal(t1, t2)
})

test_that("edge cases give equal results", {
  expect_equal(base_ttest_equalvar(1:2, 1), ttest_equalvar(1:2, 1))
  expect_equal(base_ttest_equalvar(1, 1:2), ttest_equalvar(1, 1:2))
  vals1 <- c(1.00000000000001, 1.00000000000002)
  vals2 <- c(1.00000000000001)
  expect_equal(base_ttest_equalvar(vals1, vals2), ttest_equalvar(vals1, vals2))
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("unsolvable situations produce warnings", {
  er <- "1 of the rows were essentially constant"
  expect_warning(ttest_equalvar(c(1,1,1), c(2,2,2)), er)
  er <- "2 of the rows had less than 3 total observations"
  expect_warning(ttest_equalvar(matrix(1:2, ncol=1), matrix(1:2, ncol=1)), er)
  er <-  "1 of the rows had zero 'x' observations"
  expect_warning(ttest_equalvar(NA_integer_, 1:3), er)
  er <-  "1 of the rows had zero 'y' observations"
  expect_warning(ttest_equalvar(1:3, NA_integer_), er)
})

################################################################################
############################### TEST PARAMETERS ################################
################################################################################

################################### X and Y ####################################

test_that("x cannot be missing", {
  er <- 'argument "x" is missing, with no default'
  expect_error(ttest_equalvar(), er)
  expect_error(ttest_equalvar(y=1:2), er)
})

test_that("y cannot be missing", {
  er <- 'argument "y" is missing, with no default'
  expect_error(ttest_equalvar(1:2), er)
})


test_that("invalid x values produce errors", {
  mat <- matrix(1:10, ncol=2)
  er  <- "'x' must be a numeric matrix or a vector"
  expect_error(ttest_equalvar(iris, mat), er)
  expect_error(ttest_equalvar(matrix(LETTERS[1:10], ncol=2), mat), er)
  expect_error(ttest_equalvar(list(1:10), mat), er)
  expect_error(ttest_equalvar(NA, mat), er)
  expect_error(ttest_equalvar(NULL, mat), er)
  expect_error(ttest_equalvar(matrix(nrow=0, ncol=5), mat), er)
  expect_error(ttest_equalvar(matrix(nrow=5, ncol=0), mat), er)
})

test_that("invalid y values produce errors", {
  mat <- matrix(1:10, ncol=2)
  er  <- "'y' must be a numeric matrix or a vector"
  expect_error(ttest_equalvar(mat, iris), er)
  expect_error(ttest_equalvar(mat, matrix(LETTERS[1:10], ncol=2)), er)
  expect_error(ttest_equalvar(mat, list(1:10)), er)
  expect_error(ttest_equalvar(mat, NA), er)
  expect_error(ttest_equalvar(mat, NULL), er)
  expect_error(ttest_equalvar(mat, matrix(nrow=0, ncol=5)), er)
  expect_error(ttest_equalvar(mat, matrix(nrow=5, ncol=0)), er)
})

test_that("x and y have the same number of rows", {
  mat <- matrix(1:10, nrow=5)
  er <- "'x' and 'y' must have the same number of rows"
  expect_error(ttest_equalvar(mat, mat[1:2,]), er)
})

test_that("data.frame values pass when they are all numeric", {
  expect_error(ttest_equalvar(iris[,-5], iris[,-5]), NA)
})

################################# alternative ##################################

test_that("invalid 'alternative' length", {
  er <- "'alternative' must be a character vector with length 1 or nrow\\(x\\)"
  expect_error(ttest_equalvar(1:2, 1:2, alternative=NULL), er)
  expect_error(ttest_equalvar(1:2, 1:2, alternative=c("t", "g")), er)
})

test_that("invalid 'alternative' values", {
  er <- "all 'alternative' values must be in: two.sided, less, greater"
  expect_error(ttest_equalvar(1:2, 1:2, alternative=NA), er)
  expect_error(ttest_equalvar(1:2, 1:2, alternative="w"), er)
  expect_error(ttest_equalvar(1:2, 1:2, alternative=0), er)
})

###################################### mu ######################################

test_that("invalid 'mu' type or length", {
  er <- "'mu' must be a numeric vector with length 1 or nrow\\(x\\)"
  expect_error(ttest_equalvar(1:2, 1:2, mu=NULL), er)
  expect_error(ttest_equalvar(1:2, 1:2, mu=NA), er)
  expect_error(ttest_equalvar(1:2, 1:2, mu="a"), er)
  expect_error(ttest_equalvar(1:2, 1:2, mu=c(NA, 1)), er)
})

test_that("edge cases for mu pass without errors", {
  expect_error(ttest_equalvar(1:2, 1:2, mu=Inf), NA)
  expect_error(ttest_equalvar(1:2, 1:2, mu=-Inf), NA)
  expect_error(ttest_equalvar(1:2, 1:2, mu=NA_integer_), NA)
})

################################## conf.level ##################################

test_that("invalid 'conf.level' type or length", {
  er <- "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)"
  expect_error(ttest_equalvar(1:2, 1:2, conf.level=NULL), er)
  expect_error(ttest_equalvar(1:2, 1:2, conf.level=NA), er)
  expect_error(ttest_equalvar(1:2, 1:2, conf.level="a"), er)
  expect_error(ttest_equalvar(1:2, 1:2, conf.level=c(NA, 1)), er)
})

test_that("invalid 'conf.level' values", {
  er <- "all 'conf.level' values must be between: 0 and 1"
  expect_error(ttest_equalvar(1:2, 1:2, conf.level=NA_integer_), er)
  expect_error(ttest_equalvar(1:2, 1:2, conf.level=1.1), er)
})

