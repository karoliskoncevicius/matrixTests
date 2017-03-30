context("ttest_onegroup")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_ttest_onegroup <- function(mat, alt="two.sided", mu=0, conf=0.95) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)
  res <- Map(t.test, split(mat, row(mat)), alternative=alt, mu=mu, conf=conf)

  data.frame(mean.x=unname(sapply(res, "[[", "estimate")),
             var.x=apply(mat, 1, var, na.rm=TRUE),
             obs.x=rowSums(!is.na(mat)),
             t.statistic=unname(sapply(res, "[[", "statistic")),
             p.value=unname(sapply(res, "[[", "p.value")),
             ci.low=unname(sapply(res, "[[", "conf.int")[1,]),
             ci.high=unname(sapply(res, "[[", "conf.int")[2,]),
             stderr=unname((sapply(res, "[[", "estimate")-mu) / sapply(res, "[[", "statistic")),
             df=unname(sapply(res, "[[", "parameter")),
             mean.null=unname(sapply(res, "[[", "null.value")),
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
  X[sample(length(X), 100)] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  mus  <- seq(-1, 1, length.out=nrow(X))
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- base_ttest_onegroup(X, alts, mus, cfs)
  t2 <- ttest_onegroup(X, alts, mus, cfs)

  expect_equal(t1, t2)
})

test_that("edge cases give equal results", {
  expect_equal(base_ttest_onegroup(1:2), ttest_onegroup(1:2))
  vals <- c(1.00000000000001, 1.00000000000002)
  expect_equal(base_ttest_onegroup(vals), ttest_onegroup(vals))
  expect_equal(base_ttest_onegroup(1:10^5), ttest_onegroup(1:10^5))
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("unsolvable situations produce warnings", {
  er <- "1 of the rows had less than 2 'x' observations"
  expect_warning(ttest_onegroup(NA_integer_), er)
  er <- "2 of the rows had less than 2 'x' observations"
  expect_warning(ttest_onegroup(matrix(1:2, ncol=1)), er)
  er <- "1 of the rows were essentially constant"
  expect_warning(ttest_onegroup(c(1,1,1,1,1,1)), er)
})

################################################################################
############################### TEST PARAMETERS ################################
################################################################################

###################################### X #######################################

test_that("x cannot be missing", {
  er <- 'argument "x" is missing, with no default'
  expect_error(ttest_onegroup(), er)
})

test_that("invalid x values produce errors", {
  er <- "'x' must be a numeric matrix or a vector"
  expect_error(ttest_onegroup(iris), er)
  expect_error(ttest_onegroup(matrix(LETTERS[1:10], ncol=2)), er)
  expect_error(ttest_onegroup(list(1:10)), er)
  expect_error(ttest_onegroup(NA), er)
  expect_error(ttest_onegroup(NULL), er)
  expect_error(ttest_onegroup(matrix(nrow=0, ncol=5)), er)
  expect_error(ttest_onegroup(matrix(nrow=5, ncol=0)), er)
})

test_that("data.frame values pass when they are all numeric", {
  expect_error(ttest_onegroup(iris[,-5]), NA)
})

################################# alternative ##################################

test_that("invalid 'alternative' length", {
  er <- "'alternative' must be a character vector with length 1 or nrow\\(x\\)"
  expect_error(ttest_onegroup(1:2, alternative=NULL), er)
  expect_error(ttest_onegroup(1:2, alternative=c("t", "g")), er)
})

test_that("invalid 'alternative' values", {
  er <- "all 'alternative' values must be in: two.sided, less, greater"
  expect_error(ttest_onegroup(1:2, alternative=NA), er)
  expect_error(ttest_onegroup(1:2, alternative="w"), er)
  expect_error(ttest_onegroup(1:2, alternative=0), er)
})

###################################### mu ######################################

test_that("invalid 'mu' type or length", {
  er <- "'mu' must be a numeric vector with length 1 or nrow\\(x\\)"
  expect_error(ttest_onegroup(1:2, mu=NULL), er)
  expect_error(ttest_onegroup(1:2, mu=NA), er)
  expect_error(ttest_onegroup(1:2, mu="a"), er)
  expect_error(ttest_onegroup(1:2, mu=c(NA, 1)), er)
})

test_that("edge cases for mu pass without errors", {
  expect_error(ttest_onegroup(1:2, mu=Inf), NA)
  expect_error(ttest_onegroup(1:2, mu=-Inf), NA)
  expect_error(ttest_onegroup(1:2, mu=NA_integer_), NA)
})

################################## conf.level ##################################

test_that("invalid 'conf.level' type or length", {
  er <- "'conf.level' must be a numeric vector with length 1 or nrow\\(x\\)"
  expect_error(ttest_onegroup(1:2, conf.level=NULL), er)
  expect_error(ttest_onegroup(1:2, conf.level=NA), er)
  expect_error(ttest_onegroup(1:2, conf.level="a"), er)
  expect_error(ttest_onegroup(1:2, conf.level=c(NA, 1)), er)
})

test_that("invalid 'conf.level' values", {
  er <- "all 'conf.level' values must be between: 0 and 1"
  expect_error(ttest_onegroup(1:2, conf.level=NA_integer_), er)
  expect_error(ttest_onegroup(1:2, conf.level=1.1), er)
})


