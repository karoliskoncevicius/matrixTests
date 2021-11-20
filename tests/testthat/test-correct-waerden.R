context("correctness of waerden")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

pmcmr_waerden <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad    <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ot <- og <- df <- stat <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- mat[i,!bad]
    grp <- factor(groups[!bad])
    res <- PMCMRplus::vanWaerdenTest(vec, grp)

    df[i]   <- res$parameter
    stat[i] <- res$statistic
    p[i]    <- res$p.value
    ot[i]   <- length(vec)
    og[i]   <- length(unique(grp))
  }

  data.frame(obs.tot=ot, obs.groups=og, df=df, statistic=stat, pvalue=p)
}

################################################################################
########################### TEST ON A RANDOM SAMPLE ############################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  # 25 observations per group
  set.seed(14)
  X <- matrix(rnorm(100000), ncol=100)
  X[sample(length(X), 100)] <- NA
  groups <- sample(c("a","b","c","d"), 100, replace=TRUE)
  groups[sample(length(groups), 10)] <- NA

  t1 <- pmcmr_waerden(X, groups)
  t2 <- suppressWarnings(row_waerden(X, groups))

  expect_equal(t1, t2)

  # 10 observations per group and more NAs
  set.seed(14)
  X <- matrix(rnorm(100000), ncol=40)
  X[sample(length(X), 1000)] <- NA
  groups <- sample(c("a","b","c","d"), 40, replace=TRUE)

  t1 <- pmcmr_waerden(X, groups)
  t2 <- row_waerden(X, groups)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
         100000000000003, 100000000000002, 100000000000003, 100000000000000
         )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)

  # small numbers
  x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
         1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
         )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)
})


test_that("equal variance in all groups give equal results", {
  # two groups with equal variances
  x <- c(1,2,3,2,3,4); g <- c("a","a","a","b","b","b")
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)

  # five groups with equal variances
  x <- c(1:3,2:4,3:5,4:6,5:7); g <- rep(1:5, each=3)
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)
})


test_that("minimal allowed sample size gives equal results", {
  # two groups one with two values and another with one value
  x <- c(1,2,3)
  g <- c("a","a","b")
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)

  # two groups both with one value
  x <- c(1,2)
  g <- c("a","b")
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)

  # five groups one with two values and others with one value
  x <- rnorm(6)
  g <- c("a","a","b","c","d","e")
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)

  # five groups all with one value
  x <- rnorm(5)
  g <- c("a","b","c","d","e")
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)

  # with NAs
  x <- c(1,2,NA,4,5)
  g <- c("a","a","b","b",NA)
  t1 <- pmcmr_waerden(x, g)
  t2 <- suppressWarnings(row_waerden(x, g))
  expect_equal(t1, t2)
})


test_that("groups with one remaining member give equal results", {
  # many groups - one has one remaining element
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[1:2] <- NA
  t1 <- pmcmr_waerden(x, g)
  t2 <- suppressWarnings(row_waerden(x, g))
  expect_equal(t1, t2)

  # many groups - all groups have one remaining element
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[c(1:2,4:5,7:8,10:11)] <- NA
  t1 <- pmcmr_waerden(x, g)
  t2 <- suppressWarnings(row_waerden(x, g))
  expect_equal(t1, t2)
})


test_that("constant values give equal results", {
  # all values are constant
  x <- c(1,1,1,1,1,1); g <- c("a","a","a","b","b","b")
  t1 <- pmcmr_waerden(x, g)
  t2 <- suppressWarnings(row_waerden(x, g))
  expect_equal(t1, t2)

  # within group values are constant
  x <- c(1,1,1,2,2,2); g <- c("a","a","a","b","b","b")
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)

  # one group's values are constant
  x <- c(1,1,1,2,3,4); g <- c("a","a","a","b","b","b")
  t1 <- pmcmr_waerden(x, g)
  t2 <- row_waerden(x, g)
  expect_equal(t1, t2)
})


test_that("various cases of ties give equal results", {
  x <- matrix(rnorm(100*50), nrow=100)
  x[11:100,] <- round(runif(90*50, 0, 10))
  x[1,] <- rep(x[1,1], 50)
  x[2,1:49] <- rep(x[2,1], 49)
  x[3,] <- rep(x[3,1:5])
  x[4,] <- rep(x[3,1:5], each=10)
  x[5,] <- sample(x[5,1:2], ncol(x), replace=TRUE)
  x[6,] <- sample(x[6,1:3], ncol(x), replace=TRUE)
  x[7,] <- sample(x[7,1:4], ncol(x), replace=TRUE)
  x[8,] <- sample(x[8,1:5], ncol(x), replace=TRUE)
  x[9,1:30] <- rep(x[9,1:3], each=10)
  x[10,1:40] <- rep(x[10,1:4], each=10)
  g <- rep(letters[1:5], each=10)
  t1 <- pmcmr_waerden(x, g)
  t2 <- suppressWarnings(row_waerden(x, g))
  expect_equal(t1, t2)
})


test_that("infinite values give equal results", {
  # Inf in each group.
  x <- c(1,1,Inf,2,3,-Inf)
  g <- c("a","a","a","b","b","b")
  t1 <- pmcmr_waerden(x, g)
  t2 <- suppressWarnings(row_waerden(x, g))
  expect_equal(t1, t2)

  # Only infinities
  x <- rep(Inf, 6)
  g <- c("a","a","a","b","b","b")
  t1 <- pmcmr_waerden(x, g)
  t2 <- suppressWarnings(row_waerden(x, g))
  expect_equal(t1, t2)

  # Opposite infinities
  x <- rep(c(-Inf, Inf), each=3)
  g <- c("a","a","a","b","b","b")
  t1 <- pmcmr_waerden(x, g)
  t2 <- suppressWarnings(row_waerden(x, g))
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '3 columns dropped due to missing group information'

  # 2 NAs
  expect_warning(res <- row_waerden(1:11, c(1,1,1,1,NA,NA,NA,2,2,2,2)), wrn, all=TRUE)
  expect_equal(res$obs.tot, 8)
  expect_equal(res$obs.groups, 2)

  # 4 groups with one group dropped because of missing x values
  x <- rnorm(15); x[c(1,2,3)] <- NA
  g <- c(1,1,1,2,2,2,NA,NA,NA,3,3,3,4,4,4)
  expect_warning(res <- row_waerden(x, g), wrn)
  expect_equal(res$obs.tot, 9)
  expect_equal(res$obs.groups, 3)
})


test_that("warning when a row has less than 2 total observations", {
  wrn <- 'row_waerden: 1 of the rows had less than 2 total observations\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # 1 observations in total
  x <- 1; g <- "a"
  expect_warning(res <- row_waerden(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 1)
  expect_equal(res$obs.groups, 1)

  # 1 observation after removing NAs
  x <- 1:10; g <- rep(1:5, each=2)
  x[1:5] <- NA; g[7:10] <- NA
  expect_warning(res <- row_waerden(x, g), wrn, all=FALSE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 1)
  expect_equal(res$obs.groups, 1)
})


test_that("warning when a rows has less than 2 groups", {
  wrn <- 'row_waerden: 1 of the rows had less than 2 groups with enough observations\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # one group
  x <- 1:10; g <- rep(1, 10)
  expect_warning(res <- row_waerden(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 1)

  # two groups but one only has NAs
  x <- 1:10; x[6:10] <- NA
  g <- c(rep(1,5), rep(2,5))
  expect_warning(res <- row_waerden(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 5)
  expect_equal(res$obs.groups, 1)
})


test_that("warning when all values are constant", {
  wrn <- 'row_waerden: 1 of the rows were essentially constant\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # all values are the same
  x <- rep(1,10); g <- rep(1:2, each=5)
  expect_warning(res <- row_waerden(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 2)

  # all values are the same Inf
  x <- rep(Inf,10); g <- rep(1:2, each=5)
  expect_warning(res <- row_waerden(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 2)

  # all values are the same after removing NA groups
  x <- c(rep(1,10),2,3); g <- c(rep(1:5, each=2),NA,NA)
  expect_warning(res <- row_waerden(x, g), wrn, all=FALSE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 5)
})

