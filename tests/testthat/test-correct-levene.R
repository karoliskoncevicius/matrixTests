context("correctness of levene")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

car_levene <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ot <- og <- dft <- dfr <- fst <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- mat[i,!bad]
    grp <- factor(groups[!bad])
    res <- car::leveneTest(vec ~ grp, center="mean")

    dft[i] <- res[["Df"]][1]
    dfr[i] <- res[["Df"]][2]
    fst[i] <- res[["F value"]][1]
    p[i]   <- res[["Pr(>F)"]][1]
    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
  }

  data.frame(obs.tot=ot, obs.groups=og, df.between=dft, df.within=dfr,
             statistic=fst, pvalue=p
             )
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

  t1 <- car_levene(X, groups)
  t2 <- suppressWarnings(row_levene(X, groups))

  expect_equal(t1, t2)

  # 10 observations per group and more NAs
  set.seed(14)
  X <- matrix(rnorm(100000), ncol=40)
  X[sample(length(X), 1000)] <- NA
  groups <- sample(c("a","b","c","d"), 40, replace=TRUE)

  t1 <- car_levene(X, groups)
  t2 <- row_levene(X, groups)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("warning when a infinite values are removed", {
  wrn <- 'row_levene: 1 of the rows had infinite observations that were removed\\.\nFirst occurrence at row 1'

  # -Inf and Inf among observations
  expect_warning(res <- row_levene(c(1,1,Inf,2,3,4), c(1,1,1,2,2,2)), wrn, all=TRUE)
  expect_equal(res$obs.tot, 5)
  expect_equal(res$obs.groups, 2)
})

test_that("extreme numbers give equal results", {
  # big numbers
  x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
         100000000000003, 100000000000002, 100000000000003, 100000000000000
         )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- car_levene(x, g)
  t2 <- row_levene(x, g)
  expect_equal(t1, t2)

  # small numbers
  x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
         1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
         )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)
})


test_that("equal variance in all groups give equal results", {
  # two groups with equal variances
  x <- c(1,2,3,2,3,4); g <- c("a","a","a","b","b","b")
  t1 <- car_levene(x, g)
  t2 <- row_levene(x, g)
  expect_equal(t1, t2)

  # five groups with equal variances
  x <- c(1:3,2:4,3:5,4:6,5:7); g <- rep(1:5, each=3)
  t1 <- car_levene(x, g)
  t2 <- row_levene(x, g)
  expect_equal(t1, t2)
})


test_that("minimal allowed sample size gives equal results", {
  # two groups one with three values and another with one value
  x <- c(1,2,3,4)
  g <- c("a","a","a","b")
  t1 <- car_levene(x, g)
  t2 <- row_levene(x, g)
  expect_equal(t1, t2)

  # two groups one with three values and another with two values
  x <- c(1,2,3,4,5)
  g <- c("a","a","a","b","b")
  t1 <- car_levene(x, g)
  t2 <- row_levene(x, g)
  expect_equal(t1, t2)

  # five groups one with three values and others with one
  x <- c(rnorm(3), 1, 2, 3, 4)
  g <- c("a","a","a","b","c","d","e")
  t1 <- car_levene(x, g)
  t2 <- row_levene(x, g)
  expect_equal(t1, t2)

  # five groups one with three values and others with two or one value
  x <- rnorm(9)
  g <- c("a","a","a","b","b","c","c","d","e")
  t1 <- car_levene(x, g)
  t2 <- row_levene(x, g)
  expect_equal(t1, t2)

  # with NAs
  x <- c(1,2,3,NA,4,5)
  g <- c("a","a","a","b","b",NA)
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)
})


test_that("groups with one remaining member give equal results", {
  # many groups - one has one remaining element
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[1:2] <- NA
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)

  # many groups - all groups except one have one remaining element
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[c(4:5,7:8,10:11)] <- NA
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)
})


test_that("constant values give equal results", {
  # all values are constant
  x <- c(1,1,1,1,1,1); g <- c("a","a","a","b","b","b")
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)

  # within group values are constant
  x <- c(1,1,1,2,2,2); g <- c("a","a","a","b","b","b")
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)

  # one group's values are constant
  x <- c(1,1,1,2,3,4); g <- c("a","a","a","b","b","b")
  t1 <- car_levene(x, g)
  t2 <- row_levene(x, g)
  expect_equal(t1, t2)
})


test_that("cases when some groups are constant after residuals give equal results", {
  # two groups - one has all constant values
  x <- c(rnorm(3), 3, 3, 3)
  g <- rep(letters[1:2], each=3)
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)

  # many groups - all but one have constant values
  x <- c(rnorm(3), 2, 3,3, 4,4,4, 5,5,5,5)
  g <- rep(letters[1:5], c(3,1,2,3,4))
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)

  # many groups - all but one have constant values (with NAs)
  x <- c(rnorm(3), 2,NA, 3,3,NA, 4,4,4,NA, 5,5,5,5,NA)
  g <- rep(letters[1:5], c(3,2,3,4,5))
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)

  # two groups - one has all constant values after residuals
  x <- c(rnorm(4), 3, 3, 2, 2)
  g <- rep(letters[1:2], each=4)
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)

  # many groups - all but one have constant values after residuals
  x <- c(rnorm(3), 2, 0,3, 4,4,4, 6,6,5,5)
  g <- rep(letters[1:5], c(3,1,2,3,4))
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)

  # many groups - all but one have constant values after residuals (with NAs)
  x <- c(rnorm(3), 2,NA, 0,3,NA, 4,4,4,NA, 6,6,5,5,NA)
  g <- rep(letters[1:5], c(3,2,3,4,5))
  t1 <- car_levene(x, g)
  t2 <- suppressWarnings(row_levene(x, g))
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '3 columns dropped due to missing group information'

  # 2 NAs
  expect_warning(res <- row_levene(1:11, c(1,1,1,1,NA,NA,NA,2,2,2,2)), wrn, all=TRUE)
  expect_equal(res$obs.tot, 8)
  expect_equal(res$obs.groups, 2)

  # 4 groups with one group dropped because of missing x values
  x <- rnorm(15); x[c(1,2,3)] <- NA
  g <- c(1,1,1,2,2,2,NA,NA,NA,3,3,3,4,4,4)
  expect_warning(res <- row_levene(x, g), wrn)
  expect_equal(res$obs.tot, 9)
  expect_equal(res$obs.groups, 3)
})


test_that("warning when a rows has less than 2 groups", {
  wrn <- 'row_levene: 1 of the rows had less than 2 groups with enough observations\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # one group
  x <- 1:10; g <- rep(1, 10)
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 1)

  # two groups but one only has NAs
  x <- 1:10; x[6:10] <- NA
  g <- c(rep(1,5), rep(2,5))
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 5)
  expect_equal(res$obs.groups, 1)
})


test_that("warning when a row has less than 3 observations per group", {
  wrn <- 'row_levene: 1 of the rows had no groups with at least 3 observations\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # 10 groups 10 observations
  x <- 1:10; g <- 1:10
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 10)

  # 5 groups 2 observations each
  x <- 1:10; g <- rep(1:5, each=2)
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 5)

  # two groups 4 observations but one is NA
  x <- rnorm(4); x[4] <- NA
  g <- c("a", "b", "b","b")
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 3)
  expect_equal(res$obs.groups, 2)
})


test_that("warning when all values are constant after taking residuals", {
  wrn <- 'row_levene: 1 of the rows had zero within group variance of absolute residuals from the mean\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # two groups - all constant values
  x <- c(0,0,0,0,0,0); g <- c(1,1,1,2,2,2)
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 6)
  expect_equal(res$obs.groups, 2)

  # three groups - constant values + NAs
  x <- c(1,1,1,NA,2,2,2,NA,3,3,3,NA); g <- rep(1:3, each=4)
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 9)
  expect_equal(res$obs.groups, 3)

  # values not constant, but become constant after residuals
  x <- c(1,1,2,2,3,3,4,4,5,5,6,6); g <- rep(1:3, each=4)
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 12)
  expect_equal(res$obs.groups, 3)

  # same as above but different sizes of groups
  x <- c(1,1,1,2,2,2,3,3,4,4,5,6,7); g <- rep(1:4, c(6,4,2,1))
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 13)
  expect_equal(res$obs.groups, 4)

  # same but with each group having additional NA value
  x <- c(1,1,1,2,2,2,NA,3,3,4,4,NA,5,6,NA,7,NA); g <- rep(1:4, c(7,5,3,2))
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 13)
  expect_equal(res$obs.groups, 4)

  # same but with each group having non-constant value that is removed after NA
  x <- c(1,1,1,2,2,2,9,3,3,4,4,9,5,6,9,7,9); g <- rep(1:4, c(7,5,3,2)); g[c(7,12,15,17)] <- NA
  expect_warning(res <- row_levene(x, g), wrn)
  expect_equal(res$obs.tot, 13)
  expect_equal(res$obs.groups, 4)
})


test_that("warning when variance is close to machine precision", {
  wrn <- 'row_levene: 1 of the rows had essentially constant absolute residuals from the mean: results might be unreliable\\.\nFirst occurrence at row 1'

  # two groups - constant values within group
  x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
         1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
         )
  g <- c(rep("a", 4), rep("b", 4))
  expect_warning(res <- row_levene(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 8)
  expect_equal(res$obs.groups, 2)
})

