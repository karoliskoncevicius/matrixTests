context("correctness of oneway_welch")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_oneway_welch <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ot <- og <- dft <- dfr <- fst <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])

    badgr <- names(which(table(grp)==1))
    bad   <- grp %in% badgr
    vec   <- vec[!bad]
    grp   <- grp[!bad]

    res <- oneway.test(vec ~ grp)

    dft[i] <- res$parameter[1]
    # fix to return infinity when there are groups with 0 variance
    dfr[i] <- ifelse(all(tapply(vec, grp, var)==0), Inf, res$parameter[2])
    fst[i] <- res$statistic
    p[i]   <- res$p.value
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
  set.seed(14)
  X <- matrix(rnorm(100000), ncol=100)
  X[sample(length(X), 100)] <- NA
  groups <- sample(c("a","b","c","d"), 100, replace=TRUE)
  groups[sample(length(groups), 10)] <- NA

  t1 <- base_oneway_welch(X, groups)
  t2 <- suppressWarnings(row.oneway.welch(X, groups))

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
  t1 <- base_oneway_welch(x, g)
  t2 <- row.oneway.welch(x, g)
  expect_equal(t1, t2)

  # small numbers
  x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
         1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
         )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- base_oneway_welch(x, g)
  t2 <- row.oneway.welch(x, g)
  expect_equal(t1, t2)
})


test_that("constant values give equal results", {
  # all values are constant
  x <- c(1,1,1,1); g <- c("a","a","b","b")
  t1 <- base_oneway_welch(x, g)
  t2 <- suppressWarnings(row.oneway.welch(x, g))
  expect_equal(t1, t2)

  # within group values are constant
  x <- c(1,1,2,2); g <- c("a","a","b","b")
  t1 <- base_oneway_welch(x, g)
  t2 <- suppressWarnings(row.oneway.welch(x, g))
  expect_equal(t1, t2)

  # one group's values are constant
  x <- c(1,1,2,3); g <- c("a","a","b","b")
  t1 <- base_oneway_welch(x, g)
  t1$df.within <- 1; t1$statistic <- Inf; t1$pvalue <- 0 # NOTE: fix to match
  t2 <- suppressWarnings(row.oneway.welch(x, g))
  expect_equal(t1, t2)
})


test_that("minimal allowed sample size gives equal results", {
  # two groups with two values each
  x <- rnorm(4); g <- c("a","a","b","b")
  t1 <- base_oneway_welch(x, g)
  t2 <- suppressWarnings(row.oneway.welch(x, g))
  expect_equal(t1, t2)

  # with NAs
  x <- c(1,2,3,4,NA,5)
  g <- c("a","a","b","b","b",NA)
  t1 <- base_oneway_welch(x, g)
  t2 <- suppressWarnings(row.oneway.welch(x, g))
  expect_equal(t1, t2)
})


test_that("groups with one element remaining are dropped correctly", {
  # dropping one group
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[1:2] <- NA
  t1 <- base_oneway_welch(x, g)
  t2 <- suppressWarnings(row.oneway.welch(x, g))
  expect_equal(t1, t2)

  # dropping two groups with NAs
  x <- rnorm(12); x[5] <- NA
  g <- rep(letters[1:4], each=3); g[c(1,2,4)] <- NA
  t1 <- base_oneway_welch(x, g)
  t2 <- suppressWarnings(row.oneway.welch(x, g))
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '2 columns dropped due to missing group information'

  # 2 NAs
  expect_warning(res <- row.oneway.welch(1:10, c(1,1,1,1,NA,NA,2,2,2,2)), wrn, all=TRUE)
  expect_equal(res$obs.tot, 8)
  expect_equal(res$obs.groups, 2)

  # 4 groups with NA values
  x <- rnorm(10); x[c(1,2)] <- NA
  g <- c(1,1,1,1,NA,NA,3,3,4,4)
  expect_warning(res <- row.oneway.welch(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 6)
  expect_equal(res$obs.groups, 3)
})


test_that("warning when a rows has less than 2 groups", {
  wrn <- '1 of the rows had less than 2 groups with enough observations\\. First occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # one group
  x <- 1:10; g <- rep(1, 10)
  expect_warning(res <- row.oneway.welch(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 1)

  # two groups but one only has NAs
  x <- 1:10; x[6:10] <- NA
  g <- c(rep(1,5), rep(2,5))
  expect_warning(res <- row.oneway.welch(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 5)
  expect_equal(res$obs.groups, 1)
})


test_that("warning when a row has groups with less than 2 observations", {
  wrn <- '1 of the rows had groups with less than 2 observations: those groups were removed\\. First occurrence at row 1'

  # 3 groups with one having only one observation
  x <- rnorm(5); g <- c(1,1,2,2,3)
  expect_warning(res <- row.oneway.welch(x, g), wrn, all=TRUE)
  expect_equal(suppressWarnings(row.oneway.welch(x, g)), row.oneway.welch(x[-5], g[-5]))
  expect_equal(res$obs.tot, 4)
  expect_equal(res$obs.groups, 2)

  # 4 groups with 1 having no observations
  x <- rnorm(10); x[9:10] <- NA; g <- c(1,1,1,2,2,2,3,3,4,4)
  expect_warning(res <- row.oneway.welch(x, g), wrn, all=TRUE)
  expect_equal(suppressWarnings(row.oneway.welch(x, g)), row.oneway.welch(x[-(9:10)], g[-(9:10)]))
  expect_equal(res$obs.tot, 8)
  expect_equal(res$obs.groups, 3)
})


test_that("warning when all values within each group are constant", {
  wrn <- '1 of the rows had zero variance in all of the groups\\. First occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # two groups - all values are constant
  x <- c(1,1,1,1); g <- c(1,1,2,2)
  expect_warning(res <- row.oneway.welch(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 4)
  expect_equal(res$obs.groups, 2)

  # two groups - constant values within each group
  x <- c(3,3,0,0); g <- c(1,1,2,2)
  expect_warning(res <- row.oneway.welch(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 4)
  expect_equal(res$obs.groups, 2)
})


test_that("warning when one of the groups has constant values", {
  wrn <- '1 of the rows had groups with zero variance: result might be unreliable\\. First occurrence at row 1'

  # two groups - one with constant values
  x <- c(1,2,1,1); g <- c(1,1,2,2)
  expect_warning(res <- row.oneway.welch(x, g), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic))

  # three groups - constant values within each group + NAs
  x <- c(1,2,3,4,5,5,NA); g <- c("a","a","b","b","c","c","c")
  expect_warning(res <- row.oneway.welch(x, g), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic))
})


