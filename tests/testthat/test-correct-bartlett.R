context("correctness of bartlett")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_bartlett <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ng <- nt <- vt <- ks <- p <- df <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])

    badgr <- names(which(table(grp)==1))
    bad   <- grp %in% badgr
    vec   <- vec[!bad]
    grp   <- grp[!bad]

    res <- bartlett.test(vec, grp)

    ng[i] <- length(unique(grp))
    nt[i] <- length(vec)
    vt[i] <- sum((tapply(vec, grp, length)-1) * tapply(vec, grp, var, na.rm=TRUE)) / (nt[i]-ng[i])
    ks[i] <- res$statistic
    p[i]  <- res$p.value
    df[i]  <- res$parameter
  }

  data.frame(obs.tot=nt, obs.groups=ng, var.pooled=vt, df=df, statistic.chsq=ks,
             p.value=p, stringsAsFactors=FALSE
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

  t1 <- base_bartlett(X, groups)
  t2 <- row.bartlett(X, groups)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  x <- cbind(100000000000004, 100000000000002, 100000000000003, 100000000000000,
             100000000000003, 100000000000002, 100000000000003, 100000000000000
             )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- base_bartlett(x, g)
  t2 <- row.bartlett(x, g)
  expect_equal(t1, t2)

  # small numbers
  x <- cbind(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
             1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
             )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- base_bartlett(x, g)
  t2 <- row.bartlett(x, g)
  expect_equal(t1, t2)
})


test_that("constant values give equal results", {
  # all values are constant
  x <- c(1,1,1,1); g <- c("a","a","b","b")
  t1 <- base_bartlett(x, g)
  t2 <- suppressWarnings(row.bartlett(x, g))
  expect_equal(t1, t2)

  # within group values are constant
  x <- c(1,1,2,2); g <- c("a","a","b","b")
  t1 <- base_bartlett(x, g)
  t2 <- suppressWarnings(row.bartlett(x, g))
  expect_equal(t1, t2)

  # one group's values are constant
  x <- c(1,1,2,3); g <- c("a","a","b","b")
  t1 <- base_bartlett(x, g)
  t2 <- suppressWarnings(row.bartlett(x, g))
  expect_equal(t1, t2)
})


test_that("minimum allowed sample sizes give equal results", {
  # 2 samples 2 groups
  x <- rnorm(4); g <- c("a", "a", "b", "b")
  t1 <- base_bartlett(x, g)
  t2 <- row.bartlett(x, g)
  expect_equal(t1, t2)

  # 2 samples 2 groups with NAs
  x <- rnorm(6); x[c(3,4)] <- NA
  g <- c("a","a","a","b","b","b")
  t1 <- base_bartlett(x, g)
  t2 <- row.bartlett(x, g)
  expect_equal(t1, t2)
})


test_that("groups with one element remaining are dropped correctly", {
  # dropping one group
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[1:2] <- NA
  t1 <- base_bartlett(x, g)
  t2 <- suppressWarnings(row.bartlett(x, g))
  expect_equal(t1, t2)

  # dropping two groups with NAs
  x <- rnorm(12); x[5] <- NA
  g <- rep(letters[1:4], each=3); g[c(1,2,4)] <- NA
  t1 <- base_bartlett(x, g)
  t2 <- suppressWarnings(row.bartlett(x, g))
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning when rows have less than 2 groups", {
  wrn <- '1 of the rows had less than 2 groups with enough observations\\. First occurrence at row 1'
  nacolumns <- c("statistic.chsq", "p.value")

  # one observation per group
  expect_warning(res <- row.bartlett(1:10, 1:10), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 0)
  expect_equal(res$obs.tot, 0)

  # one complete observation per group and other are NAs
  x <- c(1,NA,2,NA,3,NA,4,NA,5,NA)
  g <- rep(1:5, each=2)
  expect_warning(res <- row.bartlett(x, g), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 0)

  # all observations in one group
  expect_warning(res <- row.bartlett(1:10, rep(1,10)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)

  # All groups except one have 1 element and one NA
  x <- c(1, 2, 3, NA, 4, NA, 5, NA)
  g <- rep(c("A","B","C","D"), each=2)
  expect_warning(res <- row.bartlett(c(x), g), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)
  expect_equal(res$obs.tot, 2)

  # some of the values in groups vector are NA
  x <- c(1, 2, 3, 4, 4, NA, 5, NA)
  g <- rep(c("A","B","C","D"), each=2); g[1] <- NA
  expect_warning(res <- row.bartlett(c(x), g), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)
  expect_equal(res$obs.tot, 2)
})


test_that("warning when some groups have less than 2 observations", {
  wrn <- '1 of the rows had groups with less than 2 observations: those groups were removed\\. First occurrence at row 1'

  # one group is dropped because of one remaining value
  expect_warning(res <- row.bartlett(rnorm(5), c(1,1,2,2,3)), wrn, all=TRUE)
  expect_equal(res$obs.groups, 2)
  expect_equal(res$obs.tot, 4)
  expect_true(all(!is.na(res)))

  # one group is dropped because of zero remaining values
  x <- rnorm(6); x[5:6] <- NA
  g <- c("A","A","B","B","C","C")
  expect_warning(res <- row.bartlett(x, g), wrn, all=TRUE)
  expect_equal(res$obs.groups, 2)
  expect_equal(res$obs.tot, 4)
  expect_true(all(!is.na(res)))

  # 2 groups are dropped because of NA values
  x <- rnorm(8); x[1] <- NA
  g <- rep(c("A","B","C","D"), each=2); g[3] <- NA
  expect_warning(res <- row.bartlett(x, g), wrn)
  expect_equal(res$obs.groups, 2)
  expect_equal(res$obs.tot, 4)
  expect_true(all(!is.na(res)))
})

test_that("warning when none of the groups have variance", {
  wrn <- '1 of the rows had zero variance in all of the groups\\. First occurrence at row 1'
  nacolumns <- c("statistic.chsq", "p.value")

  # all values are constant
  expect_warning(res <- row.bartlett(c(1,1,1,1,1,1), c(1,1,2,2,3,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # all values within one group are constant
  expect_warning(res <- row.bartlett(c(1,1,2,2,1.5,1.5), c(1,1,2,2,3,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # all values are constant except a few for which the group is NA
  expect_warning(res <- row.bartlett(c(1,1,3,1,1,1,2,1), c(1,1,NA,1,2,2,NA,2)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # when one group is dropped because of not enough observations
  expect_warning(res <- row.bartlett(c(1,1,1,1,2), c(1,1,2,2,3)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
})


test_that("warning when some of the groups have no variance", {
  wrn <- '1 of the rows had groups with zero variance: result might be unreliable\\. First occurrence at row 1'

  # one group out of 3
  expect_warning(res <- row.bartlett(c(1,1,1,2,3,2), c(1,1,2,2,3,3)), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic.chsq))

  # two groups out of 3
  expect_warning(res <- row.bartlett(c(1,1,2,2,3,4), c(1,1,2,2,3,3)), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic.chsq))

  # all values within one group plus NA
  expect_warning(res <- row.bartlett(c(1,1,NA,2,3,4,0), c(1,1,1,2,2,3,3)), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic.chsq))

  # when one group is dropped because of not enough observations
  expect_warning(res <- row.bartlett(c(1,1,2,3,4), c(1,1,2,2,3)), wrn)
  expect_true(is.infinite(res$statistic.chsq))
})

