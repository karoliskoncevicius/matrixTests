context("Warnings")

################################################################################
################################ TTEST_ONEGROUP ################################
################################################################################

test_that("rows should have more than 1 available observation", {
  wrn <- '1 of the rows had less than 2 "x" observations'
  expect_warning(tst <- ttest_onegroup(1)$t.statistic, wrn)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_onegroup(NA_integer_)$t.statistic, wrn)
  expect_true(is.na(tst))
  wrn <- '10 of the rows had less than 2 "x" observations'
  expect_warning(tst <- ttest_onegroup(matrix(1:10, ncol=1))$t.statistic, wrn)
  expect_true(all(is.na(tst)))
  mat <- matrix(rnorm(100), ncol=10)
  mat[4:6,1:9] <- NA
  wrn <- '3 of the rows had less than 2 "x" observations'
  expect_warning(tst <- ttest_onegroup(mat)$t.statistic, wrn)
  expect_true(all(is.na(tst[4:6])))
})


test_that("rows should not be constant", {
  wrn <- '1 of the rows were essentially constant'
  expect_warning(tst <- ttest_onegroup(c(1,1))$t.statistic, wrn)
  expect_true(is.na(tst))
  wrn <- '10 of the rows were essentially constant'
  expect_warning(tst <- ttest_onegroup(matrix(1:10, ncol=10, nrow=10))$t.statistic, wrn)
  expect_true(all(is.na(tst)))
  mat <- matrix(1:10, ncol=10, nrow=10)
  mat[4:6,1:8] <- NA
  mat[c(1:3,7:10),10] <- 11
  wrn <- '3 of the rows were essentially constant'
  expect_warning(tst <- ttest_onegroup(mat)$t.statistic, wrn)
  expect_true(all(is.na(tst[4:6])))
})

################################################################################
################################ TTEST_EQUALVAR ################################
################################################################################

test_that("rows should have at least 3 observations", {
  wrn <- '1 of the rows had less than 3 total observations'
  expect_warning(tst <- ttest_equalvar(1, 1)$t.statistic, wrn)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_equalvar(NA_integer_, NA_integer_)$t.statistic, wrn)
  expect_true(is.na(tst))
  wrn <- '10 of the rows had less than 3 total observations'
  expect_warning(tst <- ttest_equalvar(matrix(1:10, ncol=1), matrix(1:10, ncol=1))$t.statistic, wrn)
  expect_true(all(is.na(tst)))
  matX <- matrix(rnorm(20), ncol=2)
  matY <- matrix(rnorm(20), ncol=2)
  matX[4:6,2] <- NA
  matY[4:6,1] <- NA
  matY[1:3,] <- NA
  wrn <- '6 of the rows had less than 3 total observations'
  expect_warning(tst <- ttest_equalvar(matX, matY)$t.statistic, wrn)
  expect_true(all(is.na(tst[1:6])))
  expect_warning(tst <- ttest_equalvar(matY, matX)$t.statistic, wrn)
  expect_true(all(is.na(tst[1:6])))
})

test_that("one group should have at least one observation", {
  wrnX <- '1 of the rows had zero "x" observations'
  wrnY <- '1 of the rows had zero "y" observations'
  expect_warning(tst <- ttest_equalvar(NA_integer_, c(1,1,1))$t.statistic, wrnX)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_equalvar(c(1,1,1), NA_integer_)$t.statistic, wrnY)
  expect_true(is.na(tst))
  matX <- matrix(rnorm(30), ncol=3)
  matY <- matrix(rnorm(20), ncol=2)
  matY[4:6,1] <- NA
  matY[1:3,] <- NA
  wrnX <- '3 of the rows had zero "x" observations'
  wrnY <- '3 of the rows had zero "y" observations'
  expect_warning(tst <- ttest_equalvar(matY, matX)$t.statistic, wrnX)
  expect_true(all(is.na(tst[1:3])))
  expect_warning(tst <- ttest_equalvar(matX, matY)$t.statistic, wrnY)
  expect_true(all(is.na(tst[1:3])))
})

test_that("rows should not be constant", {
  wrn <- '1 of the rows were essentially constant'
  expect_warning(tst <- ttest_equalvar(c(1,1), c(2,2))$t.statistic, wrn)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_equalvar(c(1,1,1), c(1,1,1))$t.statistic, wrn)
  expect_true(is.na(tst))
  matX <- matrix(1:10, ncol=10, nrow=10)
  matY <- matrix(10:1, ncol=10, nrow=10)
  matX[4:6,1:8] <- NA
  matX[c(1:3,7:10),10] <- 11
  wrn <- '3 of the rows were essentially constant'
  expect_warning(tst <- ttest_equalvar(matX, matY)$t.statistic, wrn)
  expect_true(all(is.na(tst[4:6])))
  expect_warning(tst <- ttest_equalvar(matY, matX)$t.statistic, wrn)
  expect_true(all(is.na(tst[4:6])))
})

################################################################################
################################# TTEST_WELCH ##################################
################################################################################

test_that("rows should have at least 2 observations in each group", {
  wrnX <- '1 of the rows had less than 2 "x" observations'
  wrnY <- '1 of the rows had less than 2 "y" observations'
  expect_warning(tst <- ttest_welch(1, 1)$t.statistic, wrnX)
  expect_warning(tst <- ttest_welch(1, 1)$t.statistic, wrnY)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_welch(1, 1:2)$t.statistic, wrnX)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_welch(1:2, 1)$t.statistic, wrnY)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_welch(NA_integer_, NA_integer_)$t.statistic, wrnX)
  expect_warning(tst <- ttest_welch(NA_integer_, NA_integer_)$t.statistic, wrnY)
  expect_true(is.na(tst))
  wrnX <- '10 of the rows had less than 2 "x" observations'
  wrnY <- '10 of the rows had less than 2 "y" observations'
  matX <- matrix(1:20, ncol=2)
  matY <- matrix(1:10, ncol=1)
  expect_warning(tst <- ttest_welch(matY, matY)$t.statistic, wrnX)
  expect_warning(tst <- ttest_welch(matY, matY)$t.statistic, wrnY)
  expect_true(all(is.na(tst)))
  expect_warning(tst <- ttest_welch(matY, matX)$t.statistic, wrnX)
  expect_true(all(is.na(tst)))
  expect_warning(tst <- ttest_welch(matX, matY)$t.statistic, wrnY)
  expect_true(all(is.na(tst)))
  wrnX <- '6 of the rows had less than 2 "x" observations'
  wrnY <- '3 of the rows had less than 2 "y" observations'
  matX <- matrix(rnorm(20), ncol=2)
  matY <- matrix(rnorm(20), ncol=2)
  matX[4:6,1] <- NA
  matX[1:3,] <- NA
  matY[4:6,2] <- NA
  expect_warning(tst <- ttest_welch(matX, matY)$t.statistic, wrnX)
  expect_warning(tst <- ttest_welch(matX, matY)$t.statistic, wrnY)
  expect_true(all(is.na(tst[1:6])))
})

test_that("rows should not be constant", {
  wrn <- '1 of the rows were essentially constant'
  expect_warning(tst <- ttest_welch(c(1,1), c(2,2))$t.statistic, wrn)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_welch(c(1,1,1), c(1,1,1))$t.statistic, wrn)
  expect_true(is.na(tst))
  matX <- matrix(1:10, ncol=10, nrow=10)
  matY <- matrix(10:1, ncol=10, nrow=10)
  matX[4:6,1:8] <- NA
  matX[c(1:3,7:10),10] <- 11
  wrn <- '3 of the rows were essentially constant'
  expect_warning(tst <- ttest_welch(matX, matY)$t.statistic, wrn)
  expect_true(all(is.na(tst[4:6])))
  expect_warning(tst <- ttest_welch(matY, matX)$t.statistic, wrn)
  expect_true(all(is.na(tst[4:6])))
})

################################################################################
################################# TTEST_PAIRED #################################
################################################################################

test_that("rows should have at least 2 paired observations", {
  wrn <- '1 of the rows had less than 2 paired observations'
  expect_warning(tst <- ttest_paired(1, 1)$t.statistic, wrn)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_paired(NA_integer_, NA_integer_)$t.statistic, wrn)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_paired(c(1,NA,2), c(3,4,NA))$t.statistic, wrn)
  expect_true(is.na(tst))
  wrn <- '10 of the rows had less than 2 paired observations'
  expect_warning(tst <- ttest_paired(matrix(1:10, ncol=1), matrix(1:10, ncol=1))$t.statistic, wrn)
  expect_true(all(is.na(tst)))
  matX <- matrix(rnorm(20), ncol=2)
  matY <- matrix(rnorm(20), ncol=2)
  matX[4:6,2] <- NA
  matY[4:6,1] <- NA
  matY[1:3,] <- NA
  wrn <- '6 of the rows had less than 2 paired observations'
  expect_warning(tst <- ttest_paired(matX, matY)$t.statistic, wrn)
  expect_true(all(is.na(tst[1:6])))
  expect_warning(tst <- ttest_paired(matY, matX)$t.statistic, wrn)
  expect_true(all(is.na(tst[1:6])))
})

test_that("rows should not be constant", {
  wrn <- '1 of the rows were essentially constant'
  expect_warning(tst <- ttest_paired(c(1,1), c(2,2))$t.statistic, wrn)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_paired(c(1,1,1), c(1,1,1))$t.statistic, wrn)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_paired(c(1:5), c(1:5))$t.statistic, wrn)
  expect_true(is.na(tst))
  expect_warning(tst <- ttest_paired(c(1:5), c(1,2,NA,4,5))$t.statistic, wrn)
  expect_true(is.na(tst))
  matX <- matrix(1:10, ncol=10, nrow=10)
  matY <- matrix(10:1, ncol=10, nrow=10)
  matX[4:6,1:8] <- NA
  matX[c(1:3,7:10),10] <- 11
  wrn <- '3 of the rows were essentially constant'
  expect_warning(tst <- ttest_paired(matX, matY)$t.statistic, wrn)
  expect_true(all(is.na(tst[4:6])))
  expect_warning(tst <- ttest_paired(matY, matX)$t.statistic, wrn)
  expect_true(all(is.na(tst[4:6])))
})

################################################################################
############################### ONEWAY_EQUALVAR ################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '2 columns dropped due to missing group information'
  expect_warning(obs <- oneway_equalvar(1:10, c(1,1,1,1,NA,NA,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, 8)
  expect_warning(obs <- oneway_equalvar(1:10, c(1,1,1,1,NaN,NaN,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, 8)
  expect_warning(obs <- oneway_equalvar(matrix(1:20, nrow=2), c(1,1,1,1,NA,NaN,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, c(8,8))
})

test_that("rows should have at least 2 groups", {
  wrn <- '1 of the rows had less than 2 groups with enough observations'
  expect_warning(fst <- oneway_equalvar(1:10, rep(1,10))$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_equalvar(1:10, c(rep(1,9),NA))$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_equalvar(c(1:8,NA,NA), c(rep(1,8),2,2))$F.statistic, wrn)
  expect_true(is.na(fst))
  wrn <- '6 of the rows had less than 2 groups with enough observations'
  mat <- matrix(rnorm(100), ncol=10)
  mat[1:2,c(1:3,7)] <- NA
  mat[3:4,c(1:6,7)] <- NA
  mat[5:6,c(3:6,7)] <- NA
  grp <- c(1,1,1,2,2,2,3,NA,NA,NA)
  expect_warning(fst <- oneway_equalvar(mat, grp)$F.statistic, wrn)
  expect_true(all(is.na(fst[1:6])))
})

test_that("rows should have more than 1 observation per group", {
  wrn <- '1 of the rows had one observation per group'
  expect_warning(fst <- oneway_equalvar(1:2, 1:2)$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_equalvar(c(1,2,3), c(1,2,NA))$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_equalvar(c(1,2,NA), c(1,2,2))$F.statistic, wrn)
  expect_true(is.na(fst))
  wrn <- '3 of the rows had one observation per group'
  mat <- matrix(rnorm(100), ncol=10)
  mat[1,c(1,2,4,6)] <- NA
  mat[2,c(1,2,3,4,6)] <- NA
  mat[3,c(2,3,4,5,7)] <- NA
  grp <- c(1,1,1,2,2,3,3,NA,NA,NA)
  expect_warning(fst <- oneway_equalvar(mat, grp)$F.statistic, wrn)
  expect_true(all(is.na(fst[1:3])))
})

test_that("all values within one group should not be constant", {
  wrn <- '1 of the rows had essentially perfect fit'
  expect_warning(fst <- oneway_equalvar(c(1,1,2,2), c(1,1,2,2))$F.statistic, wrn)
  expect_true(is.infinite(fst))
  expect_warning(fst <- oneway_equalvar(c(1,1,1,1), c(1,1,2,2))$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_equalvar(c(1,1,3,4), c(1,1,NA,2))$F.statistic, wrn)
  expect_true(is.infinite(fst))
})

################################################################################
############################### ONEWAY_EQUALVAR ################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '2 columns dropped due to missing group information'
  expect_warning(obs <- oneway_welch(1:10, c(1,1,1,1,NA,NA,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, 8)
  expect_warning(obs <- oneway_welch(1:10, c(1,1,1,1,NaN,NaN,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, 8)
  expect_warning(obs <- oneway_welch(matrix(1:20, nrow=2), c(1,1,1,1,NA,NaN,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, c(8,8))
})

test_that("rows should have at least 2 groups", {
  wrn <- '1 of the rows had less than 2 groups with enough observations'
  expect_warning(fst <- oneway_welch(1:10, rep(1,10))$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_welch(1:10, c(rep(1,9),NA))$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_welch(c(1:8,NA,NA), c(rep(1,8),2,2))$F.statistic, wrn)
  expect_true(is.na(fst))
  wrn <- '6 of the rows had less than 2 groups with enough observations'
  mat <- matrix(rnorm(90), ncol=9)
  mat[1:2,c(1:3,7)] <- NA
  mat[3:4,c(1:6,7)] <- NA
  mat[5:6,c(3:6,7)] <- NA
  grp <- c(1,1,1,2,2,2,NA,NA,NA)
  expect_warning(fst <- oneway_welch(mat, grp)$F.statistic, wrn)
  expect_true(all(is.na(fst[1:6])))
})

test_that("rows should have more than 1 observation in each group", {
  wrn <- '1 of the rows had less than 2 observations per group'
  expect_warning(fst <- oneway_welch(1:2, 1:2)$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_welch(c(1,2,3), c(1,2,2))$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_welch(c(1,2,3,NA), c(1,1,2,2))$F.statistic, wrn)
  expect_true(is.na(fst))
  wrn <- '4 of the rows had less than 2 observations per group'
  mat <- matrix(rnorm(100), ncol=10)
  mat[1,c(1,2,4,6)] <- NA
  mat[2,c(1,2,3,4,6)] <- NA
  mat[3,c(2,3,4,5,7)] <- NA
  mat[4,6] <- NA
  grp <- c(1,1,1,2,2,3,3,NA,NA,NA)
  expect_warning(fst <- oneway_welch(mat, grp)$F.statistic, wrn)
  expect_true(all(is.na(fst[1:3])))
})

test_that("all values within one group should not be constant", {
  wrn <- '1 of the rows had essentially perfect fit'
  expect_warning(fst <- oneway_welch(c(1,1,2,2), c(1,1,2,2))$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_welch(c(1,1,1,1), c(1,1,2,2))$F.statistic, wrn)
  expect_true(is.na(fst))
  expect_warning(fst <- oneway_welch(c(1,1,3,3,5), c(1,1,2,2,NA))$F.statistic, wrn)
  expect_true(is.na(fst))
})

################################################################################
################################ KRUSKAL_WALLIS ################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '2 columns dropped due to missing group information'
  expect_warning(obs <- kruskalwallis(1:10, c(1,1,1,1,NA,NA,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, 8)
  expect_warning(obs <- kruskalwallis(1:10, c(1,1,1,1,NaN,NaN,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, 8)
  expect_warning(obs <- kruskalwallis(matrix(1:20, nrow=2), c(1,1,1,1,NA,NaN,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, c(8,8))
})

test_that("rows should have at least 2 observations", {
  wrn <- '1 of the rows had less than 2 total observations'
  expect_warning(kst <- kruskalwallis(1, "a")$chsq, wrn)
  expect_true(is.na(kst))
  expect_warning(kst <- kruskalwallis(c(1,NA,NA,NA), c(1,1,2,2))$chsq, wrn)
  expect_true(is.na(kst))
  expect_warning(kst <- kruskalwallis(c(1,2,3,NA), c(NA,NA,1,2))$chsq, wrn)
  expect_true(is.na(kst))
  expect_warning(kst <- kruskalwallis(1, 1)$chsq, wrn)
  expect_true(is.na(kst))
  wrn <- '4 of the rows had less than 2 total observations'
  mat <- matrix(rnorm(100), ncol=10)
  mat[1:2,c(2:4)] <- NA
  mat[3:4,c(1:6)] <- NA
  grp <- c(1,1,1,2,NA,NA,NA,NA,NA,NA)
  expect_warning(kst <- kruskalwallis(mat, grp)$chsq, wrn)
  expect_true(all(is.na(kst[1:4])))
})

test_that("rows should have at least 2 groups", {
  wrn <- '1 of the rows had less than 2 groups with enough observations'
  expect_warning(kst <- kruskalwallis(1:10, rep(1,10))$chsq, wrn)
  expect_true(is.na(kst))
  expect_warning(kst <- kruskalwallis(1:10, c(rep(1,9),NA))$chsq, wrn)
  expect_true(is.na(kst))
  expect_warning(kst <- kruskalwallis(c(1:8,NA,NA), c(rep(1,8),2,2))$chsq, wrn)
  expect_true(is.na(kst))
  wrn <- '6 of the rows had less than 2 groups with enough observations'
  mat <- matrix(rnorm(100), ncol=10)
  mat[1:2,c(4:8)] <- NA
  mat[3:4,c(1:3,7:8)] <- NA
  mat[5:6,c(1:6)] <- NA
  grp <- c(1,1,1,2,2,2,3,3,NA,NA)
  expect_warning(kst <- kruskalwallis(mat, grp)$chsq, wrn)
  expect_true(all(is.na(kst[1:6])))
})

################################################################################
################################### BARTLETT ###################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '2 columns dropped due to missing group information'
  expect_warning(obs <- bartlett(1:10, c(1,1,1,1,NA,NA,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, 8)
  expect_warning(obs <- bartlett(1:10, c(1,1,1,1,NaN,NaN,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, 8)
  expect_warning(obs <- bartlett(matrix(1:20, nrow=2), c(1,1,1,1,NA,NaN,2,2,2,2))$obs.tot, wrn)
  expect_equal(obs, c(8,8))
})

test_that("rows should have at least 2 groups", {
  wrn <- '1 of the rows had less than 2 groups with enough observations'
  expect_warning(kst <- bartlett(1:10, rep(1,10))$chsq.statistic, wrn)
  expect_true(is.na(kst))
  expect_warning(kst <- bartlett(1:10, c(rep(1,9),NA))$chsq.statistic, wrn)
  expect_true(is.na(kst))
  expect_warning(kst <- bartlett(c(1:8,NA,NA), c(rep(1,8),2,2))$chsq.statistic, wrn)
  expect_true(is.na(kst))
  wrn <- '6 of the rows had less than 2 groups with enough observations'
  mat <- matrix(rnorm(100), ncol=10)
  mat[1:2,c(1:3,7)] <- NA
  mat[3:4,c(1:6,7)] <- NA
  mat[5:6,c(3:6,7)] <- NA
  grp <- c(1,1,1,2,2,2,3,3,NA,NA)
  expect_warning(kst <- bartlett(mat, grp)$chsq.statistic, wrn)
  expect_true(all(is.na(kst[1:6])))
})

test_that("all groups should have at least 2 observations", {
  wrn <- '1 of the rows had groups with less than 2 observations'
  expect_warning(ngr <- bartlett(rnorm(5), c(1,1,2,2,3))$obs.groups, wrn)
  expect_equal(ngr, 2)
  expect_warning(ngr <- bartlett(c(rnorm(5),NA), c(1,1,2,2,3,3))$obs.groups, wrn)
  expect_equal(ngr, 2)
  wrn <- '4 of the rows had groups with less than 2 observations'
  mat <- matrix(rnorm(100), ncol=10)
  grp <- c(1,1,1,2,2,2,3,3,4,4)
  mat[1,c(1,2)] <- NA
  mat[2,c(4,5)] <- NA
  mat[3,c(7)] <- NA
  mat[4,c(7,9)] <- NA
  expect_warning(ngr <- bartlett(mat, grp)$obs.groups, wrn)
  expect_equal(ngr, c(3,3,3,2,4,4,4,4,4,4))
})

test_that("all groups should have some variability", {
  wrn <- '1 of the rows had zero variance in all of the groups'
  expect_warning(ksq <- bartlett(c(1,1,2,2), c(1,1,2,2))$chsq.statistic, wrn)
  expect_true(is.na(ksq))
  expect_warning(ksq <- bartlett(c(1,1,2,2,3,3,4), c(1,1,2,2,3,3,NA))$chsq.statistic, wrn)
  expect_true(is.na(ksq))
  wrn <- '1 of the rows had groups with zero variance'
  expect_warning(ksq <- bartlett(c(1,1,2,3), c(1,1,2,2))$chsq.statistic, wrn)
  expect_true(is.infinite(ksq))
  expect_warning(ksq <- bartlett(c(1,1,4,2,3), c(1,1,NA,2,2))$chsq.statistic, wrn)
  expect_true(is.infinite(ksq))
})

################################################################################
#################################### IEVORA ####################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '2 columns dropped due to missing group information'
  expect_warning(res <- ievora(1:10, c(1,1,1,1,NA,NA,2,2,2,2)), wrn)
  expect_equal(res$obs.0+res$obs.1, 8)
  expect_warning(res <- ievora(1:10, c(1,1,1,1,NaN,NaN,2,2,2,2)), wrn)
  expect_equal(res$obs.0+res$obs.1, 8)
  expect_warning(res <- ievora(matrix(1:20, nrow=2), c(1,1,1,1,NA,NaN,2,2,2,2)), wrn)
  expect_equal(res$obs.0+res$obs.1, c(8,8))
})

