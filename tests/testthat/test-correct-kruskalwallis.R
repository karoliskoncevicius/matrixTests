context("correctness of kruskalwallis")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_kruskal <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ot <- og <- df <- chs <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])
    res <- kruskal.test(vec, factor(grp))

    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
    df[i]  <- res$parameter
    chs[i] <- res$statistic
    p[i]   <- res$p.value
  }

  data.frame(obs.tot=ot, obs.groups=og, df, statistic.chsq=chs, p.value=p)
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

  t1 <- base_kruskal(X, factor(groups))
  t2 <- suppressWarnings(row.kruskalwallis(X, groups))

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
  t1 <- base_kruskal(x, g)
  t2 <- row.kruskalwallis(x, g)
  expect_equal(t1, t2)

  # small numbers
  x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
         1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
         )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- base_kruskal(x, g)
  t2 <- row.kruskalwallis(x, g)
  expect_equal(t1, t2)
})


test_that("constant values give equal results", {
  # all values are constant
  x <- c(1,1,1,1); g <- c("a","a","b","b")
  t1 <- base_kruskal(x, g)
  t2 <- suppressWarnings(row.kruskalwallis(x, g))
  expect_equal(t1, t2)

  # within group values are constant
  x <- c(1,1,2,2); g <- c("a","a","b","b")
  t1 <- base_kruskal(x, g)
  t2 <- suppressWarnings(row.kruskalwallis(x, g))
  expect_equal(t1, t2)

  # one group's values are constant
  x <- c(1,1,2,3); g <- c("a","a","b","b")
  t1 <- base_kruskal(x, g)
  t2 <- suppressWarnings(row.kruskalwallis(x, g))
  expect_equal(t1, t2)
})


test_that("minimal allowed sample size gives equal results", {
  # two groups with one value per group
  x <- c(1,2); g <- c("a","b")
  t1 <- base_kruskal(x, g)
  t2 <- row.kruskalwallis(x, g)
  expect_equal(t1, t2)

  # two groups with one value per group and NAs
  x <- c(1,2,NA,3)
  g <- c("a","b","a",NA)
  t1 <- base_kruskal(x, g)
  t2 <- suppressWarnings(row.kruskalwallis(x, g))
  expect_equal(t1, t2)
})


test_that("groups with one remaining member give equal results", {
  # one group has one remaining element
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[1:2] <- NA
  t1 <- base_kruskal(x, g)
  t2 <- suppressWarnings(row.kruskalwallis(x, g))
  expect_equal(t1, t2)

  # all groups have one remaining element
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[c(1:2,4:5,7:8,10:11)] <- NA
  t1 <- base_kruskal(x, g)
  t2 <- suppressWarnings(row.kruskalwallis(x, g))
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '2 columns dropped due to missing group information'

  # 2 NAs
  expect_warning(res <- row.kruskalwallis(1:10, c(1,1,1,1,NA,NA,2,2,2,2)), wrn, all=TRUE)
  expect_equal(res$obs.tot, 8)
  expect_equal(res$obs.groups, 2)

  # 4 groups with one group dropped because of missing x values
  x <- rnorm(10); x[c(1,2)] <- NA
  g <- c(1,1,2,2,NA,NA,3,3,4,4)
  expect_warning(res <- row.kruskalwallis(x, g), wrn)
  expect_equal(res$obs.tot, 6)
  expect_equal(res$obs.groups, 3)
})


test_that("warning when row has less than 2 complete observations", {
  wrn <- '1 of the rows had less than 2 total observations\\. First occurrence at row 1'
  nacolumns <- c("statistic.chsq", "p.value")

  # one value one group
  expect_warning(res <- row.kruskalwallis(1, "a"), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)
  expect_equal(res$obs.tot, 1)

  # one value one group with NAs
  expect_warning(res <- row.kruskalwallis(c(1,NA,NA,NA), c(1,1,2,2)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)
  expect_equal(res$obs.tot, 1)
})


test_that("warning when rows have less than 2 groups", {
  wrn <- '1 of the rows had less than 2 groups with enough observations\\. First occurrence at row 1'
  nacolumns <- c("statistic.chsq", "p.value")

  # single group with 10 observations
  expect_warning(res <- row.kruskalwallis(1:10, rep(1,10)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)
  expect_equal(res$obs.tot, 10)

  # 3 groups but others have only NA values
  x <- rnorm(10); x[9:10] <- NA
  g <- c(rep("a", 8), "b", "c")
  expect_warning(res <- row.kruskalwallis(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)
  expect_equal(res$obs.tot, 8)
})

test_that("warning when rows consist of single value", {
  wrn <- '1 of the rows were essentially constant\\. First occurrence at row 1'
  nacolumns <- c("statistic.chsq", "p.value")

  # two groups
  expect_warning(res <- row.kruskalwallis(rep(1, 3), c("a","a","b")), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 2)
  expect_equal(res$obs.tot, 3)

  # two groups with NA values
  x <- c(rep(1,8),2,NA,NA)
  g <- c(rep("a", 4), rep("b", 4), NA, "c", "d")
  expect_warning(res <- row.kruskalwallis(x, g), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 2)
  expect_equal(res$obs.tot, 8)
})

