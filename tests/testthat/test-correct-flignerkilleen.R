context("correctness of flignerkilleen")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_fligner <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad    <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ng <- nt <- ks <- p <- df <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])

    res <- fligner.test(vec, factor(grp))

    ng[i] <- length(unique(grp))
    nt[i] <- length(vec)
    ks[i] <- res$statistic
    p[i]  <- res$p.value
    df[i] <- res$parameter
  }

  data.frame(obs.tot=nt, obs.groups=ng, df=df, statistic=ks, pvalue=p,
             stringsAsFactors=FALSE
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

  t1 <- base_fligner(X, groups)
  t2 <- row_flignerkilleen(X, groups)

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
  t1 <- base_fligner(x, g)
  t2 <- row_flignerkilleen(x, g)
  expect_equal(t1, t2)

  # small numbers
  x <- cbind(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
             1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
             )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- base_fligner(x, g)
  t2 <- row_flignerkilleen(x, g)
  expect_equal(t1, t2)
})


test_that("constant values give equal results", {
  # all values are constant
  x <- c(1,1,1,1); g <- c("a","a","b","b")
  t1 <- base_fligner(x, g)
  t2 <- suppressWarnings(row_flignerkilleen(x, g))
  expect_equal(t1, t2)

  # within group values are constant
  x <- c(1,1,2,2); g <- c("a","a","b","b")
  t1 <- base_fligner(x, g)
  t2 <- suppressWarnings(row_flignerkilleen(x, g))
  expect_equal(t1, t2)

  # one group's values are constant
  x <- c(1,1,2,3); g <- c("a","a","b","b")
  t1 <- base_fligner(x, g)
  t2 <- suppressWarnings(row_flignerkilleen(x, g))
  expect_equal(t1, t2)
})


test_that("minimum allowed sample sizes give equal results", {
  # 3 samples 2 groups
  x <- rnorm(3); g <- c("a", "a", "b")
  t1 <- base_fligner(x, g)
  t2 <- row_flignerkilleen(x, g)
  expect_equal(t1, t2)

  # 3 samples 2 groups with NAs
  x <- rnorm(6); x[c(3:5)] <- NA
  g <- c("a","a","a","b","b","b")
  t1 <- base_fligner(x, g)
  t2 <- row_flignerkilleen(x, g)
  expect_equal(t1, t2)
})


################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning when rows have less than 2 groups", {
  wrn <- 'row_flignerkilleen: 1 of the rows had less than 2 groups with enough observations\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # zero groups - all are NA
  expect_warning(res <- row_flignerkilleen(1:10, rep(NA,10)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 0)

  # one gorup - others are NA
  expect_warning(res <- row_flignerkilleen(1:10, c(1, rep(NA,9))), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)

  # one group with non-NA elements
  expect_warning(res <- row_flignerkilleen(c(1,2,NA), c("a","a","b")), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)

  # all observations in one group
  expect_warning(res <- row_flignerkilleen(1:10, rep(1,10)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.groups, 1)
})


test_that("warning when all groups have only one element", {
  wrn <- 'row_flignerkilleen: 1 of the rows had one observation per group\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # 5 elements with 5 groups
  expect_warning(res <- row_flignerkilleen(c(1:5), c(1:5)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # with some NA groups
  expect_warning(res <- row_flignerkilleen(c(1:7), c(NA, 1:5, NA)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # with some NA elements that force groups to contain one element
  expect_warning(res <- row_flignerkilleen(c(NA,1:5,NA), c(1, 1:5, 5)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # when both elements and groups have NA values
  expect_warning(res <- row_flignerkilleen(c(NA,1:4,NA), c(1, 1, NA, NA, 2, 2)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
})


test_that("warning when none of the groups have variance", {
  wrn <- 'row_flignerkilleen: 1 of the rows had zero variance in all of the groups\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue")

  # all values are constant
  expect_warning(res <- row_flignerkilleen(c(1,1,1,1,1,1), c(1,1,2,2,3,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # all values within one group are constant
  expect_warning(res <- row_flignerkilleen(c(1,1,2,2,1.5,1.5), c(1,1,2,2,3,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # all values within one group are constant with NAs
  expect_warning(res <- row_flignerkilleen(c(1,1,NA,3,1,1,1,2,1,NA), c(1,1,1,NA,1,2,2,NA,2,2)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # all values are constant except a few for which the group is NA
  expect_warning(res <- row_flignerkilleen(c(1,1,3,1,1,1,2,1), c(1,1,NA,1,2,2,NA,2)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
})

