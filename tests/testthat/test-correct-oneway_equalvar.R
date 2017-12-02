context("correctness of oneway_equalvar")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_oneway_equalvar <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  st <- sr <- mt <- mr <- ot <- og <- dft <- dfr <- fst <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- mat[i,!bad]
    grp <- factor(groups[!bad])
    res <- summary(aov(vec ~ grp))[[1]]

    st[i]  <- res[1,2]
    sr[i]  <- res[2,2]
    mt[i]  <- res[1,3]
    mr[i]  <- res[2,3]
    dft[i] <- res[1,1]
    dfr[i] <- res[2,1]
    fst[i] <- res[1,4]
    p[i]   <- res[1,5]
    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
  }

  data.frame(obs.tot=ot, obs.groups=og, sum.sq.treatment=st, sum.sq.residuals=sr,
             mean.sq.treatment=mt, mean.sq.residuals=mr, df.treatment=dft,
             df.residuals=dfr, statistic.F=fst, p.value=p
             )
}

# NOTE: fall-back version for when the first one has rounding errors.
base_oneway_equalvar2 <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  st <- sr <- mt <- mr <- ot <- og <- dft <- dfr <- fst <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- mat[i,!bad]
    grp <- factor(groups[!bad])
    res <- oneway.test(vec ~ grp, var.equal=TRUE)

    dft[i] <- res$parameter[1]
    dfr[i] <- res$parameter[2]
    fst[i] <- res$statistic
    p[i]   <- res$p.value
    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
    st[i]  <- sum(tapply(vec, grp, length) * (tapply(vec, grp, mean)-mean(vec))^2)
    sr[i]  <- sum((tapply(vec, grp, length)-1) * tapply(vec, grp, var))
    mt[i]  <- st[i] / dft[i]
    mr[i]  <- sr[i] / dfr[i]
  }

  data.frame(obs.tot=ot, obs.groups=og, sum.sq.treatment=st, sum.sq.residuals=sr,
             mean.sq.treatment=mt, mean.sq.residuals=mr, df.treatment=dft,
             df.residuals=dfr, statistic.F=fst, p.value=p
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

  t1 <- base_oneway_equalvar(X, groups)
  t2 <- suppressWarnings(row.oneway.equalvar(X, groups))

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # NOTE: using fall-back version
  # big numbers
  x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
         100000000000003, 100000000000002, 100000000000003, 100000000000000
         )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- base_oneway_equalvar2(x, g)
  t2 <- row.oneway.equalvar(x, g)
  expect_equal(t1, t2)

  # small numbers
  x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
         1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
         )
  g <- c(rep("a", 4), rep("b", 4))
  t1 <- base_oneway_equalvar2(x, g)
  t2 <- row.oneway.equalvar(x, g)
  expect_equal(t1, t2)
})


test_that("constant values give equal results", {
  # all values are constant
  x <- c(1,1,1,1); g <- c("a","a","b","b")
  t1 <- base_oneway_equalvar(x, g)
  t2 <- suppressWarnings(row.oneway.equalvar(x, g))
  expect_equal(t1, t2)

  # within group values are constant
  # NOTE: using fall-back version
  x <- c(1,1,2,2); g <- c("a","a","b","b")
  t1 <- base_oneway_equalvar2(x, g)
  t2 <- suppressWarnings(row.oneway.equalvar(x, g))
  expect_equal(t1, t2)

  # one group's values are constant
  x <- c(1,1,2,3); g <- c("a","a","b","b")
  t1 <- base_oneway_equalvar(x, g)
  t2 <- suppressWarnings(row.oneway.equalvar(x, g))
  expect_equal(t1, t2)
})


test_that("minimal allowed sample size gives equal results", {
  # two groups one with two values and another with one value
  x <- c(1,2,3)
  g <- c("a","a","b")
  t1 <- base_oneway_equalvar(x, g)
  t2 <- suppressWarnings(row.oneway.equalvar(x, g))
  expect_equal(t1, t2)

  # with NAs
  x <- c(1,2,3,NA,4)
  g <- c("a","a","b","b",NA)
  t1 <- base_oneway_equalvar(x, g)
  t2 <- suppressWarnings(row.oneway.equalvar(x, g))
  expect_equal(t1, t2)
})


test_that("groups with one remaining member give equal results", {
  # two groups - one has one remaining element
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[1:2] <- NA
  t1 <- base_oneway_equalvar(x, g)
  t2 <- suppressWarnings(row.oneway.equalvar(x, g))
  expect_equal(t1, t2)

  # many groups - all groups except one have one remaining element
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3); g[c(1,4:5,7:8,10:11)] <- NA
  t1 <- base_oneway_equalvar(x, g)
  t2 <- suppressWarnings(row.oneway.equalvar(x, g))
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning is shown when columns are removed because of NA groups", {
  wrn <- '2 columns dropped due to missing group information'

  # 2 NAs
  expect_warning(res <- row.oneway.equalvar(1:10, c(1,1,1,1,NA,NA,2,2,2,2)), wrn, all=TRUE)
  expect_equal(res$obs.tot, 8)
  expect_equal(res$obs.groups, 2)

  # 4 groups with one group dropped because of missing x values
  x <- rnorm(10); x[c(1,2)] <- NA
  g <- c(1,1,2,2,NA,NA,3,3,4,4)
  expect_warning(res <- row.oneway.equalvar(x, g), wrn)
  expect_equal(res$obs.tot, 6)
  expect_equal(res$obs.groups, 3)
})


test_that("warning when a rows has less than 2 groups", {
  wrn <- '1 of the rows had less than 2 groups with enough observations\\. First occurrence at row 1'
  nacolumns <- c("statistic.F", "p.value")

  # one group
  x <- 1:10; g <- rep(1, 10)
  expect_warning(res <- row.oneway.equalvar(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 1)

  # two groups but one only has NAs
  x <- 1:10; x[6:10] <- NA
  g <- c(rep(1,5), rep(2,5))
  expect_warning(res <- row.oneway.equalvar(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 5)
  expect_equal(res$obs.groups, 1)
})


test_that("warning when a row has 1 observation per group", {
  wrn <- '1 of the rows had one observation per group\\. First occurrence at row 1'
  nacolumns <- c("statistic.F", "p.value")

  # 10 groups 10 observations
  x <- 1:10; g <- 1:10
  expect_warning(res <- row.oneway.equalvar(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 10)
  expect_equal(res$obs.groups, 10)

  # two groups 3 observations but one is NA
  x <- rnorm(3); x[3] <- NA
  g <- c("a", "b", "b")
  expect_warning(res <- row.oneway.equalvar(x, g), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 2)
  expect_equal(res$obs.groups, 2)
})

test_that("warning when all values are constant", {
  wrn <- '1 of the rows had essentially constant values\\. First occurrence at row 1'
  nacolumns <- c("statistic.F", "p.value")

  # two groups - all constant values
  x <- c(0,0,0,0); g <- c(1,1,2,2)
  expect_warning(res <- row.oneway.equalvar(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 4)
  expect_equal(res$obs.groups, 2)

  # three groups - constant values + NAs
  x <- c(3,3,NA,3,3,NA,3,3,NA); g <- c(1,1,1,2,2,2,3,3,3)
  expect_warning(res <- row.oneway.equalvar(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 6)
  expect_equal(res$obs.groups, 3)
})



test_that("warning when all values within each group are constant", {
  wrn <- '1 of the rows had zero within group variance: result might be unreliable\\. First occurrence at row 1'

  # two groups - constant values within group
  x <- c(1,1,0,0); g <- c(1,1,2,2)
  expect_warning(res <- row.oneway.equalvar(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 4)
  expect_equal(res$obs.groups, 2)

  # three groups - constant values within each group + NAs
  x <- c(3,3,NA,0,0,NA,-1,-1,NA); g <- c(1,1,1,2,2,2,3,3,3)
  expect_warning(res <- row.oneway.equalvar(x, g), wrn, all=TRUE)
  expect_equal(res$obs.tot, 6)
  expect_equal(res$obs.groups, 3)
})


