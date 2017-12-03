context("correctness of ttest_welch")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_ttest_welch <- function(mat1, mat2, alt="two.sided", mu=0, conf=0.95) {
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  if(length(alt)==1) alt <- rep(alt, nrow(mat1))
  if(length(mu)==1) mu <- rep(mu, nrow(mat1))
  if(length(conf)==1) conf <- rep(conf, nrow(mat1))

  mx <- my <- md <- vx <- vy <- vp <- nx <- ny <- nt <- tst <- p <- cl <- ch <-
    se <- df <- m0 <- cnf <- numeric(nrow(mat1))
  al <- character(nrow(mat1))
  for(i in 1:nrow(mat1)) {
    vec1 <- na.omit(mat1[i,])
    vec2 <- na.omit(mat2[i,])
    res <- t.test(vec1, vec2, alternative=alt[i], mu=mu[i], conf.level=conf[i])

    vx[i]  <- var(vec1)
    vy[i]  <- var(vec2)
    nx[i]  <- length(vec1)
    ny[i]  <- length(vec2)
    nt[i]  <- nx[i] + ny[i]
    mx[i]  <- res$estimate[1]
    my[i]  <- res$estimate[2]
    md[i]  <- mx[i]-my[i]
    tst[i] <- res$statistic
    p[i]   <- res$p.value
    cl[i]  <- res$conf.int[1]
    ch[i]  <- res$conf.int[2]
    df[i]  <- res$parameter
    m0[i]  <- res$null.value
    al[i]  <- res$alternative
    cnf[i] <- attr(res$conf.int, "conf.level")
    se[i]  <- sqrt(sqrt(vx[i]/nx[i])^2 + sqrt(vy[i]/ny[i])^2)
  }

  data.frame(obs.x=nx, obs.y=ny, obs.tot=nt, mean.x=mx, mean.y=my, mean.diff=md,
             var.x=vx, var.y=vy, stderr=se, df=df, statistic.t=tst, p.value=p,
             conf.low=cl, conf.high=ch, alternative=al, mean.null=m0,
             conf.level=cnf, stringsAsFactors=FALSE
             )
}

################################################################################
########################### TEST ON A RANDOM SAMPLE ############################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(10000), ncol=10)
  Y <- matrix(rnorm(10000), ncol=10)
  X[sample(length(X), 100)] <- NA
  Y[sample(length(Y), 100)] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  mus  <- seq(-1, 1, length.out=nrow(X))
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- base_ttest_welch(X, Y, alts, mus, cfs)
  t2 <- row.t.welch(X, Y, alts, mus, cfs)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  vals1 <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
  vals2 <- c(100000000000003, 100000000000002, 100000000000003, 100000000000000)
  t1 <- base_ttest_welch(vals1, vals2)
  t2 <- row.t.welch(vals1, vals2)
  expect_equal(t1, t2)

  # small numbers
  vals1 <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0)
  vals2 <- c(0.00000000000003, 0.00000000000002, 0.00000000000003, 0)
  t1 <- base_ttest_welch(vals1, vals2)
  t2 <- row.t.welch(vals1, vals2)
  expect_equal(t1, t2)
})


test_that("minumum allowed sample sizes give equal results", {
  # two observations in both groups
  x <- matrix(rnorm(6), ncol=2); y <- matrix(rnorm(6), ncol=2)
  alt <- c("two.sided", "greater", "less")
  t1 <- base_ttest_welch(x, y, alt)
  t2 <- row.t.welch(x, y, alt)
  expect_equal(t1, t2)

  # three observations in both groups - one is NA
  x <- matrix(rnorm(9), ncol=3); x[,1] <- NA
  y <- matrix(rnorm(9), ncol=3); y[,2] <- NA
  alt <- c("two.sided", "greater", "less")
  t1 <- base_ttest_welch(x, y, alt)
  t2 <- row.t.welch(x, y, alt)
  expect_equal(t1, t2)
})


test_that("parameter edge cases give equal results", {
  set.seed(14)
  alt <- c("l", "t", "g")
  mus <- c(-Inf, -1, 0, 1, Inf)
  cfs <- c(0, 0.95, 1)
  pars <- expand.grid(alt, mus, cfs, stringsAsFactors=FALSE)
  X1 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X2 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X1[sample(length(X1), nrow(pars)*2)] <- NA
  X2[sample(length(X2), nrow(pars)*2)] <- NA

  t1 <- base_ttest_welch(X1, X2, pars[,1], pars[,2], pars[,3])
  t2 <- row.t.welch(X1, X2, pars[,1], pars[,2], pars[,3])
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning when a row has less than 2 observations in one group", {
  wrnX <- '1 of the rows had less than 2 "x" observations\\. First occurrence at row 1'
  wrnY <- '1 of the rows had less than 2 "y" observations\\. First occurrence at row 1'
  nacolumns <- c("statistic.t", "p.value", "conf.low", "conf.high")

  # 0 observations in both
  expect_warning(res <- row.t.welch(NA_integer_, NA_integer_), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 0)
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 0)

  # 0 observations in X
  expect_warning(res <- row.t.welch(NA_integer_, 1), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 1)
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 1)

  # 0 observations in Y
  expect_warning(res <- row.t.welch(1, NA_integer_), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 1)
  expect_equal(res$obs.x, 1)
  expect_equal(res$obs.y, 0)

  # 2 observations in X but one is NA
  expect_warning(res <- row.t.welch(c(2,NA), c(1,1)), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 3)
  expect_equal(res$obs.x, 1)
  expect_equal(res$obs.y, 2)
})


test_that("warning when row has constant values", {
  wrn <- '1 of the rows were essentially constant\\. First occurrence at row 1'
  nacolumns <- c("statistic.t", "p.value", "conf.low", "conf.high")

  # all values are equal
  expect_warning(res <- row.t.welch(c(1,1,1), c(1,1,1)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # values equal in their group
  expect_warning(res <- row.t.welch(c(1,1), c(2,2)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # equal values in each group + some NAs
  expect_warning(res <- row.t.welch(c(1,1,1,NA), c(NA,2,2,2)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
})


