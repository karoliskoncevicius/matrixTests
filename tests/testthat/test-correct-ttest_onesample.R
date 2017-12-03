context("corectness of ttest_onesample")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_ttest_onesample <- function(mat, alt="two.sided", mu=0, conf=0.95) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)
  if(length(alt)==1) alt <- rep(alt, nrow(mat))
  if(length(mu)==1) mu <- rep(mu, nrow(mat))
  if(length(conf)==1) conf <- rep(conf, nrow(mat))

  mx <- vx <- nx <- tst <- p <- cl <- ch <- se <- df <- m0 <- cnf <- numeric(nrow(mat))
  al <- character(nrow(mat))
  for(i in 1:nrow(mat)) {
    vec <- na.omit(mat[i,])
    res <- t.test(vec, alternative=alt[i], mu=mu[i], conf.level=conf[i])

    vx[i]  <- var(vec)
    nx[i]  <- length(vec)
    mx[i]  <- res$estimate
    tst[i] <- res$statistic
    p[i]   <- res$p.value
    cl[i]  <- res$conf.int[1]
    ch[i]  <- res$conf.int[2]
    se[i]  <- sqrt(vx[i]) / sqrt(nx[i])
    df[i]  <- res$parameter
    m0[i]  <- res$null.value
    al[i]  <- res$alternative
    cnf[i] <- attr(res$conf.int, "conf.level")
  }

  data.frame(obs.x=nx, mean.x=mx, var.x=vx, stderr=se, df=df, statistic.t=tst,
             p.value=p, conf.low=cl, conf.high=ch, alternative=al, mean.null=m0,
             conf.level=cnf, stringsAsFactors=FALSE
             )
}

################################################################################
########################### TEST ON A RANDOM SAMPLE ############################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(10000), ncol=10)
  X[sample(length(X), nrow(X))] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  mus  <- runif(nrow(X), -1,1)
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- base_ttest_onesample(X, alts, mus, cfs)
  t2 <- row.t.onesample(X, alts, mus, cfs)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  vals <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
  t1 <- base_ttest_onesample(vals)
  t2 <- row.t.onesample(vals)
  expect_equal(t1, t2)

  # small numbers
  vals <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0)
  t1 <- base_ttest_onesample(vals)
  t2 <- row.t.onesample(vals)
  expect_equal(t1, t2)
})


test_that("minumum allowed sample sizes give equal results", {
  # two numbers
  x <- matrix(rnorm(6), ncol=2)
  alt <- c("two.sided", "greater", "less")
  t1 <- base_ttest_onesample(x, alt)
  t2 <- suppressWarnings(row.t.onesample(x, alt))
  expect_equal(t1, t2)

  # three numbers with NAs
  x <- matrix(rnorm(9), ncol=3); x[,3] <- NA
  alt <- c("two.sided", "greater", "less")
  t1 <- base_ttest_onesample(x, alt)
  t2 <- suppressWarnings(row.t.onesample(x, alt))
  expect_equal(t1, t2)
})


test_that("parameter edge cases give equal results", {
  set.seed(14)
  alt <- c("l", "t", "g")
  mus <- c(-Inf, -1, 0, 1, Inf)
  cfs <- c(0, 0.5, 1)
  pars <- expand.grid(alt, mus, cfs, stringsAsFactors=FALSE)
  X <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X[sample(length(X), nrow(pars)*2)] <- NA

  t1 <- base_ttest_onesample(X, pars[,1], pars[,2], pars[,3])
  t2 <- row.t.onesample(X, pars[,1], pars[,2], pars[,3])
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warnign when row has less than 2 available observations", {
  wrn <- '1 of the rows had less than 2 "x" observations\\. First occurrence at row 1'
  nacolumns <- c("statistic.t", "p.value", "conf.low", "conf.high")

  # single observation
  expect_warning(res <- row.t.onesample(1), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 1)

  # single observation with some NA values
  expect_warning(res <- row.t.onesample(c(0,NA,NA)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 1)

  # zero observations
  expect_warning(res <- row.t.onesample(NA_integer_), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 0)
})


test_that("warning when a row has all constant values", {
  wrn <- '1 of the rows were essentially constant\\. First occurrence at row 1'
  nacolumns <- c("statistic.t", "p.value", "conf.low", "conf.high")

  # two equal observations
  expect_warning(res <- row.t.onesample(c(1,1)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 2)
  expect_equal(res$var.x, 0)

  # three observations with some NA values
  expect_warning(res <- row.t.onesample(c(0,0,0,NA,NA)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 3)
  expect_equal(res$var.x, 0)
})

