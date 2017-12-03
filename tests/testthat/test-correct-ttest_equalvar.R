context("correctness of ttest_equalvar")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_ttest_equalvar <- function(mat1, mat2, alt="two.sided", mu=0, conf=0.95) {
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
    res <- t.test(vec1, vec2, alternative=alt[i], mu=mu[i], conf.level=conf[i],
                  var.equal=TRUE
                  )

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
    # pooled variance
    vp <- rep(0, nrow(mat1))
    vp <- ifelse(nx > 1, vp + (nx-1) * vx, vp)
    vp <- ifelse(ny > 1, vp + (ny-1) * vy, vp)
    vp <- vp/(nx+ny-2)
    se[i]  <- sqrt(vp[i] * (1/nx[i] + 1/ny[i]))
  }

  data.frame(obs.x=nx, obs.y=ny, obs.tot=nt, mean.x=mx, mean.y=my, mean.diff=md,
             var.x=vx, var.y=vy, var.pooled=vp, stderr=se, df=df, stat.t=tst,
             pval=p, conf.low=cl, conf.high=ch, alternative=al, mean.null=m0,
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
  X[sample(length(X), nrow(X))] <- NA
  Y[sample(length(Y), nrow(X))] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  mus  <- seq(-1, 1, length.out=nrow(X))
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- base_ttest_equalvar(X, Y, alts, mus, cfs)
  t2 <- row.t.equalvar(X, Y, alts, mus, cfs)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  vals1 <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
  vals2 <- c(100000000000003, 100000000000002, 100000000000003, 100000000000000)
  t1 <- base_ttest_equalvar(vals1, vals2)
  t2 <- row.t.equalvar(vals1, vals2)
  expect_equal(t1, t2)

  # small numbers
  vals1 <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0)
  vals2 <- c(0.00000000000003, 0.00000000000002, 0.00000000000003, 0)
  t1 <- base_ttest_equalvar(vals1, vals2)
  t2 <- row.t.equalvar(vals1, vals2)
  expect_equal(t1, t2)
})


test_that("minumum allowed sample sizes give equal results", {
  # two numbers in one group and 1 in other group
  x <- matrix(rnorm(6), ncol=2); y <- matrix(rnorm(3), ncol=1)
  alt <- c("two.sided", "greater", "less")
  t1 <- base_ttest_equalvar(x, y, alt)
  t2 <- row.t.equalvar(x, y, alt)
  expect_equal(t1, t2)

  # two numbers in one group and 1 in other group with NAs
  x <- matrix(rnorm(6), ncol=2); y <- matrix(rnorm(6), ncol=2); x[,2] <- NA
  alt <- c("two.sided", "greater", "less")
  t1 <- base_ttest_equalvar(x, y, alt)
  t2 <- row.t.equalvar(x, y, alt)
  expect_equal(t1, t2)
})


test_that("parameter edge cases give equal results", {
  set.seed(14)
  alt <- c("l", "t", "g")
  mus <- c(-Inf, -1, 0, 1, Inf)
  cfs <- c(0, 0.5, 1)
  pars <- expand.grid(alt, mus, cfs, stringsAsFactors=FALSE)
  X1 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X2 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X1[sample(length(X1), nrow(pars)*2)] <- NA
  X2[sample(length(X2), nrow(pars)*2)] <- NA

  t1 <- base_ttest_equalvar(X1, X2, pars[,1], pars[,2], pars[,3])
  t2 <- row.t.equalvar(X1, X2, pars[,1], pars[,2], pars[,3])
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warnign when a row has less than 3 observations", {
  wrn <- '1 of the rows had less than 3 total observations\\. First occurrence at row 1'
  nacolumns <- c("stat.t", "pval", "conf.low", "conf.high")

  # no observations
  expect_warning(res <- row.t.equalvar(NA_integer_, NA_integer_), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 0)

  # one observation in each group
  expect_warning(res <- row.t.equalvar(1, 1), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 2)

  # one group has 2 observations and another 0
  expect_warning(res <- row.t.equalvar(c(1,2), NA_integer_), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 2)
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 0)

  # both have 2 observations but one is NA in each
  expect_warning(res <- row.t.equalvar(c(1,NA), c(3,NA)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 2)
})


test_that("warning when one group has no observations", {
  wrnX <- '1 of the rows had zero "x" observations\\. First occurrence at row 1'
  wrnY <- '1 of the rows had zero "y" observations\\. First occurrence at row 1'
  nacolumns <- c("stat.t", "pval", "conf.low", "conf.high")

  # no observations in X
  expect_warning(res <- row.t.equalvar(NA_integer_, c(1,1,1)), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 3)
  expect_equal(res$obs.y, 3)

  # no observations in Y
  expect_warning(res <- row.t.equalvar(c(1,1,1), NA_integer_), wrnY, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.tot, 3)
  expect_equal(res$obs.x, 3)
})


test_that("warning when row has constant values", {
  wrn <- '1 of the rows were essentially constant\\. First occurrence at row 1'
  nacolumns <- c("stat.t", "pval", "conf.low", "conf.high")

  # all values are equal
  expect_warning(res <- row.t.equalvar(c(1,1,1), c(1,1,1)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # values equal in their group
  expect_warning(res <- row.t.equalvar(c(1), c(2,2)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # equal values in each group + some NAs
  expect_warning(res <- row.t.equalvar(c(1,1,1,NA), c(NA,2,2,2)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
})


