context("corectness of ttest_onegroup")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_ttest_onegroup <- function(mat, alt="two.sided", mu=0, conf=0.95) {
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

  data.frame(mean.x=mx, var.x=vx, obs.x=nx, t.statistic=tst, p.value=p,
             ci.low=cl, ci.high=ch, stderr=se, df=df,
             mean.null=m0, conf.level=cnf, alternative=al,
             stringsAsFactors=FALSE
             )
}

################################################################################
######################### TEST CONSISTENCY WITH T.TEST #########################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(1000), ncol=10)
  X[sample(length(X), 100)] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  mus  <- runif(nrow(X), -1,1)
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- base_ttest_onegroup(X, alts, mus, cfs)
  t2 <- ttest_onegroup(X, alts, mus, cfs)

  expect_equal(t1, t2)
})

################################## EDGE CASES ##################################

test_that("weird numbers give equal results", {
  vals <- c(1.00000000000001, 1.00000000000002, 1.00000000000003)
  expect_equal(base_ttest_onegroup(vals), ttest_onegroup(vals))
  expect_equal(base_ttest_onegroup(1:10^5), ttest_onegroup(1:10^5))
})

test_that("minumum allowed sample sizes give equal results", {
  expect_equal(base_ttest_onegroup(1:2), ttest_onegroup(1:2))
  expect_equal(base_ttest_onegroup(c(1,2,NA,NA,NA)), ttest_onegroup(c(1,2,NA,NA,NA)))
})

test_that("parameter edge cases give equal results", {
  set.seed(14)
  alt <- c("l", "t", "g")
  mus <- c(-Inf, 0, Inf)
  cfs <- c(0, 0.95, 1)
  pars <- expand.grid(alt, mus, cfs, stringsAsFactors=FALSE)
  X <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X[sample(length(X), nrow(pars)*2)] <- NA

  t1 <- base_ttest_onegroup(X, pars[,1], pars[,2], pars[,3])
  t2 <- ttest_onegroup(X, pars[,1], pars[,2], pars[,3])
  expect_equal(t1, t2)
})

