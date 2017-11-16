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

  data.frame(mean.x=mx, mean.y=my, mean.diff=md, var.x=vx, var.y=vy,
             obs.x=nx, obs.y=ny, obs.tot=nt, t.statistic=tst, p.value=p,
             ci.low=cl, ci.high=ch, stderr=se, df=df, mean.null=m0,
             conf.level=cnf, alternative=al, stringsAsFactors=FALSE
             )
}

################################################################################
######################### TEST CONSISTENCY WITH T.TEST #########################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(1000), ncol=10)
  Y <- matrix(rnorm(1000), ncol=10)
  X[sample(length(X), 100)] <- NA
  Y[sample(length(Y), 100)] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  mus  <- seq(-1, 1, length.out=nrow(X))
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- base_ttest_welch(X, Y, alts, mus, cfs)
  t2 <- ttest_welch(X, Y, alts, mus, cfs)

  expect_equal(t1, t2)
})

################################## EDGE CASES ##################################

test_that("weird numbers give equal results", {
  vals1 <- c(1.00000000000001, 1.00000000000002, 1.00000000000003)
  vals2 <- c(1.00000000000003, 1.00000000000002, 1.00000000000003)
  expect_equal(base_ttest_welch(vals1, vals2), ttest_welch(vals1, vals2))
  expect_equal(base_ttest_welch(1:10^5, vals1), ttest_welch(1:10^5, vals1))
})

test_that("minumum allowed sample sizes give equal results", {
  expect_equal(base_ttest_welch(1:2, 2:3), ttest_welch(1:2, 2:3))
  expect_equal(base_ttest_welch(1:10, c(NA,1,2)), ttest_welch(1:10, c(NA,1,2)))
})

test_that("parameter edge cases give equal results", {
  set.seed(14)
  alt <- c("l", "t", "g")
  mus <- c(-Inf, 0, Inf)
  cfs <- c(0, 0.95, 1)
  pars <- expand.grid(alt, mus, cfs, stringsAsFactors=FALSE)
  X1 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X2 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X1[sample(length(X1), nrow(pars)*2)] <- NA
  X2[sample(length(X2), nrow(pars)*2)] <- NA

  t1 <- base_ttest_welch(X1, X2, pars[,1], pars[,2], pars[,3])
  t2 <- ttest_welch(X1, X2, pars[,1], pars[,2], pars[,3])
  expect_equal(t1, t2)
})

