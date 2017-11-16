context("correctness of ttest_paired")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_ttest_paired <- function(mat1, mat2, alt="two.sided", mu=0, conf=0.95) {
  stopifnot(ncol(mat1)==ncol(mat2))
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  if(length(alt)==1) alt <- rep(alt, nrow(mat1))
  if(length(mu)==1) mu <- rep(mu, nrow(mat1))
  if(length(conf)==1) conf <- rep(conf, nrow(mat1))

  mx <- my <- md <- vx <- vy <- vd <- nx <- ny <- nt <- tst <- p <- cl <- ch <-
    se <- df <- m0 <- cnf <- numeric(nrow(mat1))
  al <- character(nrow(mat1))
  for(i in 1:nrow(mat1)) {
    vec1 <- mat1[i,]
    vec2 <- mat2[i,]
    res <- t.test(vec1, vec2, alternative=alt[i], mu=mu[i], conf.level=conf[i],
                  paired=TRUE
                  )

    mx[i]  <- mean(na.omit(vec1))
    my[i]  <- mean(na.omit(vec2))
    md[i]  <- res$estimate
    vx[i]  <- var(na.omit(vec1))
    vy[i]  <- var(na.omit(vec2))
    vd[i]  <- var(na.omit(vec1-vec2))
    nx[i]  <- length(na.omit(vec1))
    ny[i]  <- length(na.omit(vec2))
    nt[i]  <- length(na.omit(vec1-vec2))
    tst[i] <- res$statistic
    p[i]   <- res$p.value
    cl[i]  <- res$conf.int[1]
    ch[i]  <- res$conf.int[2]
    df[i]  <- res$parameter
    m0[i]  <- res$null.value
    al[i]  <- res$alternative
    cnf[i] <- attr(res$conf.int, "conf.level")
    se[i]  <- sqrt(vd[i]/nt[i])
  }

  data.frame(mean.x=mx, mean.y=my, mean.diff=md, var.x=vx, var.y=vy,
             var.diff=vd, obs.x=nx, obs.y=ny, obs.pair=nt, t.statistic=tst,
             p.value=p, ci.low=cl, ci.high=ch, stderr=se, df=df, mean.null=m0,
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

  t1 <- base_ttest_paired(X, Y, alts, mus, cfs)
  t2 <- ttest_paired(X, Y, alts, mus, cfs)

  expect_equal(t1, t2)
})

################################## EDGE CASES ##################################

test_that("weird numbers give equal results", {
  vals1 <- c(1.00000000000001, 1.00000000000002, 1.00000000000003)
  vals2 <- c(1.00000000000003, 1.00000000000002, 1.00000000000003)
  expect_equal(base_ttest_paired(vals1, vals2), ttest_paired(vals1, vals2))
})

test_that("minumum allowed sample sizes give equal results", {
  expect_equal(base_ttest_paired(1:2, 3:2), ttest_paired(1:2, 3:2))
  expect_equal(base_ttest_paired(c(1:2,NA,3), c(3:2,3,NA)), ttest_paired(c(1:2,NA,3), c(3:2,3,NA)))
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

  t1 <- base_ttest_paired(X1, X2, pars[,1], pars[,2], pars[,3])
  t2 <- ttest_paired(X1, X2, pars[,1], pars[,2], pars[,3])
  expect_equal(t1, t2)
})

