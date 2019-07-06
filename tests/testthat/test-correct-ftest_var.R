context("correctness of f_var")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_f_var <- function(mat1, mat2, rat=1, alt="two.sided", conf=0.95) {
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  if(length(rat)==1) rat <- rep(rat, nrow(mat1))
  if(length(alt)==1) alt <- rep(alt, nrow(mat1))
  if(length(conf)==1) conf <- rep(conf, nrow(mat1))

  nx <- ny <- nt <- vx <- vy <- vrat <- df1 <- df2 <- stat <- p <- cl <- ch <-
    r0 <- cnf <- numeric(nrow(mat1))
  al <- character(nrow(mat1))
  for(i in 1:nrow(mat1)) {
    vec1 <- mat1[i,]
    vec2 <- mat2[i,]
    res <- var.test(vec1, vec2, ratio=rat[i], alternative=alt[i], conf.level=conf[i])

    nx[i]  <- length(na.omit(vec1))
    ny[i]  <- length(na.omit(vec2))
    nt[i]  <- length(na.omit(vec1)) + length(na.omit(vec2))
    vx[i]  <- var(vec1, na.rm=TRUE)
    vy[i]  <- var(vec2, na.rm=TRUE)
    vrat[i]  <- res$estimate
    df1[i]  <- res$parameter[1]
    df2[i]  <- res$parameter[2]
    stat[i]  <- res$statistic
    p[i] <- res$p.value
    cl[i]  <- res$conf.int[1]
    ch[i]  <- res$conf.int[2]
    r0[i]  <- res$null.value
    al[i]  <- res$alternative
    cnf[i] <- attr(res$conf.int, "conf.level")
  }

  data.frame(obs.x=nx, obs.y=ny, obs.tot=nt, var.x=vx, var.y=vy, var.ratio=vrat,
             df.num=df1, df.denom=df2, statistic=stat, pvalue=p, conf.low=cl,
             conf.high=ch, ratio.null=r0, alternative=al, conf.level=cnf,
             stringsAsFactors=FALSE
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
  rats <- seq(0.001, 1000, length.out=nrow(X))
  cfs  <- seq(0.001, 0.999, length.out=nrow(X))

  t1 <- base_f_var(X, Y, rats, alts, cfs)
  t2 <- row_f_var(X, Y, rats, alts, cfs)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  vals1 <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
  vals2 <- c(100000000000003, 100000000000002, 100000000000003, 100000000000000)
  t1 <- base_f_var(vals1, vals2)
  t2 <- row_f_var(vals1, vals2)
  expect_equal(t1, t2)

  # small numbers
  vals1 <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0)
  vals2 <- c(0.00000000000003, 0.00000000000002, 0.00000000000003, 0)
  t1 <- base_f_var(vals1, vals2)
  t2 <- row_f_var(vals1, vals2)
  expect_equal(t1, t2)
})


test_that("minumum allowed sample sizes give equal results", {
  # two numbers in both groups
  x <- matrix(rnorm(6), ncol=2); y <- matrix(rnorm(6), ncol=2)
  alt <- c("two.sided", "greater", "less")
  t1 <- base_f_var(x, y, alt=alt)
  t2 <- row_f_var(x, y, alternative=alt)
  expect_equal(t1, t2)

  # two numbers in both groups + one NA
  x <- matrix(rnorm(9), ncol=3); y <- matrix(rnorm(6), ncol=2); x[,2] <- NA
  alt <- c("two.sided", "greater", "less")
  t1 <- base_f_var(x, y, alt=alt)
  t2 <- row_f_var(x, y, alt=alt)
  expect_equal(t1, t2)
})

test_that("parameter edge cases give equal results", {
  set.seed(14)
  alt <- c("l", "t", "g")
  rat <- c(0.000001, 0.5, 1, 2, 99999)
  cfs <- c(0.000001, 0.5, 0.9999)
  pars <- expand.grid(rat, alt, cfs, stringsAsFactors=FALSE)
  X1 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X2 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X1[sample(length(X1), nrow(pars)*2)] <- NA
  X2[sample(length(X2), nrow(pars)*2)] <- NA

  t1 <- base_f_var(X1, X2, pars[,1], pars[,2], pars[,3])
  t2 <- row_f_var(X1, X2, pars[,1], pars[,2], pars[,3])
  expect_equal(t1, t2)
})

test_that("parameter edge cases give equal results (when x is constant)", {
  set.seed(14)
  alt <- c("l", "t", "g")
  rat <- c(0.000001, 0.5, 1, 2, 99999)
  cfs <- c(0.000001, 0.5, 0.9999)
  pars <- expand.grid(rat, alt, cfs, stringsAsFactors=FALSE)
  X1 <- matrix(1, nrow=nrow(pars), ncol=10)
  X2 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X1[sample(length(X1), nrow(pars)*2)] <- NA
  X2[sample(length(X2), nrow(pars)*2)] <- NA

  t1 <- base_f_var(X1, X2, pars[,1], pars[,2], pars[,3])
  t2 <- suppressWarnings(row_f_var(X1, X2, pars[,1], pars[,2], pars[,3]))
  expect_equal(t1, t2)
  expect_true(all(t2$var.x==0))
})

test_that("parameter edge cases give equal results (when y is constant)", {
  set.seed(14)
  alt <- c("l", "t", "g")
  rat <- c(0.000001, 0.5, 1, 2, 99999)
  cfs <- c(0.000001, 0.5, 0.9999)
  pars <- expand.grid(rat, alt, cfs, stringsAsFactors=FALSE)
  X1 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X2 <- matrix(1, nrow=nrow(pars), ncol=10)
  X1[sample(length(X1), nrow(pars)*2)] <- NA
  X2[sample(length(X2), nrow(pars)*2)] <- NA

  t1 <- base_f_var(X1, X2, pars[,1], pars[,2], pars[,3])
  t2 <- suppressWarnings(row_f_var(X1, X2, pars[,1], pars[,2], pars[,3]))
  expect_equal(t1, t2)
  expect_true(all(t2$var.y==0))
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning when a infinite values are removed", {
  wrnX <- 'row_f_var: 1 of the rows had infinite "x" observations that were removed.'
  wrnY <- 'row_f_var: 1 of the rows had infinite "y" observations that were removed.'

  # -Inf and Inf among x observations
  expect_warning(res <- row_f_var(c(-Inf,5,7,Inf), 1:4), wrnX, all=TRUE)
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 4)
  expect_equal(res$obs.tot, 6)

  # -Inf and Inf among y observations
  expect_warning(res <- row_f_var(1:4, c(-Inf,5,7,Inf)), wrnY, all=TRUE)
  expect_equal(res$obs.x, 4)
  expect_equal(res$obs.y, 2)
  expect_equal(res$obs.tot, 6)

  # -Inf and Inf among both
  expect_warning(res <- row_f_var(c(1,2,3,Inf), c(-Inf,7,6,5)), wrnX)
  expect_warning(res <- row_f_var(c(1,2,3,Inf), c(-Inf,7,6,5)), wrnY)
  expect_equal(res$obs.x, 3)
  expect_equal(res$obs.y, 3)
  expect_equal(res$obs.tot, 6)
})


test_that("warning when some groups have less than 2 observations", {
  wrnX <- 'row_f_var: 1 of the rows had less than 1 "x" observation\\.\nFirst occurrence at row 1'
  wrnY <- 'row_f_var: 1 of the rows had less than 1 "y" observation\\.\nFirst occurrence at row 1'
  nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

  # x has 0 observations
  expect_warning(res <- row_f_var(numeric(), 1:2), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 2)
  expect_equal(res$obs.tot, 2)

  # y has 0 observations
  expect_warning(res <- row_f_var(1:2, numeric()), wrnY, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 0)
  expect_equal(res$obs.tot, 2)

  # x has 1 observations
  expect_warning(res <- row_f_var(3, 1:2), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 1)
  expect_equal(res$obs.y, 2)
  expect_equal(res$obs.tot, 3)

  # y has 1 observations
  expect_warning(res <- row_f_var(1:2, 3), wrnY, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 1)
  expect_equal(res$obs.tot, 3)

  # x has only NA
  expect_warning(res <- row_f_var(NA_integer_, 1:2), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 2)
  expect_equal(res$obs.tot, 2)

  # y has only NA
  expect_warning(res <- row_f_var(1:2, NA_integer_), wrnY, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 0)
  expect_equal(res$obs.tot, 2)

  # x has observation + NA
  expect_warning(res <- row_f_var(c(0,NA), 1:2), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 1)
  expect_equal(res$obs.y, 2)
  expect_equal(res$obs.tot, 3)

  # y has observation + NA
  expect_warning(res <- row_f_var(1:2, c(0,NA)), wrnY, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 1)
  expect_equal(res$obs.tot, 3)

  # x has observation + Infs and NAs
  expect_warning(res <- row_f_var(c(0,Inf,NA,-Inf), 1:2), wrnX)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 1)
  expect_equal(res$obs.y, 2)
  expect_equal(res$obs.tot, 3)

  # y has observation + Infs and NAs
  expect_warning(res <- row_f_var(1:2, c(0,Inf,NA,-Inf)), wrnY)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 1)
  expect_equal(res$obs.tot, 3)

})

test_that("warning when one of the groups have 0 variance", {
  wrnX <- 'row_f_var: 1 of the rows had zero variance in "x"\\.\nFirst occurrence at row 1'
  wrnY <- 'row_f_var: 1 of the rows had zero variance in "y"\\.\nFirst occurrence at row 1'

  # x has constant values
  expect_warning(res <- row_f_var(rep(1,5), rnorm(5)), wrnX, all=TRUE)
  expect_equal(res$var.x, 0)
  expect_equal(res$var.ratio, 0)

  # y has constant values
  expect_warning(res <- row_f_var(rnorm(5), rep(1,5)), wrnY, all=TRUE)
  expect_equal(res$var.y, 0)
  expect_equal(res$var.ratio, Inf)

  # x has constant values with Inf and NAs
  expect_warning(res <- row_f_var(c(rep(1,5), NA, Inf, -Inf), rnorm(5)), wrnX)
  expect_equal(res$var.x, 0)
  expect_equal(res$var.ratio, 0)

  # y has constant values with Inf and NAs
  expect_warning(res <- row_f_var(rnorm(5), c(rep(1,5), NA, Inf, -Inf)), wrnY)
  expect_equal(res$var.y, 0)
  expect_equal(res$var.ratio, Inf)

  # both x and y have constant values
  expect_warning(res <- row_f_var(rep(2,5), rep(1,5)), wrnX)
  expect_warning(res <- row_f_var(rep(2,5), rep(1,5)), wrnY)
  expect_equal(res$var.x, 0)
  expect_equal(res$var.y, 0)
  expect_equal(res$var.ratio, NaN)
})

