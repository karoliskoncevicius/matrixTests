context("correctness of cor_pearson")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_cor_pearson <- function(mat1, mat2, alt="two.sided", conf=0.95) {
  stopifnot(ncol(mat1)==ncol(mat2))
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  if(length(alt)==1) alt <- rep(alt, nrow(mat1))
  if(length(conf)==1) conf <- rep(conf, nrow(mat1))

  np <- cor <- tst <- p <- cl <- ch <- df <- mu <- cnf <- numeric(nrow(mat1))
  al <- character(nrow(mat1))
  for(i in 1:nrow(mat1)) {
    good <- complete.cases(mat1[i,], mat2[i,])
    vec1 <- mat1[i,good]
    vec2 <- mat2[i,good]
    res <- cor.test(vec1, vec2, alternative=alt[i], conf.level=conf[i],
                    method="pearson"
                    )

    np[i]  <- length(vec1)
    cor[i] <- res$estimate
    tst[i] <- res$statistic
    p[i]   <- res$p.value
    cl[i]  <- ifelse(is.null(res$conf.int), NA, res$conf.int[1])
    ch[i]  <- ifelse(is.null(res$conf.int), NA, res$conf.int[2])
    cnf[i] <- ifelse(is.null(res$conf.int), conf, attr(res$conf.int, "conf.level"))
    df[i]  <- res$parameter
    mu[i]  <- res$null.value
    al[i]  <- res$alternative
  }

  data.frame(obs.paired=np, cor=cor, df=df, statistic.t=tst,
             p.value=p, ci.low=cl, ci.high=ch, alternative=al,
             cor.null=mu, conf.level=cnf, stringsAsFactors=FALSE
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
  Y[sample(length(Y), nrow(Y))] <- NA
  alts <- rep(c("t", "g", "l"), length.out=nrow(X))
  cfs  <- seq(0, 1, length.out=nrow(X))

  t1 <- base_cor_pearson(X, Y, alts, cfs)
  t2 <- row.cor.pearson(X, Y, alts, cfs)

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  vals1 <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
  vals2 <- c(100000000000003, 100000000000002, 100000000000003, 100000000000000)
  t1 <- base_cor_pearson(vals1, vals2)
  t2 <- row.cor.pearson(vals1, vals2)
  expect_equal(t1, t2)

  # small numbers
  vals1 <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1)
  vals2 <- c(1.00000000000003, 1.00000000000002, 1.00000000000003, 1)
  t1 <- base_cor_pearson(vals1, vals2)
  t2 <- row.cor.pearson(vals1, vals2)
  expect_equal(t1, t2)
})

test_that("perfect correlations give equal results", {
  # positive case
  t1 <- base_cor_pearson(1:4, 1:4)
  t2 <- suppressWarnings(row.cor.pearson(1:4, 1:4))
  expect_equal(t1, t2)

  # negative case
  t1 <- base_cor_pearson(1:4, 4:1)
  t2 <- suppressWarnings(row.cor.pearson(1:4, 4:1))
  expect_equal(t1, t2)
})

test_that("minumum allowed sample sizes give equal results", {
  # three numbers
  x <- matrix(rnorm(9), ncol=3); y <- matrix(rnorm(9), ncol=3)
  alt <- c("two.sided", "greater", "less")
  t1 <- base_cor_pearson(x, y)
  t2 <- suppressWarnings(row.cor.pearson(x, y))
  expect_equal(t1, t2)

  # three numbers with NAs
  x <- matrix(rnorm(15), ncol=5); x[,4] <- NA
  y <- matrix(rnorm(15), ncol=5); y[,5] <- NA
  t1 <- base_cor_pearson(x, y)
  t2 <- suppressWarnings(row.cor.pearson(x, y))
  expect_equal(t1, t2)

  # four numbers (will produce confidence intervals)
  x <- matrix(rnorm(12), ncol=4); y <- matrix(rnorm(12), ncol=4)
  alt <- c("two.sided", "greater", "less")
  t1 <- base_cor_pearson(x, y)
  t2 <- row.cor.pearson(x, y)
  expect_equal(t1, t2)
})

test_that("parameter edge cases give equal results", {
  set.seed(14)
  alt <- c("l", "t", "g")
  cfs <- c(0, 0.5, 1)
  pars <- expand.grid(alt, cfs, stringsAsFactors=FALSE)
  X1 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X2 <- matrix(rnorm(10*nrow(pars)), ncol=10)
  X1[sample(length(X1), nrow(pars)*2)] <- NA
  X2[sample(length(X2), nrow(pars)*2)] <- NA

  t1 <- base_cor_pearson(X1, X2, pars[,1], pars[,2])
  t2 <- suppressWarnings(row.cor.pearson(X1, X2, pars[,1], pars[,2]))

  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warning when rows have exactly 3 complete observations", {
  wrn <- '1 of the rows had exactly 3 complete observations: no confidence intervals produced\\. First occurrence at row 1'
  nacolumns <- c("ci.low", "ci.high")

  # standard case
  expect_warning(res <- row.cor.pearson(c(1,1,2), c(1,2,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_true(all(!is.na(res[,!colnames(res) %in% nacolumns])))

  # with NAs present
  expect_warning(res <- row.cor.pearson(c(1,1,4,NA,4), c(1,2,1,3,NA)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_true(all(!is.na(res[,!colnames(res) %in% nacolumns])))
})


test_that("warning when rows have less than 3 complete observations", {
  wrn <- '1 of the rows had less than 3 complete observations\\. First occurrence at row 1'
  nacolumns <- c("statistic.t", "p.value", "ci.low", "ci.high")

  # 1 observations
  expect_warning(res <- row.cor.pearson(1, 1), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # 2 observations
  expect_warning(res <- row.cor.pearson(c(1,2), c(1,3)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # with NAs present
  expect_warning(res <- row.cor.pearson(c(2,1,NA,4), c(1,2,1,NA)), wrn, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
})

test_that("warning when one of the variables has zero standard deviation", {
  wrnX <- '1 of the rows had zero standard deviation in x\\. First occurrence at row 1'
  wrnY <- '1 of the rows had zero standard deviation in y\\. First occurrence at row 1'
  nacolumns <- c("statistic.t", "p.value", "ci.low", "ci.high")

  # only x
  expect_warning(res <- row.cor.pearson(c(1,1,1,1), c(1,2,3,4)), wrnX, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # only y
  expect_warning(res <- row.cor.pearson(c(1,1,2,2), c(2,2,2,2)), wrnY, all=TRUE)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # both
  expect_warning(res <- row.cor.pearson(c(1,1,1), c(2,2,2)), wrnX)
  expect_warning(res <- row.cor.pearson(c(1,1,1), c(2,2,2)), wrnY)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))

  # with NAs present
  expect_warning(res <- row.cor.pearson(c(1,1,1,NA,4), c(2,2,2,4,NA)), wrnX)
  expect_warning(res <- row.cor.pearson(c(1,1,1,NA,4), c(2,2,2,4,NA)), wrnY)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
})

test_that("warning when correlation is perfect", {
  wrn <- '1 of the rows had essentially perfect fit: results might be unreliable for small sample sizes\\. First occurrence at row 1'

  # positive case
  expect_warning(res <- row.cor.pearson(c(1:4), c(1:4)), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic.t))

  # negative case
  expect_warning(res <- row.cor.pearson(c(1:4), c(4:1)), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic.t))

  # with NAs present
  expect_warning(res <- row.cor.pearson(c(1,2,3,4,NA,1), c(1,2,3,4,1,NA)), wrn, all=TRUE)
  expect_true(is.infinite(res$statistic.t))
})


