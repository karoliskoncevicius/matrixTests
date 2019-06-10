context("correctness of wilcoxon_twosample")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_wilcoxon_twosample <- function(mat1, mat2, alt="two.sided", mu=0, exact=NA, correct=TRUE) {
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  if(length(alt)==1) alt <- rep(alt, nrow(mat1))
  if(length(mu)==1) mu <- rep(mu, nrow(mat1))
  if(length(exact)==1) exact <- rep(exact, nrow(mat1))
  if(length(correct)==1) correct <- rep(correct, nrow(mat1))

  nx <- ny <- nxy <- stat <- p <- m0 <- numeric(nrow(mat1))
  al <- character(nrow(mat1))
  ext <- cor <- logical(nrow(mat1))
  for(i in 1:nrow(mat1)) {
    ex  <- if(is.na(exact[i])) NULL else exact[i]
    res <- wilcox.test(mat1[i,], mat2[i,], alternative=alt[i], mu=mu[i], exact=ex, correct=correct[i])

    mat1[i,!is.finite(mat1[i,])] <- NA
    mat2[i,!is.finite(mat2[i,])] <- NA
    nx[i]   <- sum(!is.na(mat1[i,]))
    ny[i]   <- sum(!is.na(mat2[i,]))
    nxy[i]  <- nx[i] + ny[i]
    stat[i] <- res$statistic
    p[i]    <- res$p.value
    al[i]   <- res$alternative
    m0[i]   <- res$null.value

    cor[i] <- grepl("continuity correction", res$method)
    ext[i] <- if(is.na(exact[i])) (nx[i] < 50 & ny[i] < 50) else exact[i]
    ext[i] <- if(length(unique(na.omit(c(mat1[i,]-m0[i], mat2[i,])))) != nxy[i]) FALSE else ext[i]
  }

  data.frame(obs.x=nx, obs.y=ny, obs.tot=nxy, statistic=stat, pvalue=p,
             alternative=al, location.null=m0, exact=ext, corrected=cor,
             stringsAsFactors=FALSE
             )
}

################################################################################
########################### TEST ON A RANDOM SAMPLE ############################
################################################################################

############################## SMALL SAMPLE SIZE ###############################

test_that("monte-carlo random testing gives equal results 1", {
  set.seed(14)
  X <- matrix(rnorm(10000), ncol=10)
  Y <- matrix(rnorm(20000), ncol=20)
  X[sample(length(X), nrow(X))] <- NA
  Y[sample(length(Y), 2*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

############################### BIG SAMPLE SIZE ################################

test_that("monte-carlo random testing gives equal results 2", {
  set.seed(14)
  X <- matrix(rnorm(70000), ncol=70)
  Y <- matrix(rnorm(60000), ncol=60)
  X[sample(length(X), 7*nrow(X))] <- NA
  Y[sample(length(Y), 8*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

######################### BORDERLINE SAMPLE SIZE OF 50 #########################

test_that("monte-carlo random testing gives equal results 3", {
  set.seed(14)
  X <- matrix(rnorm(52000), ncol=52)
  Y <- matrix(rnorm(52000), ncol=52)
  X[sample(length(X), 5*nrow(X))] <- NA
  Y[sample(length(Y), 5*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

################### ONE GROUPS IS ABOVE 50 AND ANOTER BELOW ####################

test_that("monte-carlo random testing gives equal results 4", {
  set.seed(14)
  X <- matrix(rnorm(70000), ncol=70)
  Y <- matrix(rnorm(15000), ncol=15)
  X[sample(length(X), 7*nrow(X))] <- NA
  Y[sample(length(Y), 1*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

######################### SMALL SAMPLE SIZE WITH TIES ##########################

test_that("monte-carlo random testing gives equal results 5", {
  set.seed(14)
  X <- matrix(round(runif(10000, -15, 15)), ncol=10)
  Y <- matrix(round(runif(11000, -15, 15)), ncol=11)
  X[sample(length(X), nrow(X))] <- NA
  Y[sample(length(Y), nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

########################## BIG SAMPLE SIZE WITH TIES ###########################

test_that("monte-carlo random testing gives equal results 6", {
  set.seed(14)
  X <- matrix(round(runif(60000, -60, 60)), ncol=60)
  Y <- matrix(round(runif(60000, -60, 60)), ncol=60)
  X[sample(length(X), 10*nrow(X))] <- NA
  Y[sample(length(Y), 10*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

####################### BORDERLINE SAMPLE SIZE WITH TIES #######################

test_that("monte-carlo random testing gives equal results 7", {
  set.seed(14)
  X <- matrix(round(runif(52000, -50, 50)), ncol=52)
  Y <- matrix(round(runif(53000, -50, 50)), ncol=53)
  X[sample(length(X), 5*nrow(X))] <- NA
  Y[sample(length(Y), 5*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

################### ONE BIG GROUP ONE SMALL GROUP WITH TIES ####################

test_that("monte-carlo random testing gives equal results 8", {
  set.seed(14)
  X <- matrix(round(runif(15000, -20, 20)), ncol=15)
  Y <- matrix(round(runif(70000, -15, 15)), ncol=70)
  X[sample(length(X), 2*nrow(X))] <- NA
  Y[sample(length(Y), 7*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_twosample(X, Y, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  vals1 <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
  vals2 <- c(100000000000005, 100000000000006, 100000000000009)
  mat1  <- matrix(vals1, nrow=nrow(pars), ncol=length(vals1), byrow=TRUE)
  mat2  <- matrix(vals2, nrow=nrow(pars), ncol=length(vals2), byrow=TRUE)
  t1 <- base_wilcoxon_twosample(mat1, mat2, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_twosample(mat1, mat2, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)

  # small numbers
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  vals1 <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0.00000000000001)
  vals2 <- c(0.00000000000005, 0.00000000000006, 0.00000000000000, 0.00000000000007)
  mat1  <- matrix(vals1, nrow=nrow(pars), ncol=length(vals1), byrow=TRUE)
  mat2  <- matrix(vals2, nrow=nrow(pars), ncol=length(vals2), byrow=TRUE)
  t1 <- base_wilcoxon_twosample(mat1, mat2, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_twosample(mat1, mat2, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)
})


test_that("minumum allowed sample sizes give equal results", {
  # single number in each group
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x1  <- matrix(rnorm(nrow(pars)), ncol=1)
  x2  <- matrix(rnorm(nrow(pars)), ncol=1)
  t1 <- base_wilcoxon_twosample(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_twosample(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)
  expect_true(all(t2$obs.x==1))
  expect_true(all(t2$obs.y==1))
  expect_true(all(t2$obs.tot==2))

  # 3 observation in both, but only one non NA or Inf
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x1  <- cbind(matrix(rnorm(nrow(pars)), ncol=1), NA, Inf)
  x2  <- cbind(-Inf, NA, matrix(rnorm(nrow(pars)), ncol=1))
  t1 <- base_wilcoxon_twosample(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- suppressWarnings(row_wilcoxon_twosample(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3]))
  expect_equal(t1, t2)
  expect_true(all(t2$obs.x==1))
  expect_true(all(t2$obs.y==1))
  expect_true(all(t2$obs.tot==2))

  # three numbers in both grouops, all equal
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x1  <- matrix(rnorm(nrow(pars)), ncol=1)
  x1  <- cbind(x1,x1,x1)
  x2  <- x1
  t1 <- suppressWarnings(base_wilcoxon_twosample(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3]))
  t2 <- suppressWarnings(row_wilcoxon_twosample(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3]))
  expect_equal(t1, t2)

  # three numbers all equal after subtracting mu from x
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x1 <- matrix(rnorm(nrow(pars)), ncol=1)
  x1 <- cbind(x1, x1, x1)
  x2 <- x1 - 1
  t1 <- suppressWarnings(base_wilcoxon_twosample(x1, x2, pars[,1], mu=1, exact=pars[,2], correct=pars[,3]))
  t2 <- suppressWarnings(row_wilcoxon_twosample(x1, x2, pars[,1], mu=1, exact=pars[,2], correct=pars[,3]))
  expect_equal(t1, t2)
})


test_that("parameter edge cases give equal results", {
  set.seed(14)
  alt <- c("l", "t", "g")
  mus <- c(-10000, 0, 10000)
  ext <- c(TRUE, FALSE, NA)
  cor <- c(TRUE, FALSE)
  pars <- expand.grid(alt, mus, ext, cor, stringsAsFactors=FALSE)
  X <- matrix(round(runif(10*nrow(pars), -25, 25)), ncol=10)
  X[sample(length(X), nrow(pars)*2)] <- NA
  Y <- matrix(round(runif(8*nrow(pars), -25, 25)), ncol=8)
  Y[sample(length(Y), nrow(pars)*2)] <- NA

  t1 <- suppressWarnings(base_wilcoxon_twosample(X, Y, pars[,1], pars[,2], pars[,3], pars[,4]))
  t2 <- suppressWarnings(row_wilcoxon_twosample(X, Y, pars[,1], pars[,2], pars[,3], pars[,4]))
  expect_equal(t1, t2)
})

################################################################################
#################### INTERACTION BETWEEN EXACT AND CORRECT #####################
################################################################################

test_that("small sample size automatically turns exact=TRUE", {
  # 49 samples in both groups
  x  <- rnorm(49)
  y  <- rnorm(49)
  t1 <- base_wilcoxon_twosample(x, y)
  t2 <- row_wilcoxon_twosample(x, y, exact=NA)
  t3 <- row_wilcoxon_twosample(x, y, exact=TRUE)
  expect_true(t2$exact)
  expect_equal(t2$obs.x, 49)
  expect_equal(t2$obs.y, 49)
  expect_equal(t2$obs.tot, 49*2)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 50 samples in both groups - one with NA
  x  <- c(rnorm(49), NA)
  y  <- c(NA, rnorm(49))
  t1 <- base_wilcoxon_twosample(x, y)
  t2 <- row_wilcoxon_twosample(x, y, exact=NA)
  t3 <- row_wilcoxon_twosample(x, y, exact=TRUE)
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t2$obs.x, 49)
  expect_equal(t2$obs.y, 49)
  expect_equal(t2$obs.tot, 2*49)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
})

test_that("big sample size automatically turns exact=FALSE", {
  # 50 samples
  x  <- rnorm(50)
  y  <- rnorm(50)
  t1 <- base_wilcoxon_twosample(x, y)
  t2 <- row_wilcoxon_twosample(x, y, exact=NA)
  t3 <- row_wilcoxon_twosample(x, y, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.x, 50)
  expect_equal(t2$obs.y, 50)
  expect_equal(t2$obs.tot, 100)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 51 samples with one NA
  x  <- c(NA, rnorm(50))
  y  <- c(rnorm(50), NA)
  t1 <- base_wilcoxon_twosample(x, y)
  t2 <- row_wilcoxon_twosample(x, y, exact=NA)
  t3 <- row_wilcoxon_twosample(x, y, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.x, 50)
  expect_equal(t2$obs.y, 50)
  expect_equal(t2$obs.tot, 100)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # One small group and another with 50 values
  x  <- rnorm(10)
  y  <- rnorm(50)
  t1 <- suppressWarnings(base_wilcoxon_twosample(x, y))
  t2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.x, 10)
  expect_equal(t2$obs.y, 50)
  expect_equal(t2$obs.tot, 60)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
})

test_that("ties automatically turn exact=FALSE, even when exact=TRUE specified", {
  # 10 samples with one tie after subtracting mu=1
  x  <- c(rnorm(9), 2)
  y  <- c(rnorm(9), 1)
  t1 <- suppressWarnings(base_wilcoxon_twosample(x, y, mu=1))
  t2 <- suppressWarnings(row_wilcoxon_twosample(x, y, mu=1, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_twosample(x, y, mu=1, exact=TRUE))
  t4 <- row_wilcoxon_twosample(x, y, mu=1, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
  # 105 samples with 5 ties
  x  <- c(rnorm(50), 2, 2, 2, 2, 2)
  y  <- c(rnorm(50))
  t1 <- suppressWarnings(base_wilcoxon_twosample(x, y))
  t2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=TRUE))
  t4 <- row_wilcoxon_twosample(x, y, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
})

test_that("correct=TRUE/FALSE only applies when exact=FALSE", {
  # big sample size so exact=FALSE
  x  <- rnorm(50)
  y  <- rnorm(50)
  t1 <- base_wilcoxon_twosample(x, y, correct=FALSE)
  t2 <- row_wilcoxon_twosample(x, y, exact=NA, correct=FALSE)
  t3 <- row_wilcoxon_twosample(x, y, exact=FALSE, correct=FALSE)
  expect_false(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # big sample size, but exact=FALSE becasue of ties
  x  <- c(rnorm(50), 2, 2)
  y  <- c(rnorm(40), 4, 4)
  t1 <- suppressWarnings(base_wilcoxon_twosample(x, y, exact=TRUE, correct=FALSE))
  t2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=NA, correct=FALSE))
  t3 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=TRUE, correct=FALSE))
  expect_false(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # small sample size but exact is turned off
  x  <- rnorm(20)
  y  <- rnorm(20)
  t1 <- base_wilcoxon_twosample(x, y, exact=FALSE, correct=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=FALSE, correct=FALSE))
  expect_false(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
})

test_that("correct=TRUE is turned to FALSE when exact=TRUE", {
  # big sample size with exact=TRUE
  x  <- rnorm(60)
  y  <- rnorm(60)
  t1 <- base_wilcoxon_twosample(x, y, exact=TRUE, correct=TRUE)
  t2 <- row_wilcoxon_twosample(x, y, exact=TRUE, correct=TRUE)
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  # big sample size but NAs make it less than 50 so exact is turned on
  x  <- c(rnorm(49), NA, NaN)
  y  <- c(rnorm(49), NA, NaN)
  t1 <- base_wilcoxon_twosample(x, y, correct=TRUE)
  t2 <- row_wilcoxon_twosample(x, y, correct=TRUE)
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  # small sample size and exact is left as NA
  x  <- rnorm(20)
  y  <- rnorm(15)
  t1 <- base_wilcoxon_twosample(x, y, exact=NA, correct=TRUE)
  t2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=NA, correct=TRUE))
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warnign when row has less than 1 available observations", {
  wrnX <- 'row_wilcoxon_twosample: 1 of the rows had less than 1 remaining finite "x" observation'
  wrnY <- 'row_wilcoxon_twosample: 1 of the rows had less than 1 remaining finite "y" observation'
  nacolumns <- c("statistic", "pvalue")

  # no x observations
  expect_warning(res <- row_wilcoxon_twosample(numeric(), 1), wrnX)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 1)
  expect_equal(res$obs.tot, 1)

  # no y observations
  expect_warning(res <- row_wilcoxon_twosample(1, numeric()), wrnY)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 1)
  expect_equal(res$obs.y, 0)
  expect_equal(res$obs.tot, 1)

  # no x no y observations
  expect_warning(res <- row_wilcoxon_twosample(numeric(), numeric()), wrnX)
  expect_warning(res <- row_wilcoxon_twosample(numeric(), numeric()), wrnY)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 0)
  expect_equal(res$obs.tot, 0)

  # only NAs in x
  expect_warning(res <- row_wilcoxon_twosample(c(NA, NaN, NA), 1:10), wrnX)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 10)
  expect_equal(res$obs.tot, 10)

  # only NAs in y
  expect_warning(res <- row_wilcoxon_twosample(1:10, c(NA, NaN, NA)), wrnY)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 10)
  expect_equal(res$obs.y, 0)
  expect_equal(res$obs.tot, 10)

  # only NAs in both
  expect_warning(res <- row_wilcoxon_twosample(numeric(), c(NA, NaN, NA)), wrnX)
  expect_warning(res <- row_wilcoxon_twosample(numeric(), c(NA, NaN, NA)), wrnY)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 0)
  expect_equal(res$obs.tot, 0)
})

test_that("warning when a infinite values are removed", {
  wrnX <- 'row_wilcoxon_twosample: 1 of the rows had infinite "x" observations that were removed.'
  wrnY <- 'row_wilcoxon_twosample: 1 of the rows had infinite "y" observations that were removed.'

  # -Inf and Inf among x observations
  expect_warning(res <- row_wilcoxon_twosample(c(-Inf,5,7,Inf), 1:4), wrnX, all=TRUE)
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 4)
  expect_equal(res$obs.tot, 6)

  # -Inf and Inf among y observations
  expect_warning(res <- row_wilcoxon_twosample(1:4, c(-Inf,5,7,Inf)), wrnY, all=TRUE)
  expect_equal(res$obs.x, 4)
  expect_equal(res$obs.y, 2)
  expect_equal(res$obs.tot, 6)

  # -Inf and Inf among both
  expect_warning(res <- row_wilcoxon_twosample(c(1,2,3,Inf), c(-Inf,7,6,5)), wrnX)
  expect_warning(res <- row_wilcoxon_twosample(c(1,2,3,Inf), c(-Inf,7,6,5)), wrnY)
  expect_equal(res$obs.x, 3)
  expect_equal(res$obs.y, 3)
  expect_equal(res$obs.tot, 6)
})

test_that("warning when ties are present", {
  wrn <- 'row_wilcoxon_twosample: 1 of the rows had ties: cannot compute exact p-values with ties.'

  # warning when exact=TRUE
  expect_warning(res <- row_wilcoxon_twosample(c(3,3,1,2), c(2,0,5,-2), exact=TRUE), wrn, all=TRUE)
  expect_equal(res$obs.x, 4)
  expect_equal(res$obs.y, 4)
  expect_equal(res$obs.tot, 8)
  expect_false(res$exact)

  # no warning when exact=FALSE
  expect_warning(res <- row_wilcoxon_twosample(c(3,3,1,2), c(2,2,5,-2), exact=FALSE), NA)
  expect_false(res$exact)

  # ties only after subtracting mu
  expect_warning(res <- row_wilcoxon_twosample(c(3.1,4,1,2), c(-1,-2,3,-4), mu=0.1, exact=TRUE), wrn, all=TRUE)
  expect_equal(res$obs.x, 4)
  expect_equal(res$obs.y, 4)
  expect_equal(res$obs.tot, 8)
  expect_false(res$exact)
})

