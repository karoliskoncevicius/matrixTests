context("correctness of wilcoxon_paired")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_wilcoxon_paired <- function(mat1, mat2, alt="two.sided", mu=0, exact=NA, correct=TRUE) {
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
    res <- wilcox.test(mat1[i,], mat2[i,], alternative=alt[i], mu=mu[i], paired=TRUE, exact=ex, correct=correct[i])

    mat1[i,!is.finite(mat1[i,])] <- NA
    mat2[i,!is.finite(mat2[i,])] <- NA
    inds <- complete.cases(mat1[i,], mat2[i,])
    vec1 <- mat1[i,inds]
    vec2 <- mat2[i,inds]
    nx[i]   <- sum(!is.na(mat1[i,]))
    ny[i]   <- sum(!is.na(mat2[i,]))
    nxy[i]  <- sum((vec1 - vec2 - mu[i]) != 0, na.rm=TRUE)
    stat[i] <- res$statistic
    p[i]    <- res$p.value
    al[i]   <- res$alternative
    m0[i]   <- res$null.value

    cor[i] <- grepl("continuity correction", res$method)
    ext[i] <- if(is.na(exact[i])) (nxy[i] < 50) else exact[i]
    ext[i] <- if(length(unique(abs(vec1-vec2-m0[i]))) != length(vec1)) FALSE else ext[i]
    ext[i] <- if(any((vec1-vec2-m0[i]) == 0)) FALSE else ext[i]
  }

  data.frame(obs.x=nx, obs.y=ny, obs.paired=nxy, statistic=stat, pvalue=p,
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
  X <- matrix(rnorm(100000), ncol=10)
  Y <- matrix(rnorm(100000), ncol=10)
  X[sample(length(X), nrow(X))] <- NA
  Y[sample(length(Y), nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

############################### BIG SAMPLE SIZE ################################

test_that("monte-carlo random testing gives equal results 2", {
  set.seed(14)
  X <- matrix(rnorm(1000000), ncol=100)
  Y <- matrix(rnorm(1000000), ncol=100)
  X[sample(length(X), 10*nrow(X))] <- NA
  Y[sample(length(Y), 10*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

######################### BORDERLINE SAMPLE SIZE OF 50 #########################

test_that("monte-carlo random testing gives equal results 3", {
  set.seed(14)
  X <- matrix(rnorm(500000), ncol=50)
  Y <- matrix(rnorm(500000), ncol=50)
  X[sample(length(X), 5*nrow(X))] <- NA
  Y[sample(length(Y), 5*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

#################### SMALL SAMPLE SIZE WITH TIES AND ZEROES ####################

test_that("monte-carlo random testing gives equal results 4", {
  set.seed(14)
  X <- matrix(round(runif(100000, -15, 15)), ncol=10)
  Y <- matrix(round(runif(100000, -15, 15)), ncol=10)
  X[sample(length(X), nrow(X))] <- NA
  Y[sample(length(Y), nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

##################### BIG SAMPLE SIZE WITH TIES AND ZEROES #####################

test_that("monte-carlo random testing gives equal results 5", {
  set.seed(14)
  X <- matrix(round(runif(1000000, -150, 15)), ncol=100)
  Y <- matrix(round(runif(1000000, -150, 15)), ncol=100)
  X[sample(length(X), 10*nrow(X))] <- NA
  Y[sample(length(Y), 10*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

################# BORDERLINE SAMPLE SIZE WITH TIES AND ZEROES ##################

test_that("monte-carlo random testing gives equal results 6", {
  set.seed(14)
  X <- matrix(round(runif(500000, -150, 15)), ncol=50)
  Y <- matrix(round(runif(500000, -150, 15)), ncol=50)
  X[sample(length(X), 5*nrow(X))] <- NA
  Y[sample(length(Y), 5*nrow(Y))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_paired(X, Y, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  vals1 <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
  vals2 <- c(100000000000002, 100000000000003, 100000000000000, 100000000000004)
  mat1  <- matrix(vals1, nrow=nrow(pars), ncol=length(vals1), byrow=TRUE)
  mat2  <- matrix(vals2, nrow=nrow(pars), ncol=length(vals2), byrow=TRUE)
  t1 <- base_wilcoxon_paired(mat1, mat2, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_paired(mat1, mat2, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)

  # small numbers
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  vals1 <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0.00000000000001)
  vals2 <- c(0.00000000000002, 0.00000000000003, 0.00000000000000, 0.00000000000005)
  mat1  <- matrix(vals1, nrow=nrow(pars), ncol=length(vals1), byrow=TRUE)
  mat2  <- matrix(vals2, nrow=nrow(pars), ncol=length(vals2), byrow=TRUE)
  t1 <- base_wilcoxon_paired(mat1, mat2, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_paired(mat1, mat2, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)
})


test_that("minumum allowed sample sizes give equal results", {
  # single number
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x1  <- matrix(rnorm(nrow(pars)), ncol=1)
  x2  <- matrix(rnorm(nrow(pars)), ncol=1)
  t1 <- base_wilcoxon_paired(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_paired(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)
  expect_true(all(t2$obs.paired==1))

  # 3 observation in both, but only one common
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x1  <- matrix(rnorm(nrow(pars)*3), ncol=3)
  x2  <- matrix(rnorm(nrow(pars)*3), ncol=3)
  x1[,2] <- NA
  x2[,3] <- NA
  t1 <- base_wilcoxon_paired(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_paired(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)
  expect_true(all(t2$obs.paired==1))

  # three numbers all equal
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x1  <- matrix(rnorm(nrow(pars)), ncol=1)
  x1  <- cbind(x1,x1,x1)
  x2  <- matrix(rnorm(nrow(pars)), ncol=1)
  x2  <- cbind(x2,x2,x2)
  t1 <- suppressWarnings(base_wilcoxon_paired(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3]))
  t2 <- suppressWarnings(row_wilcoxon_paired(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3]))
  expect_equal(t1, t2)

  # three numbers all differences equal
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x1 <- matrix(rnorm(nrow(pars)*3), ncol=3)
  x2 <- x1 + 0.2
  t1 <- suppressWarnings(base_wilcoxon_paired(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3]))
  t2 <- suppressWarnings(row_wilcoxon_paired(x1, x2, pars[,1], exact=pars[,2], correct=pars[,3]))
  expect_equal(t1, t2)

  # three numbers, two equal to mu after subtraction
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x1 <- matrix(rnorm(nrow(pars)*3), ncol=3)
  x2 <- x1 -1
  x2[,3] <- x2[,3] + rnorm(nrow(x2))
  t1 <- suppressWarnings(base_wilcoxon_paired(x1, x2, pars[,1], 1, exact=pars[,2], correct=pars[,3]))
  t2 <- suppressWarnings(row_wilcoxon_paired(x1, x2, pars[,1], 1, exact=pars[,2], correct=pars[,3]))
  expect_equal(t1, t2)
})


test_that("parameter edge cases give equal results", {
  set.seed(14)
  alt <- c("l", "t", "g")
  mus <- c(-10000, 0, 10000)
  ext <- c(TRUE, FALSE, NA)
  cor <- c(TRUE, FALSE)
  pars <- expand.grid(alt, mus, ext, cor, stringsAsFactors=FALSE)
  X <- matrix(round(runif(10*nrow(pars), -15, 15)), ncol=10)
  X[sample(length(X), nrow(pars)*2)] <- NA
  Y <- matrix(round(runif(10*nrow(pars), -15, 15)), ncol=10)
  Y[sample(length(Y), nrow(pars)*2)] <- NA

  t1 <- suppressWarnings(base_wilcoxon_paired(X, Y, pars[,1], pars[,2], pars[,3], pars[,4]))
  t2 <- suppressWarnings(row_wilcoxon_paired(X, Y, pars[,1], pars[,2], pars[,3], pars[,4]))
  expect_equal(t1, t2)
})

################################################################################
#################### INTERACTION BETWEEN EXACT AND CORRECT #####################
################################################################################

test_that("small sample size automatically turns exact=TRUE", {
  # 49 samples
  x  <- rnorm(49)
  y  <- rnorm(49)
  t1 <- base_wilcoxon_paired(x, y)
  t2 <- row_wilcoxon_paired(x, y, exact=NA)
  t3 <- row_wilcoxon_paired(x, y, exact=TRUE)
  expect_true(t2$exact)
  expect_equal(t2$obs.paired, 49)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 50 samples with one NA
  x  <- rnorm(49); x  <- c(x, NA)
  y  <- rnorm(50)
  t1 <- base_wilcoxon_paired(x, y)
  t2 <- row_wilcoxon_paired(x, y, exact=NA)
  t3 <- row_wilcoxon_paired(x, y, exact=TRUE)
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t2$obs.paired, 49)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 51 samples with Na and NaN
  x  <- rnorm(50); x  <- c(x, NA)
  y  <- rnorm(50); y  <- c(NaN, y)
  t1 <- base_wilcoxon_paired(x, y)
  t2 <- row_wilcoxon_paired(x, y, exact=NA)
  t3 <- row_wilcoxon_paired(x, y, exact=TRUE)
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t2$obs.paired, 49)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
})

test_that("big sample size automatically turns exact=FALSE", {
  # 50 samples
  x  <- rnorm(50)
  y  <- rnorm(50)
  t1 <- base_wilcoxon_paired(x, y)
  t2 <- row_wilcoxon_paired(x, y, exact=NA)
  t3 <- row_wilcoxon_paired(x, y, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.paired, 50)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 51 samples with one NA
  x  <- rnorm(51)
  y  <- c(rnorm(50), NA)
  t1 <- base_wilcoxon_paired(x, y)
  t2 <- row_wilcoxon_paired(x, y, exact=NA)
  t3 <- row_wilcoxon_paired(x, y, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.paired, 50)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 51 samples with 1 value equal to mu
  x  <- rnorm(51)
  y  <- rnorm(51)
  y[51] <- x[51]
  t1 <- suppressWarnings(base_wilcoxon_paired(x, y))
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.paired, 50)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
})

test_that("ties automatically turn exact=FALSE, even when exact=TRUE specified", {
  # 11 samples with one tie
  x  <- c(rnorm(9), 2, 2)
  y  <- c(rnorm(9), 1, 1)
  t1 <- suppressWarnings(base_wilcoxon_paired(x, y))
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=TRUE))
  t4 <- row_wilcoxon_paired(x, y, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
  # 105 samples with 5 ties
  x  <- c(rnorm(100), 2, 2, 2, 2, 2)
  y  <- c(rnorm(100), 1, 1, 1, 1, 1)
  t1 <- suppressWarnings(base_wilcoxon_paired(x, y))
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=TRUE))
  t4 <- row_wilcoxon_paired(x, y, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
})

test_that("zeroes automatically turn exact=FALSE, even when exact=TRUE specified", {
  # 11 samples with one zero
  x  <- c(rnorm(10), 2)
  y  <- c(rnorm(10), 1)
  t1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE))
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE))
  t4 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.paired, 10)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
  # 105 samples with 5 zeroes
  x  <- c(rnorm(100), 2, 2, 2, 2, 2)
  y  <- c(rnorm(100), 1, 1, 1, 1, 1)
  t1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE))
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE))
  t4 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.paired, 100)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
})

test_that("mix of zeroes and ties automatically turn exact=FALSE, even when exact=TRUE specified", {
  # 11 samples with one zero and one tie
  x  <- c(rnorm(8), 2, 2, 3)
  y  <- c(rnorm(8), 4, 4, 2)
  t1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE))
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE))
  t4 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.paired, 10)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
  # 110 samples with 5 zeroes and 5 ties
  x  <- c(rnorm(100), 2,2,2,2,2, 3,3,3,3,3)
  y  <- c(rnorm(100), 4,4,4,4,4, 2,2,2,2,2)
  t1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE))
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE))
  t4 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t2$obs.paired, 105)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
})

test_that("correct=TRUE/FALSE only applies when exact=FALSE", {
  # big sample size so exact=FALSE
  x  <- rnorm(100)
  y  <- rnorm(100)
  t1 <- base_wilcoxon_paired(x, y, correct=FALSE)
  t2 <- row_wilcoxon_paired(x, y, exact=NA, correct=FALSE)
  t3 <- row_wilcoxon_paired(x, y, exact=FALSE, correct=FALSE)
  expect_false(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # big sample size, but exact=FALSE becasue of ties and zeroes
  x  <- c(rnorm(100), 2, 2, 3)
  y  <- c(rnorm(100), 4, 4, 2)
  t1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE, correct=FALSE))
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA, correct=FALSE))
  t3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE, correct=FALSE))
  expect_false(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # small sample size but exact is turned off
  x  <- rnorm(20)
  y  <- rnorm(20)
  t1 <- base_wilcoxon_paired(x, y, exact=FALSE, correct=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=FALSE, correct=FALSE))
  expect_false(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
})

test_that("correct=TRUE is turned to FALSE when exact=TRUE", {
  # big sample size with exact=TRUE
  x  <- rnorm(100)
  y  <- rnorm(100)
  t1 <- base_wilcoxon_paired(x, y, exact=TRUE, correct=TRUE)
  t2 <- row_wilcoxon_paired(x, y, exact=TRUE, correct=TRUE)
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  # big sample size but NAs make it less than 50 so exact is turned on
  x  <- rnorm(51)
  y  <- c(rnorm(49), NaN, NA)
  t1 <- base_wilcoxon_paired(x, y, correct=TRUE)
  t2 <- row_wilcoxon_paired(x, y, correct=TRUE)
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  # small sample size and exact is left as NA
  x  <- rnorm(20)
  y  <- rnorm(20)
  t1 <- base_wilcoxon_paired(x, y, exact=NA, correct=TRUE)
  t2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=NA, correct=TRUE))
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warnign when row has less than 1 available observations", {
  wrn <- 'row_wilcoxon_paired: 1 of the rows had less than 1 remaining finite paired "x-y" observation.'
  nacolumns <- c("statistic", "pvalue")

  # no observations
  expect_warning(res <- row_wilcoxon_paired(numeric(), numeric()), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 0)
  expect_equal(res$obs.paired, 0)

  # only NAs in one side
  expect_warning(res <- row_wilcoxon_paired(1:3, c(NA, NaN, NA)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 3)
  expect_equal(res$obs.y, 0)
  expect_equal(res$obs.paired, 0)

  # no common complete observations
  expect_warning(res <- row_wilcoxon_paired(c(1, NA, 3, NA), c(NA, 2, NA, 4)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 2)
  expect_equal(res$obs.paired, 0)

  # only NA observations
  expect_warning(res <- row_wilcoxon_paired(c(NA,NaN,NA), c(NA, NA, NaN)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs.x, 0)
  expect_equal(res$obs.y, 0)
  expect_equal(res$obs.paired, 0)
})

test_that("warning when a infinite values are removed", {
  wrnX <- 'row_wilcoxon_paired: 1 of the rows had infinite "x" observations that were removed.'
  wrnY <- 'row_wilcoxon_paired: 1 of the rows had infinite "y" observations that were removed.'

  # -Inf and Inf among x observations
  expect_warning(res <- row_wilcoxon_paired(c(-Inf,5,7,Inf), 1:4), wrnX, all=TRUE)
  expect_equal(res$obs.x, 2)
  expect_equal(res$obs.y, 4)
  expect_equal(res$obs.paired, 2)

  # -Inf and Inf among y observations
  expect_warning(res <- row_wilcoxon_paired(1:4, c(-Inf,5,7,Inf)), wrnY, all=TRUE)
  expect_equal(res$obs.x, 4)
  expect_equal(res$obs.y, 2)
  expect_equal(res$obs.paired, 2)

  # -Inf and Inf among both
  expect_warning(res <- row_wilcoxon_paired(c(1,2,3,Inf), c(-Inf,7,6,5)), wrnX)
  expect_warning(res <- row_wilcoxon_paired(c(1,2,3,Inf), c(-Inf,7,6,5)), wrnY)
  expect_equal(res$obs.x, 3)
  expect_equal(res$obs.y, 3)
  expect_equal(res$obs.paired, 2)
})

test_that("warning when x=mu values are removed", {
  wrn <- 'row_wilcoxon_paired: 1 of the rows had observations with "x-y" equal "mu" that were removed.'

  # a few x = mu
  expect_warning(res <- row_wilcoxon_paired(1:4, c(0,2,4,3), mu=1, exact=FALSE), wrn, all=TRUE)
  expect_equal(res$obs.x, 4)
  expect_equal(res$obs.y, 4)
  expect_equal(res$obs.paired, 2)
})

test_that("warning when x=mu values are removed and exact cannot be used", {
  wrn1 <- 'row_wilcoxon_paired: 1 of the rows had observations with "x-y" equal "mu" that were removed.'
  wrn2 <- 'row_wilcoxon_paired: 1 of the rows had zeroes: cannot compute exact p-values with zeroes.'

  # a few x = mu
  expect_warning(res <- row_wilcoxon_paired(1:4, c(0,2,4,3), mu=1, exact=TRUE), wrn1)
  expect_warning(res <- row_wilcoxon_paired(1:4, c(0,2,4,3), mu=1, exact=TRUE), wrn2)
  expect_equal(res$obs.x, 4)
  expect_equal(res$obs.y, 4)
  expect_equal(res$obs.paired, 2)
  expect_false(res$exact)

  # no second warnng when exact=FALSE
  expect_warning(res <- row_wilcoxon_paired(1:4, c(0,2,4,3), mu=1, exact=FALSE), wrn1)
  expect_false(res$exact)
})

test_that("warning when ties are present", {
  wrn <- 'row_wilcoxon_paired: 1 of the rows had ties: cannot compute exact p-values with ties.'

  # warning when exact=TRUE
  expect_warning(res <- row_wilcoxon_paired(c(3,3,1,2), c(2,2,5,-2), exact=TRUE), wrn, all=TRUE)
  expect_equal(res$obs.x, 4)
  expect_equal(res$obs.y, 4)
  expect_equal(res$obs.paired, 4)
  expect_false(res$exact)

  # no warning when exact=FALSE
  expect_warning(res <- row_wilcoxon_paired(c(3,3,1,2), c(2,2,5,-2), exact=FALSE), NA)
  expect_false(res$exact)
})

test_that("only warning about zero when both zeroes and ties are present", {
  wrn1 <- 'row_wilcoxon_paired: 1 of the rows had observations with "x-y" equal "mu" that were removed.'
  wrn2 <- 'row_wilcoxon_paired: 1 of the rows had zeroes: cannot compute exact p-values with zeroes.'

  # no warning when exact=FALSE
  expect_warning(res <- row_wilcoxon_paired(1:5, c(1,1,2,9,9), exact=TRUE), paste0(wrn1, "|", wrn2), all=TRUE)
  expect_false(res$exact)
})

