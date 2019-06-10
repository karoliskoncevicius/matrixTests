context("correctness of wilcoxon_onesample")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_wilcoxon_onesample <- function(mat, alt="two.sided", mu=0, exact=NA, correct=TRUE) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)
  if(length(alt)==1) alt <- rep(alt, nrow(mat))
  if(length(mu)==1) mu <- rep(mu, nrow(mat))
  if(length(exact)==1) exact <- rep(exact, nrow(mat))
  if(length(correct)==1) correct <- rep(correct, nrow(mat))

  nx <- stat <- p <- m0 <- numeric(nrow(mat))
  al <- character(nrow(mat))
  ext <- cor <- logical(nrow(mat))
  for(i in 1:nrow(mat)) {
    vec <- mat[i,][is.finite(mat[i,])]
    ex  <- if(is.na(exact[i])) NULL else exact[i]
    res <- wilcox.test(vec, alternative=alt[i], mu=mu[i], exact=ex, correct=correct[i])

    nx[i]   <- sum((vec-mu[i]) != 0)
    stat[i] <- res$statistic
    p[i]    <- res$p.value
    al[i]   <- res$alternative
    m0[i]   <- res$null.value

    cor[i] <- grepl("continuity correction", res$method)
    ext[i] <- if(is.na(exact[i])) (nx[i] < 50) else exact[i]
    ext[i] <- if(length(unique(abs(vec-m0[i]))) != length(vec)) FALSE else ext[i]
    ext[i] <- if(any((vec-m0[i]) == 0)) FALSE else ext[i]
  }

  data.frame(obs=nx, statistic=stat, pvalue=p, alternative=al, location.null=m0,
             exact=ext, corrected=cor, stringsAsFactors=FALSE
             )
}

################################################################################
########################### TEST ON A RANDOM SAMPLE ############################
################################################################################

############################## SMALL SAMPLE SIZE ###############################

test_that("monte-carlo random testing gives equal results 1", {
  set.seed(14)
  X <- matrix(rnorm(10000), ncol=10)
  X[sample(length(X), nrow(X))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

############################### BIG SAMPLE SIZE ################################

test_that("monte-carlo random testing gives equal results 2", {
  set.seed(14)
  X <- matrix(rnorm(100000), ncol=100)
  X[sample(length(X), 10*nrow(X))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

######################### BORDERLINE SAMPLE SIZE OF 50 #########################

test_that("monte-carlo random testing gives equal results 3", {
  set.seed(14)
  X <- matrix(rnorm(50000), ncol=50)
  X[sample(length(X), 5*nrow(X))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- base_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor)
  t2 <- row_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor)

  expect_equal(t1, t2)
})

#################### SMALL SAMPLE SIZE WITH TIES AND ZEROES ####################

test_that("monte-carlo random testing gives equal results 4", {
  set.seed(14)
  X <- matrix(round(runif(10000, -15, 15)), ncol=10)
  X[sample(length(X), nrow(X))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

##################### BIG SAMPLE SIZE WITH TIES AND ZEROES #####################

test_that("monte-carlo random testing gives equal results 5", {
  set.seed(14)
  X <- matrix(round(runif(100000, -150, 15)), ncol=100)
  X[sample(length(X), 10*nrow(X))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

################# BORDERLINE SAMPLE SIZE WITH TIES AND ZEROES ##################

test_that("monte-carlo random testing gives equal results 6", {
  set.seed(14)
  X <- matrix(round(runif(50000, -150, 15)), ncol=50)
  X[sample(length(X), 5*nrow(X))] <- NA
  alts <- sample(rep(c("t", "g", "l"), length.out=nrow(X)))
  mus  <- runif(nrow(X), -1,1)
  mus[sample(length(mus), length(mus)/8)] <- 0
  mus[sample(length(mus), length(mus)/8)] <- 1
  cor <- sample(rep(c(TRUE, FALSE), length.out=nrow(X)))
  ext <- sample(rep(c(TRUE, FALSE, NA), length.out=nrow(X)))

  t1 <- suppressWarnings(base_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor))
  t2 <- suppressWarnings(row_wilcoxon_onesample(X, alts, mus, exact=ext, correct=cor))

  expect_equal(t1, t2)
})

################################################################################
############################### TEST EDGE CASES ################################
################################################################################

test_that("extreme numbers give equal results", {
  # big numbers
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  vals <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
  mat  <- matrix(vals, nrow=nrow(pars), ncol=length(vals), byrow=TRUE)
  t1 <- base_wilcoxon_onesample(mat, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_onesample(mat, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)

  # small numbers
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  vals <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0.00000000000001)
  mat  <- matrix(vals, nrow=nrow(pars), ncol=length(vals), byrow=TRUE)
  t1 <- base_wilcoxon_onesample(mat, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_onesample(mat, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)
})


test_that("minumum allowed sample sizes give equal results", {
  # single number
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x  <- matrix(rnorm(nrow(pars)), ncol=1)
  t1 <- base_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
  t2 <- row_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
  expect_equal(t1, t2)

  # three numbers all equal
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x  <- matrix(rnorm(nrow(pars)), ncol=1)
  x  <- cbind(x,x,x)
  t1 <- suppressWarnings(base_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3]))
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3]))
  expect_equal(t1, t2)

  # three numbers, two equal to mu
  pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
  x  <- matrix(rnorm(nrow(pars)), ncol=1)
  x  <- cbind(x, 1, 1)
  t1 <- suppressWarnings(base_wilcoxon_onesample(x, pars[,1], 1, exact=pars[,2], correct=pars[,3]))
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, pars[,1], 1, exact=pars[,2], correct=pars[,3]))
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

  t1 <- suppressWarnings(base_wilcoxon_onesample(X, pars[,1], pars[,2], pars[,3], pars[,4]))
  t2 <- suppressWarnings(row_wilcoxon_onesample(X, pars[,1], pars[,2], pars[,3], pars[,4]))
  expect_equal(t1, t2)
})

################################################################################
#################### INTERACTION BETWEEN EXACT AND CORRECT #####################
################################################################################

test_that("small sample size automatically turns exact=TRUE", {
  # 49 samples
  x  <- rnorm(49)
  t1 <- base_wilcoxon_onesample(x, exact=TRUE)
  t2 <- row_wilcoxon_onesample(x, exact=NA)
  t3 <- row_wilcoxon_onesample(x, exact=TRUE)
  expect_true(t2$exact)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 50 samples with one NA
  x  <- rnorm(49)
  x  <- c(x, NA)
  t1 <- base_wilcoxon_onesample(x, exact=TRUE)
  t2 <- row_wilcoxon_onesample(x, exact=NA)
  t3 <- row_wilcoxon_onesample(x, exact=TRUE)
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 51 samples with Inf and -Inf
  x  <- rnorm(49)
  x  <- c(-Inf, x, Inf)
  t1 <- base_wilcoxon_onesample(x, exact=TRUE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, exact=TRUE))
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
})

test_that("big sample size automatically turns exact=FALSE", {
  # 50 samples
  x  <- rnorm(50)
  t1 <- base_wilcoxon_onesample(x, exact=FALSE)
  t2 <- row_wilcoxon_onesample(x, exact=NA)
  t3 <- row_wilcoxon_onesample(x, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 51 samples with one NA
  x  <- rnorm(50)
  x  <- c(x, NA)
  t1 <- base_wilcoxon_onesample(x, exact=FALSE)
  t2 <- row_wilcoxon_onesample(x, exact=NA)
  t3 <- row_wilcoxon_onesample(x, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # 51 samples with Inf and -Inf
  x  <- rnorm(50)
  x  <- c(-Inf, x, Inf)
  t1 <- base_wilcoxon_onesample(x, exact=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
})

test_that("ties automatically turn exact=FALSE, even when exact=TRUE specified", {
  # 11 samples with one tie
  x  <- rnorm(10)
  x  <- c(x, x[5])
  t1 <- base_wilcoxon_onesample(x, exact=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, exact=TRUE))
  t4 <- row_wilcoxon_onesample(x, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
  # 105 samples with 5 ties
  x  <- rnorm(100)
  x  <- c(x, x[1:5])
  t1 <- base_wilcoxon_onesample(x, exact=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, exact=TRUE))
  t4 <- row_wilcoxon_onesample(x, exact=FALSE)
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
})

test_that("zeroes automatically turn exact=FALSE, even when exact=TRUE specified", {
  # 11 samples with one zero
  x  <- rnorm(11)
  t1 <- base_wilcoxon_onesample(x, mu=x[1], exact=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=TRUE))
  t4 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
  # 105 samples with 5 zeroes
  x  <- rnorm(100)
  x  <- c(x, 1, 1, 1, 1, 1)
  t1 <- base_wilcoxon_onesample(x, mu=1, exact=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=TRUE))
  t4 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
})

test_that("mix of zeroes and ties automatically turn exact=FALSE, even when exact=TRUE specified", {
  # 11 samples with one zero and one tie
  x  <- rnorm(11)
  x  <- c(x, x[2])
  t1 <- base_wilcoxon_onesample(x, mu=x[1], exact=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=TRUE))
  t4 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
  # 110 samples with 5 zeroes and 5 ties
  x  <- rnorm(100)
  x  <- c(x, 1, 1, 1, 1, 1)
  x  <- c(x, x[1:5])
  t1 <- base_wilcoxon_onesample(x, mu=1, exact=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=NA))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=TRUE))
  t4 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=FALSE))
  expect_false(t2$exact)
  expect_true(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  expect_equal(t3, t4)
})

test_that("correct=TRUE/FALSE only applies when exact=FALSE", {
  # big sample size so exact=FALSE
  x  <- rnorm(100)
  t1 <- base_wilcoxon_onesample(x, exact=FALSE, correct=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA, correct=FALSE))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, exact=FALSE, correct=FALSE))
  expect_false(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # big sample size, but exact=FALSE becasue of ties and zeroes
  x  <- rnorm(100)
  x  <- c(x, 1, 1, 1, 1, 1)
  x  <- c(x, x[1:5])
  t1 <- suppressWarnings(base_wilcoxon_onesample(x, mu=1, exact=TRUE, correct=FALSE))
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=NA, correct=FALSE))
  t3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=TRUE, correct=FALSE))
  expect_false(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  expect_equal(t2, t3)
  # small sample size but exact is turned off
  x  <- rnorm(20)
  t1 <- base_wilcoxon_onesample(x, exact=FALSE, correct=FALSE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=FALSE, correct=FALSE))
  expect_false(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
})

test_that("correct=TRUE is turned to FALSE when exact=TRUE", {
  # big sample size with exact=TRUE
  x  <- rnorm(100)
  t1 <- base_wilcoxon_onesample(x, exact=TRUE, correct=TRUE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=TRUE, correct=TRUE))
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  # big sample size but NAs make it less than 50 so exact is turned on
  x  <- rnorm(49)
  x  <- c(x, NaN, NA, Inf, -Inf)
  t1 <- base_wilcoxon_onesample(x, exact=TRUE, correct=TRUE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=TRUE, correct=TRUE))
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
  # small sample size and exact is left as NA
  x  <- rnorm(20)
  t1 <- base_wilcoxon_onesample(x, exact=NA, correct=TRUE)
  t2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA, correct=TRUE))
  expect_true(t2$exact)
  expect_false(t2$corrected)
  expect_equal(t1, t2)
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("warnign when row has less than 1 available observations", {
  wrn <- 'row_wilcoxon_onesample: 1 of the rows had less than 1 remaining finite "x" observation.'
  nacolumns <- c("statistic", "pvalue")

  # no observations
  expect_warning(res <- row_wilcoxon_onesample(numeric()), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs, 0)

  # only NA observations
  expect_warning(res <- row_wilcoxon_onesample(c(Inf,NA,NaN,NA)), wrn)
  expect_true(all(is.na(res[,colnames(res) %in% nacolumns])))
  expect_equal(res$obs, 0)
})

test_that("warning when a infinite values are removed", {
  wrn <- 'row_wilcoxon_onesample: 1 of the rows had infinite observations that were removed.'

  # -Inf and Inf among observations
  expect_warning(res <- row_wilcoxon_onesample(c(-Inf,2,3,Inf)), wrn, all=TRUE)
  expect_equal(res$obs, 2)
})

test_that("warning when x=mu values are removed", {
  wrn <- 'row_wilcoxon_onesample: 1 of the rows had observations with "x" equal "mu" that were removed.'

  # a few x = mu
  expect_warning(res <- row_wilcoxon_onesample(c(1,2,3,4,1), mu=1, exact=FALSE), wrn, all=TRUE)
  expect_equal(res$obs, 3)
})

test_that("warning when x=mu values are removed and exact cannot be used", {
  wrn1 <- 'row_wilcoxon_onesample: 1 of the rows had observations with "x" equal "mu" that were removed.'
  wrn2 <- 'row_wilcoxon_onesample: 1 of the rows had zeroes: cannot compute exact p-values with zeroes.'

  # a few x = mu
  expect_warning(res <- row_wilcoxon_onesample(c(1,2,3,4,1), mu=1, exact=TRUE), wrn1)
  expect_warning(res <- row_wilcoxon_onesample(c(1,2,3,4,1), mu=1, exact=TRUE), wrn2)
  expect_equal(res$obs, 3)
  expect_false(res$exact)
  # no second warnng when exact=FALSE
  expect_warning(res <- row_wilcoxon_onesample(c(1,2,3,4,1), mu=1, exact=FALSE), wrn1, all=TRUE)
  expect_equal(res$obs, 3)
  expect_false(res$exact)
})

test_that("warning when ties are present", {
  wrn <- 'row_wilcoxon_onesample: 1 of the rows had ties: cannot compute exact p-values with ties.'

  # warning when exact=TRUE
  expect_warning(res <- row_wilcoxon_onesample(c(1,2,3,4,1), exact=TRUE), wrn)
  expect_equal(res$obs, 5)
  expect_false(res$exact)

  # no warning when exact=FALSE
  expect_warning(res <- row_wilcoxon_onesample(c(1,2,3,4,1), exact=FALSE), NA)
  expect_equal(res$obs, 5)
  expect_false(res$exact)
})

test_that("only warning about zero when both zeroes and ties are present", {
  wrn1 <- 'row_wilcoxon_onesample: 1 of the rows had observations with "x" equal "mu" that were removed.'
  wrn2 <- 'row_wilcoxon_onesample: 1 of the rows had zeroes: cannot compute exact p-values with zeroes.'

  # no warning when exact=FALSE
  expect_warning(res <- row_wilcoxon_onesample(c(1,2,3,4,1), mu=2, exact=TRUE), paste0(wrn1, "|", wrn2), all=TRUE)
  expect_equal(res$obs, 4)
  expect_false(res$exact)
})

