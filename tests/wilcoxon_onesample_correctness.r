library(matrixTests)

#--- functions -----------------------------------------------------------------

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
    vec <- mat[i,][!is.na(mat[i,])]
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


#--- montecarlo ----------------------------------------------------------------

# 3 observations
x <- matrix(rnorm(3000), ncol=3)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_wilcoxon_onesample(x, alts, mus, exact=ext, correct=cor)
res2 <- row_wilcoxon_onesample(x, alts, mus, exact=ext, correct=cor)
stopifnot(all.equal(res1, res2))

# 60 observations
x <- matrix(rnorm(60000), ncol=60)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_wilcoxon_onesample(x, alts, mus, exact=ext, correct=cor)
res2 <- row_wilcoxon_onesample(x, alts, mus, exact=ext, correct=cor)
stopifnot(all.equal(res1, res2))

# 10 observations with ties and zeroes
x <- matrix(round(runif(10000, -15, 15)), ncol=10)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
mus[sample(length(mus), length(mus)/8)] <- 0
mus[sample(length(mus), length(mus)/8)] <- 1
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- suppressWarnings(base_wilcoxon_onesample(x, alts, mus, exact=ext, correct=cor))
res2 <- suppressWarnings(row_wilcoxon_onesample(x, alts, mus, exact=ext, correct=cor))
stopifnot(all.equal(res1, res2))

# 60 observations with ties and zeroes
x <- matrix(round(runif(60000, -15, 15)), ncol=60)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
mus[sample(length(mus), length(mus)/8)] <- 0
mus[sample(length(mus), length(mus)/8)] <- 1
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- suppressWarnings(base_wilcoxon_onesample(x, alts, mus, exact=ext, correct=cor))
res2 <- suppressWarnings(row_wilcoxon_onesample(x, alts, mus, exact=ext, correct=cor))
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
x    <- matrix(x, nrow=nrow(pars), ncol=length(x), byrow=TRUE)
res1 <- base_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))

# small numbers
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0.00000000000001)
x    <- matrix(x, nrow=nrow(pars), ncol=length(x), byrow=TRUE)
res1 <- base_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))

# has infinity
x <- c(Inf, rnorm(4))
res1 <- base_wilcoxon_onesample(x)
res2 <- row_wilcoxon_onesample(x)
stopifnot(all.equal(res1, res2))

# infinity is treated as highes rank
x1 <- c(Inf, 1:4)
x2 <- c(100, 1:4)
res1 <- row_wilcoxon_onesample(x1)
res2 <- row_wilcoxon_onesample(x2)
stopifnot(all.equal(res1, res2))

# has negative and positive infinities
x <- c(Inf, rnorm(4), -Inf)
res1 <- suppressWarnings(base_wilcoxon_onesample(x))
res2 <- suppressWarnings(row_wilcoxon_onesample(x))
stopifnot(all.equal(res1, res2))

# all infinities are treated as highest ranks
x1 <- c(-Inf, -Inf, 1:4, Inf, Inf)
x2 <- c(-100, -100, 1:4, 100, 100)
res1 <- suppressWarnings(row_wilcoxon_onesample(x1))
res2 <- suppressWarnings(row_wilcoxon_onesample(x2))
stopifnot(all.equal(res1, res2))


#--- minimal sample size -------------------------------------------------------

# single number
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rnorm(nrow(pars)), ncol=1)
res1 <- base_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))

# single number with NAs
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- cbind(rnorm(nrow(pars)), NA, NA)
res1 <- base_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))

# three numbers all equal
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rnorm(nrow(pars)), nrow=nrow(pars), ncol=3)
res1 <- suppressWarnings(base_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3]))
res2 <- suppressWarnings(row_wilcoxon_onesample(x, pars[,1], exact=pars[,2], correct=pars[,3]))
stopifnot(all.equal(res1, res2))

# three numbers, two equal to null
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x  <- cbind(matrix(rnorm(nrow(pars)), ncol=1), 1, 1)
res1 <- suppressWarnings(base_wilcoxon_onesample(x, pars[,1], 1, exact=pars[,2], correct=pars[,3]))
res2 <- suppressWarnings(row_wilcoxon_onesample(x, pars[,1], 1, exact=pars[,2], correct=pars[,3]))
stopifnot(all.equal(res1, res2))


#--- parameter edge cases ------------------------------------------------------

# various corner cases with NAs
alt <- c("l", "t", "g")
mus <- c(-10000, 0, 10000)  # TODO: decide if Inf is reasonable mu for this test
ext <- c(TRUE, FALSE, NA)
cor <- c(TRUE, FALSE)
pars <- expand.grid(alt, mus, ext, cor, stringsAsFactors=FALSE)
x <- matrix(round(runif(10*nrow(pars), -15, 15)), ncol=10)
x[sample(length(x), nrow(pars)*2)] <- NA
res1 <- suppressWarnings(base_wilcoxon_onesample(x, pars[,1], pars[,2], pars[,3], pars[,4]))
res2 <- suppressWarnings(row_wilcoxon_onesample(x, pars[,1], pars[,2], pars[,3], pars[,4]))
stopifnot(all.equal(res1, res2))

# null exactly equal to the median
res1 <- suppressWarnings(base_wilcoxon_onesample(c(1,2,3), mu=2))
res2 <- suppressWarnings(row_wilcoxon_onesample(c(1,2,3), mu=2))
stopifnot(all.equal(res2$pvalue, 1))
stopifnot(all.equal(res1, res2))


#--- interactions with exact and correct ---------------------------------------

# 49 samples (< 50: exact = TRUE)
x <- rnorm(49)
res1 <- base_wilcoxon_onesample(x, exact=TRUE)
res2 <- row_wilcoxon_onesample(x, exact=TRUE)
res3 <- row_wilcoxon_onesample(x, exact=NA)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 50 samples with one NA (< 50: exact = TRUE)
x  <- c(rnorm(49), NA)
res1 <- base_wilcoxon_onesample(x, exact=TRUE)
res2 <- row_wilcoxon_onesample(x, exact=TRUE)
res3 <- row_wilcoxon_onesample(x, exact=NA)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 50 samples (>= 50: exact = FALSE)
x <- rnorm(50)
res1 <- base_wilcoxon_onesample(x, exact=FALSE)
res2 <- row_wilcoxon_onesample(x, exact=FALSE)
res3 <- row_wilcoxon_onesample(x, exact=NA)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 51 samples with one NA (>= 50: exact = FALSE)
x <- c(rnorm(50), NA)
res1 <- base_wilcoxon_onesample(x, exact=FALSE)
res2 <- row_wilcoxon_onesample(x, exact=FALSE)
res3 <- row_wilcoxon_onesample(x, exact=NA)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 11 samples with one tie (ties: exact = FALSE)
x <- c(1.2,rnorm(9),1.2)
res1 <- base_wilcoxon_onesample(x, exact=FALSE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_onesample(x, exact=TRUE))
res4 <- row_wilcoxon_onesample(x, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res3$exact, FALSE))
stopifnot(all.equal(res2$corrected, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 105 samples with 5 ties (ties: exact = FALSE)
x <- c(1.2,rnorm(100),rep(1.2,4))
res1 <- base_wilcoxon_onesample(x, exact=FALSE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_onesample(x, exact=TRUE))
res4 <- row_wilcoxon_onesample(x, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res3$exact, FALSE))
stopifnot(all.equal(res2$corrected, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 11 samples with one zero (zero: exact = FALSE)
x <- rnorm(11)
res1 <- base_wilcoxon_onesample(x, mu=x[1], exact=FALSE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=NA))
res3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=TRUE))
res4 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res3$exact, FALSE))
stopifnot(all.equal(res2$corrected, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 105 samples with 5 zeroes (zero: exact = FALSE)
x <- c(1, rnorm(100), 1, 1, 1, 1)
res1 <- base_wilcoxon_onesample(x, mu=x[1], exact=FALSE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=NA))
res3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=TRUE))
res4 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res3$exact, FALSE))
stopifnot(all.equal(res2$corrected, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 11 samples with one zero and one tie (zero or tie: exact = FALSE)
x <- c(rnorm(10),1.2,1.2)
res1 <- base_wilcoxon_onesample(x, mu=x[1], exact=FALSE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=NA))
res3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=TRUE))
res4 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res3$exact, FALSE))
stopifnot(all.equal(res2$corrected, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 110 samples with 5 zeroes and 5 ties (zero or tie: exact = FALSE)
x <- c(1,1,1,1,1,rnorm(100),1.2,1.2,1.2,1.2,1.2)
res1 <- base_wilcoxon_onesample(x, mu=x[1], exact=FALSE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=NA))
res3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=TRUE))
res4 <- suppressWarnings(row_wilcoxon_onesample(x, mu=x[1], exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res3$exact, FALSE))
stopifnot(all.equal(res2$corrected, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# big sample size so exact=FALSE (exact = FALSE: correct can be specified)
x <- rnorm(100)
res1 <- base_wilcoxon_onesample(x, exact=FALSE, correct=FALSE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA, correct=FALSE))
res3 <- suppressWarnings(row_wilcoxon_onesample(x, exact=FALSE, correct=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# big sample size, exact is turned on but exact=FALSE because of ties and zeroes (exact = FALSE: correct can be specified)
x <- c(1.2,1.2,1.2,1.2,1.2,rnorm(100),1,1,1,1,1)
res1 <- suppressWarnings(base_wilcoxon_onesample(x, mu=1, exact=TRUE, correct=FALSE))
res2 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=NA, correct=FALSE))
res3 <- suppressWarnings(row_wilcoxon_onesample(x, mu=1, exact=TRUE, correct=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# small sample size but exact is turned off (exact = FALSE: correct can be specified)
x <- rnorm(20)
res1 <- base_wilcoxon_onesample(x, exact=FALSE, correct=FALSE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=FALSE, correct=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

# big sample size with exact=TRUE (exact = TRUE: correct turned to FALSE)
x <- rnorm(100)
res1 <- base_wilcoxon_onesample(x, exact=TRUE, correct=TRUE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=TRUE, correct=TRUE))
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

# big sample size but NAs make it less than 50 so exact is turned on (exact = TRUE: correct turned to FALSE)
x <- c(rnorm(49), NaN, NA)
res1 <- base_wilcoxon_onesample(x, exact=TRUE, correct=TRUE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=TRUE, correct=TRUE))
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

# small sample size and exact is left as NA (exact = TRUE: correct turned to FALSE)
x <- rnorm(20)
res1 <- base_wilcoxon_onesample(x, exact=NA, correct=TRUE)
res2 <- suppressWarnings(row_wilcoxon_onesample(x, exact=NA, correct=TRUE))
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

