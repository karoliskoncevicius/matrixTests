library(matrixTests)

#--- functions -----------------------------------------------------------------

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


#--- montecarlo ----------------------------------------------------------------

# 4 paired samples
x <- matrix(rnorm(4000), ncol=4)
y <- matrix(rnorm(4000), ncol=4)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
cor  <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext  <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_wilcoxon_paired(x, y, alts, mus, exact=ext, correct=cor)
res2 <- row_wilcoxon_paired(x, y, alts, mus, exact=ext, correct=cor)
stopifnot(all.equal(res1, res2))

# 60 paired samples
x <- matrix(rnorm(60000), ncol=60)
y <- matrix(rnorm(60000), ncol=60)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
cor  <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext  <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_wilcoxon_paired(x, y, alts, mus, exact=ext, correct=cor)
res2 <- row_wilcoxon_paired(x, y, alts, mus, exact=ext, correct=cor)
stopifnot(all.equal(res1, res2))


# 6 paired samples per group with ties and zeroes
x <- matrix(round(runif(6000, -10, 10)), ncol=6)
y <- matrix(round(runif(6000, -10, 10)), ncol=6)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
mus[sample(length(mus), length(mus)/8)] <- 0
mus[sample(length(mus), length(mus)/8)] <- 1
cor  <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext  <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, alts, mus, exact=ext, correct=cor))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, alts, mus, exact=ext, correct=cor))
stopifnot(all.equal(res1, res2))

# 60 paired samples per group with ties and zeroes
x <- matrix(round(runif(60000, -15, 15)), ncol=60)
y <- matrix(round(runif(60000, -15, 15)), ncol=60)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
mus[sample(length(mus), length(mus)/8)] <- 0
mus[sample(length(mus), length(mus)/8)] <- 1
cor  <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext  <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, alts, mus, exact=ext, correct=cor))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, alts, mus, exact=ext, correct=cor))
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
pars  <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
y <- c(100000000000002, 100000000000003, 100000000000000, 100000000000004)
x <- matrix(x, nrow=nrow(pars), ncol=length(x), byrow=TRUE)
y <- matrix(y, nrow=nrow(pars), ncol=length(y), byrow=TRUE)
res1 <- base_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))

# small numbers
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0.00000000000001)
y <- c(0.00000000000002, 0.00000000000003, 0.00000000000000, 0.00000000000005)
x <- matrix(x, nrow=nrow(pars), ncol=length(x), byrow=TRUE)
y <- matrix(y, nrow=nrow(pars), ncol=length(y), byrow=TRUE)
res1 <- base_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))

# x has infinity
x <- c(Inf, rnorm(4))
y <- rnorm(5)
res1 <- base_wilcoxon_paired(x, y)
res2 <- row_wilcoxon_paired(x, y)
stopifnot(all.equal(res1, res2))

# x has both positive and negative infinities
x <- c(Inf, rnorm(4), -Inf)
y <- rnorm(6)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y))
stopifnot(all.equal(res1, res2))

# y has infinity
x <- rnorm(5)
y <- c(Inf, rnorm(4))
res1 <- base_wilcoxon_paired(x, y)
res2 <- row_wilcoxon_paired(x, y)
stopifnot(all.equal(res1, res2))

# y has positive and negative infinities
x <- rnorm(6)
y <- c(Inf, rnorm(4), -Inf)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y))
stopifnot(all.equal(res1, res2))

# both and y have various infinities
# NOTE: base implementation doesn't work when both Inf with the same sign is present in both groups
x <- c(-Inf, rnorm(4), Inf, Inf)
y <- c(Inf, Inf, rnorm(4), -Inf)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y))
stopifnot(all.equal(res1, res2))

# infinities are treated as highest ranks
# NOTE: this scenario is a bit more complex because the ranks are calculated after taking differences.
x1 <- c(-Inf, 1:10, Inf, Inf)
y1 <- c(-Inf, Inf, 31:40, -Inf)
x2 <- c(-100, 0, 2:10, 200, 100)
y2 <- c(-100, 200, 31:39, 0, -100)
res1 <- suppressWarnings(row_wilcoxon_paired(x1, y1))
res2 <- suppressWarnings(row_wilcoxon_paired(x2, y2))
stopifnot(all.equal(res1, res2))


#--- minimal sample size -------------------------------------------------------

# single paired number
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rnorm(nrow(pars)), ncol=1)
y    <- matrix(rnorm(nrow(pars)), ncol=1)
res1 <- base_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$obs.paired, rep(1, nrow(x))))

# 3 paired numbers but only one common, others are paired with NAs
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- cbind(matrix(rnorm(nrow(pars)*2), ncol=2), NA)
y    <- cbind(NA, matrix(rnorm(nrow(pars)*2), ncol=2))
res1 <- base_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$obs.paired, rep(1, nrow(x))))

# three paired numbers all equal
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x  <- matrix(rnorm(nrow(pars)), nrow=nrow(pars), ncol=3)
y  <- matrix(rnorm(nrow(pars)), nrow=nrow(pars), ncol=3)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3]))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3]))
stopifnot(all.equal(res1, res2))

# three numbers all differences equal
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x <- matrix(rnorm(nrow(pars)*3), ncol=3)
y <- x + 0.2
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3]))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, pars[,1], exact=pars[,2], correct=pars[,3]))
stopifnot(all.equal(res1, res2))

# three numbers, two equal to mu after subtraction
# NOTE: seems to be numerically unstable
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x <- matrix(rnorm(nrow(pars)*3), ncol=3)
y <- matrix(rnorm(nrow(pars)*3), ncol=3)
y[,1:2] <- x[,1:2]-1
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, pars[,1], 1, exact=pars[,2], correct=pars[,3]))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, pars[,1], 1, exact=pars[,2], correct=pars[,3]))
stopifnot(all.equal(res1, res2))


#--- parameter edge cases ------------------------------------------------------

# various corner cases with NAs
alt  <- c("l", "t", "g")
mus  <- c(-10000, 0, 10000)  # TODO: decide if Inf is reasonable mu for this test
ext  <- c(TRUE, FALSE, NA)
cor  <- c(TRUE, FALSE)
pars <- expand.grid(alt, mus, ext, cor, stringsAsFactors=FALSE)
x <- matrix(round(runif(10*nrow(pars), -15, 15)), ncol=10)
x[sample(length(x), nrow(pars)*2)] <- NA
y <- matrix(round(runif(10*nrow(pars), -15, 15)), ncol=10)
y[sample(length(y), nrow(pars)*2)] <- NA
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, pars[,1], pars[,2], pars[,3], pars[,4]))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, pars[,1], pars[,2], pars[,3], pars[,4]))
stopifnot(all.equal(res1, res2))


#--- interactions with exact and correct ---------------------------------------

# 49 samples (< 50: exact = TRUE)
x <- rnorm(49)
y <- rnorm(49)
res1 <- base_wilcoxon_paired(x, y)
res2 <- row_wilcoxon_paired(x, y, exact=NA)
res3 <- row_wilcoxon_paired(x, y, exact=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$obs.paired, 49))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 50 samples with one NA (< 50: exact = TRUE)
x <- c(NA, rnorm(49))
y <- rnorm(50)
res1 <- base_wilcoxon_paired(x, y)
res2 <- row_wilcoxon_paired(x, y, exact=NA)
res3 <- row_wilcoxon_paired(x, y, exact=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$obs.paired, 49))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 51 samples with NA and NaN (< 50: exact = TRUE)
x <- c(rnorm(50), NA)
y <- c(NaN, rnorm(50))
res1 <- base_wilcoxon_paired(x, y)
res2 <- row_wilcoxon_paired(x, y, exact=NA)
res3 <- row_wilcoxon_paired(x, y, exact=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$obs.paired, 49))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 50 samples (>= 50: exact = FALSE)
x <- rnorm(50)
y <- rnorm(50)
res1 <- base_wilcoxon_paired(x, y)
res2 <- row_wilcoxon_paired(x, y, exact=NA)
res3 <- row_wilcoxon_paired(x, y, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.paired, 50))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 51 samples with one NA (>= 50: exact = FALSE)
x <- rnorm(51)
y <- c(rnorm(50), NA)
res1 <- base_wilcoxon_paired(x, y)
res2 <- row_wilcoxon_paired(x, y, exact=NA)
res3 <- row_wilcoxon_paired(x, y, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.paired, 50))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 50 samples with 1 value equal to mu (>= 50: exact = FALSE)
# NOTE: even thou 49 paired observations are left, base and matrixTests still treats it as 50
x <- c(rnorm(49),1.2)
y <- c(rnorm(49),1.2)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.paired, 49))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 10 samples with one tie (ties: exact = FALSE)
x <- c(rnorm(8), 2, 2)
y <- c(rnorm(8), 1, 1)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=TRUE))
res4 <- row_wilcoxon_paired(x, y, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 105 samples with 5 ties (ties: exact = FALSE)
x <- c(rnorm(100), 2, 2, 2, 2, 2)
y <- c(rnorm(100), 1, 1, 1, 1, 1)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=TRUE))
res4 <- row_wilcoxon_paired(x, y, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 10 samples with 1 zero (zeroes: exact = FALSE)
x <- c(rnorm(9), 2)
y <- c(rnorm(9), 1)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE))
res4 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.paired, 9))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 105 samples with 5 zeroes (zeroes: exact = FALSE)
x <- c(rnorm(100), 2, 2, 2, 2, 2)
y <- c(rnorm(100), 1, 1, 1, 1, 1)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE))
res4 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.paired, 100))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 11 samples with one zero and one tie (zeroes or ties: exact = FALSE)
x <- c(rnorm(8), 2, 2, 3)
y <- c(rnorm(8), 4, 4, 2)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE))
res4 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.paired, 10))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 110 samples with 5 zeroes and 5 ties
x <- c(rnorm(100), 2,2,2,2,2, 3,3,3,3,3)
y <- c(rnorm(100), 4,4,4,4,4, 2,2,2,2,2)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE))
res4 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.paired, 105))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))


# big sample size so exact=FALSE (exact = FALSE: correct can be specified)
x <- rnorm(100)
y <- rnorm(100)
res1 <- base_wilcoxon_paired(x, y, correct=FALSE)
res2 <- row_wilcoxon_paired(x, y, exact=NA, correct=FALSE)
res3 <- row_wilcoxon_paired(x, y, exact=FALSE, correct=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# big sample size, exact is turned on but exact=FALSE because of ties (exact = FALSE: correct can be specified)
x <- c(rnorm(100), 2, 2, 3)
y <- c(rnorm(100), 4, 4, 2)
res1 <- suppressWarnings(base_wilcoxon_paired(x, y, mu=1, exact=TRUE, correct=FALSE))
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=NA, correct=FALSE))
res3 <- suppressWarnings(row_wilcoxon_paired(x, y, mu=1, exact=TRUE, correct=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# small sample size but exact is turned off (exact = FALSE: correct can be specified)
x <- rnorm(20)
y <- rnorm(20)
res1 <- base_wilcoxon_paired(x, y, exact=FALSE, correct=FALSE)
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=FALSE, correct=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))


# big sample size with exact=TRUE (exact = TRUE: correct turned to FALSE)
x <- rnorm(100)
y <- rnorm(100)
res1 <- base_wilcoxon_paired(x, y, exact=TRUE, correct=TRUE)
res2 <- row_wilcoxon_paired(x, y, exact=TRUE, correct=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

# big sample size but NAs make it less than 50 so exact is turned on (exact = TRUE: correct turned to FALSE)
x <- rnorm(51)
y <- c(rnorm(49), NaN, NA)
res1 <- base_wilcoxon_paired(x, y, correct=TRUE)
res2 <- row_wilcoxon_paired(x, y, correct=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

# small sample size and exact is left as NA (exact = TRUE: correct turned to FALSE)
x <- rnorm(20)
y <- rnorm(20)
res1 <- base_wilcoxon_paired(x, y, exact=NA, correct=TRUE)
res2 <- suppressWarnings(row_wilcoxon_paired(x, y, exact=NA, correct=TRUE))
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

