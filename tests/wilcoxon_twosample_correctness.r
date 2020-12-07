library(matrixTests)

#--- functions -----------------------------------------------------------------

base_wilcoxon_twosample <- function(mat1, mat2, null=0, alternative="two.sided", exact=NA, correct=TRUE) {
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  if(length(alternative)==1) alternative <- rep(alternative, nrow(mat1))
  if(length(null)==1) null <- rep(null, nrow(mat1))
  if(length(exact)==1) exact <- rep(exact, nrow(mat1))
  if(length(correct)==1) correct <- rep(correct, nrow(mat1))

  nx <- ny <- nxy <- stat <- p <- m0 <- numeric(nrow(mat1))
  al <- character(nrow(mat1))
  ext <- cor <- logical(nrow(mat1))
  for(i in 1:nrow(mat1)) {
    ex  <- if(is.na(exact[i])) NULL else exact[i]

    res <- wilcox.test(mat1[i,], mat2[i,], alternative=alternative[i], mu=null[i], exact=ex, correct=correct[i])

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
             location.null=m0, alternative=al, exact=ext, corrected=cor,
             stringsAsFactors=FALSE
             )
}


#--- montecarlo ----------------------------------------------------------------

# 2 samples in one group and 4 in another
x <- matrix(rnorm(2000), ncol=2)
y <- matrix(rnorm(4000), ncol=4)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor)
res2 <- row_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor)
stopifnot(all.equal(res1, res2))

# 51 samples in one group and 51 in another
x <- matrix(rnorm(51000), ncol=51)
y <- matrix(rnorm(51000), ncol=51)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor)
res2 <- row_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor)
stopifnot(all.equal(res1, res2))

# 4 samples in one group and 51 in another
x <- matrix(rnorm(4000), ncol=4)
y <- matrix(rnorm(51000), ncol=51)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor)
res2 <- row_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor)
stopifnot(all.equal(res1, res2))

# 6 samples per group with ties and zeroes
x <- matrix(round(runif(6000, -15, 15)), ncol=6)
y <- matrix(round(runif(6000, -15, 15)), ncol=6)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
mus[sample(length(mus), length(mus)/8)] <- 0
mus[sample(length(mus), length(mus)/8)] <- 1
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor))
stopifnot(all.equal(res1, res2))

# 60 samples per group with ties and zeroes
x <- matrix(round(runif(60000, -15, 15)), ncol=60)
y <- matrix(round(runif(60000, -15, 15)), ncol=60)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
mus[sample(length(mus), length(mus)/8)] <- 0
mus[sample(length(mus), length(mus)/8)] <- 1
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor))
stopifnot(all.equal(res1, res2))


# 6 samples in one group and 60 in another with ties and zeroes
# 60 samples per group with ties and zeroes
x <- matrix(round(runif(6000, -15, 15)), ncol=6)
y <- matrix(round(runif(60000, -15, 15)), ncol=60)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- runif(nrow(x), -1, 1)
mus[sample(length(mus), length(mus)/8)] <- 0
mus[sample(length(mus), length(mus)/8)] <- 1
cor <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, mus, alts, exact=ext, correct=cor))
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
x    <- matrix(x, nrow=nrow(pars), ncol=length(x), byrow=TRUE)
y    <- c(100000000000005, 100000000000006, 100000000000009)
y    <- matrix(y, nrow=nrow(pars), ncol=length(y), byrow=TRUE)
res1 <- base_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))

# small numbers
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0.00000000000001)
x    <- matrix(x, nrow=nrow(pars), ncol=length(x), byrow=TRUE)
y    <- c(0.00000000000005, 0.00000000000006, 0.00000000000000)
y    <- matrix(y, nrow=nrow(pars), ncol=length(y), byrow=TRUE)
res1 <- base_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))

# x has infinity
x <- c(Inf, rnorm(4))
y <- rnorm(5)
res1 <- base_wilcoxon_twosample(x, y)
res2 <- row_wilcoxon_twosample(x, y)
stopifnot(all.equal(res1, res2))

# x has both positive and negative infinities
x <- c(Inf, rnorm(4), -Inf)
y <- rnorm(5)
res1 <- base_wilcoxon_twosample(x, y)
res2 <- row_wilcoxon_twosample(x, y)
stopifnot(all.equal(res1, res2))

# y has infinity
x <- rnorm(5)
y <- c(Inf, rnorm(4))
res1 <- base_wilcoxon_twosample(x, y)
res2 <- row_wilcoxon_twosample(x, y)
stopifnot(all.equal(res1, res2))

# y has positive and negative infinities
x <- rnorm(5)
y <- c(Inf, rnorm(4), -Inf)
res1 <- base_wilcoxon_twosample(x, y)
res2 <- row_wilcoxon_twosample(x, y)
stopifnot(all.equal(res1, res2))

# both and y have various infinities
x <- c(-Inf, rnorm(4), Inf, Inf)
y <- c(Inf, rnorm(4), -Inf)
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y))
stopifnot(all.equal(res1, res2))

# infinities are treated as highest ranks
x1 <- c(-Inf, 1:10, Inf, Inf)
y1 <- c(Inf, 31:40, -Inf)
x2 <- c(-100, 1:10, 100, 100)
y2 <- c(100, 31:40, -100)
res1 <- suppressWarnings(row_wilcoxon_twosample(x1, y1))
res2 <- suppressWarnings(row_wilcoxon_twosample(x2, y2))
stopifnot(all.equal(res1, res2))


#--- minimal sample size -------------------------------------------------------

# single number in each group
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rnorm(nrow(pars)), ncol=1)
y    <- matrix(rnorm(nrow(pars)), ncol=1)
res1 <- base_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$obs.x, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.y, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.tot, rep(2, nrow(x))))

# 3 observation in both, but only one non NA
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- cbind(matrix(rnorm(nrow(pars)), ncol=1), NA, NaN)
y    <- cbind(NaN, NA, matrix(rnorm(nrow(pars)), ncol=1))
res1 <- base_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3])
res2 <- row_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3])
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$obs.x, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.y, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.tot, rep(2, nrow(x))))

# single equal number in each group
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rep(1, nrow(pars)), ncol=1)
y    <- matrix(rep(1, nrow(pars)), ncol=1)
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3]))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3]))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$obs.x, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.y, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.tot, rep(2, nrow(x))))

# three numbers in both groups, all equal
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rnorm(nrow(pars)), ncol=1)
x    <- cbind(x,x,x)
y    <- x
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3]))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, alternative=pars[,1], exact=pars[,2], correct=pars[,3]))
stopifnot(all.equal(res1, res2))

# three numbers all equal after subtracting null from x
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rnorm(nrow(pars)), ncol=1)
x    <- cbind(x, x, x)
y    <- x - 1
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y, alternative=pars[,1], null=1, exact=pars[,2], correct=pars[,3]))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, alternative=pars[,1], null=1, exact=pars[,2], correct=pars[,3]))
stopifnot(all.equal(res1, res2))


#--- parameter edge cases ------------------------------------------------------

# various corner cases with NAs
alt <- c("l", "t", "g")
mus <- c(-10000, 0, 10000)  # TODO: decide if Inf is reasonable null for this test
ext <- c(TRUE, FALSE, NA)
cor <- c(TRUE, FALSE)
pars <- expand.grid(mus, alt, ext, cor, stringsAsFactors=FALSE)
x <- matrix(round(runif(10*nrow(pars), -25, 25)), ncol=10)
x[sample(length(x), nrow(pars)*2)] <- NA
y <- matrix(round(runif(8*nrow(pars), -25, 25)), ncol=8)
y[sample(length(y), nrow(pars)*2)] <- NA
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y, pars[,1], pars[,2], pars[,3], pars[,4]))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, pars[,1], pars[,2], pars[,3], pars[,4]))
stopifnot(all.equal(res1, res2))


#--- interactions with exact and correct ---------------------------------------

# 49 samples (< 50: exact = TRUE)
x <- rnorm(49)
y <- rnorm(49)
res1 <- base_wilcoxon_twosample(x, y)
res2 <- row_wilcoxon_twosample(x, y, exact=NA)
res3 <- row_wilcoxon_twosample(x, y, exact=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$obs.x, 49))
stopifnot(all.equal(res2$obs.y, 49))
stopifnot(all.equal(res2$obs.tot, 49*2))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 50 samples with one NA (< 50: exact = TRUE)
x <- c(rnorm(49), NA)
y <- c(NA, rnorm(49))
res1 <- base_wilcoxon_twosample(x, y)
res2 <- row_wilcoxon_twosample(x, y, exact=NA)
res3 <- row_wilcoxon_twosample(x, y, exact=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$obs.x, 49))
stopifnot(all.equal(res2$obs.y, 49))
stopifnot(all.equal(res2$obs.tot, 49*2))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 50 samples (>= 50: exact = FALSE)
x <- rnorm(50)
y <- rnorm(50)
res1 <- base_wilcoxon_twosample(x, y)
res2 <- row_wilcoxon_twosample(x, y, exact=NA)
res3 <- row_wilcoxon_twosample(x, y, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.x, 50))
stopifnot(all.equal(res2$obs.y, 50))
stopifnot(all.equal(res2$obs.tot, 100))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 51 samples with one NA (>= 50: exact = FALSE)
x <- c(NA, rnorm(50))
y <- c(rnorm(50), NA)
res1 <- base_wilcoxon_twosample(x, y)
res2 <- row_wilcoxon_twosample(x, y, exact=NA)
res3 <- row_wilcoxon_twosample(x, y, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.x, 50))
stopifnot(all.equal(res2$obs.y, 50))
stopifnot(all.equal(res2$obs.tot, 100))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# One small group and another with 50 values (>= 50: exact = FALSE)
x <- rnorm(5)
y <- rnorm(50)
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$obs.x, 5))
stopifnot(all.equal(res2$obs.y, 50))
stopifnot(all.equal(res2$obs.tot, 55))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 10 samples with one tie after subtracting null=1 (ties: exact = FALSE)
x <- c(rnorm(9), 2)
y <- c(rnorm(9), 1)
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y, null=1))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, null=1, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_twosample(x, y, null=1, exact=TRUE))
res4 <- row_wilcoxon_twosample(x, y, null=1, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res3$exact, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# 105 samples with 5 ties (ties: exact = FALSE)
x <- c(rnorm(50), 2, 2, 2, 2, 2)
y <- c(rnorm(50))
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=NA))
res3 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=TRUE))
res4 <- row_wilcoxon_twosample(x, y, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res3$exact, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res3, res4))

# big sample size so exact=FALSE (exact = FALSE: correct can be specified)
x <- rnorm(50)
y <- rnorm(50)
res1 <- base_wilcoxon_twosample(x, y, correct=FALSE)
res2 <- row_wilcoxon_twosample(x, y, exact=NA, correct=FALSE)
res3 <- row_wilcoxon_twosample(x, y, exact=FALSE, correct=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# big sample size, exact is turned on but exact=FALSE because of ties (exact = FALSE: correct can be specified)
x <- c(rnorm(50), 2, 2)
y <- c(rnorm(40), 4, 4)
res1 <- suppressWarnings(base_wilcoxon_twosample(x, y, exact=TRUE, correct=FALSE))
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=NA, correct=FALSE))
res3 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=TRUE, correct=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# small sample size but exact is turned off (exact = FALSE: correct can be specified)
x <- rnorm(20)
y <- rnorm(20)
res1 <- base_wilcoxon_twosample(x, y, exact=FALSE, correct=FALSE)
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=FALSE, correct=FALSE))
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

# big sample size with exact=TRUE (exact = TRUE: correct turned to FALSE)
x <- rnorm(60)
y <- rnorm(60)
res1 <- base_wilcoxon_twosample(x, y, exact=TRUE, correct=TRUE)
res2 <- row_wilcoxon_twosample(x, y, exact=TRUE, correct=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

# big sample size but NAs make it less than 50 so exact is turned on (exact = TRUE: correct turned to FALSE)
x <- c(rnorm(49), NA, NaN)
y <- c(rnorm(49), NA, NaN)
res1 <- base_wilcoxon_twosample(x, y, correct=TRUE)
res2 <- row_wilcoxon_twosample(x, y, correct=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

# small sample size and exact is left as NA (exact = TRUE: correct turned to FALSE)
x <- rnorm(20)
y <- rnorm(15)
res1 <- base_wilcoxon_twosample(x, y, exact=NA, correct=TRUE)
res2 <- suppressWarnings(row_wilcoxon_twosample(x, y, exact=NA, correct=TRUE))
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res2$corrected, FALSE))
stopifnot(all.equal(res1, res2))

