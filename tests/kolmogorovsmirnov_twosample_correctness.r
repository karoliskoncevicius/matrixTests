library(matrixTests)

#--- functions -----------------------------------------------------------------

base_kolmogorovsmirnov_twosample <- function(mat1, mat2, alternative="t", exact=NA) {
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  if(length(alternative)==1) alternative <- rep(alternative, nrow(mat1))
  if(length(exact)==1) exact <- rep(exact, nrow(mat1))

  nx  <- ny <- nxy <- stat <- p <- numeric(nrow(mat1))
  al  <- character(nrow(mat1))
  ext <- logical(nrow(mat1))
  for(i in 1:nrow(mat1)) {
    ex  <- if(is.na(exact[i])) NULL else exact[i]

    res <- ks.test(mat1[i,], mat2[i,], alternative=alternative[i], exact=ex)

    nx[i]   <- sum(!is.na(mat1[i,]))
    ny[i]   <- sum(!is.na(mat2[i,]))
    nxy[i]  <- nx[i] + ny[i]
    stat[i] <- res$statistic
    p[i]    <- res$p.value
    al[i]   <- c(l="less", g="greater", t="two.sided")[alternative[i]]
    ext[i]  <- res$exact
  }

  data.frame(obs.x=nx, obs.y=ny, obs.tot=nxy, statistic=stat, pvalue=p,
             alternative=al, exact=ext, stringsAsFactors=FALSE
             )
}


#--- montecarlo ----------------------------------------------------------------

# 25 samples per group (not exceeding exact cutoff)
x <- matrix(rnorm(2000), ncol=2)
y <- matrix(rnorm(4000), ncol=4)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_kolmogorovsmirnov_twosample(x, y, alts, exact=ext)
res2 <- row_kolmogorovsmirnov_twosample(x, y, alts, exact=ext)
stopifnot(all.equal(res1, res2))

# 25 samples per group with ties (not exceeding exact cutoff)
x <- matrix(round(runif(25000, -15, 15)), ncol=25)
y <- matrix(round(runif(25000, -15, 15)), ncol=25)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- suppressWarnings(base_kolmogorovsmirnov_twosample(x, y, alts, exact=ext))
res2 <- suppressWarnings(row_kolmogorovsmirnov_twosample(x, y, alts, exact=ext))
stopifnot(all.equal(res1, res2))

# 2 samples in one group and 101 in another (not exceeding exact cutoff)
x <- matrix(rnorm(2000), ncol=2)
y <- matrix(rnorm(101000), ncol=101)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_kolmogorovsmirnov_twosample(x, y, alts, exact=ext)
res2 <- row_kolmogorovsmirnov_twosample(x, y, alts, exact=ext)
stopifnot(all.equal(res1, res2))

# 101 samples per group (exceeding exact cutoff)
x <- matrix(rnorm(101000), ncol=101)
y <- matrix(rnorm(101000), ncol=101)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_kolmogorovsmirnov_twosample(x, y, alts, exact=ext)
res2 <- row_kolmogorovsmirnov_twosample(x, y, alts, exact=ext)
stopifnot(all.equal(res1, res2))

# 101 samples per group with ties (exceeding exact cutoff)
x <- matrix(round(runif(101000, -15, 15)), ncol=101)
y <- matrix(round(runif(101000, -15, 15)), ncol=101)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- suppressWarnings(base_kolmogorovsmirnov_twosample(x, y, alts, exact=ext))
res2 <- suppressWarnings(row_kolmogorovsmirnov_twosample(x, y, alts, exact=ext))
stopifnot(all.equal(res1, res2))

# 2 samples in one group and 5001 in another (exceeding exact cutoff)
x <- matrix(rnorm(2000), ncol=2)
y <- matrix(rnorm(5001000), ncol=5001)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
ext <- sample(c(TRUE, FALSE, NA), nrow(x), replace=TRUE)
res1 <- base_kolmogorovsmirnov_twosample(x, y, alts, exact=ext)
res2 <- row_kolmogorovsmirnov_twosample(x, y, alts, exact=ext)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
x    <- matrix(x, nrow=nrow(pars), ncol=length(x), byrow=TRUE)
y    <- c(100000000000005, 100000000000006, 100000000000009)
y    <- matrix(y, nrow=nrow(pars), ncol=length(y), byrow=TRUE)
res1 <- base_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2])
res2 <- row_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2])
stopifnot(all.equal(res1, res2))

# small numbers
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0.00000000000001)
x    <- matrix(x, nrow=nrow(pars), ncol=length(x), byrow=TRUE)
y    <- c(0.00000000000005, 0.00000000000006, 0.00000000000000)
y    <- matrix(y, nrow=nrow(pars), ncol=length(y), byrow=TRUE)
res1 <- base_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2])
res2 <- row_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2])
stopifnot(all.equal(res1, res2))

# large sample (only non-exact, because of computation time)
pars <- expand.grid(c("t","g","l"), stringsAsFactors=FALSE)
x <- rnorm(10^6)
y <- rnorm(10^6)
x <- matrix(x, nrow=nrow(pars), ncol=length(x), byrow=TRUE)
y <- matrix(y, nrow=nrow(pars), ncol=length(y), byrow=TRUE)
res1 <- base_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1])
res2 <- row_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1])
stopifnot(all.equal(res1, res2))

# x and y have infinities (exact)
x <- matrix(rnorm(100), ncol=10)
y <- matrix(rnorm(100), ncol=10)
x[1,1] <- Inf
x[2,1] <- -Inf
y[3,1] <- Inf
y[4,1] <- -Inf
x[5,1:2] <- c(Inf, -Inf)
y[6,1:2] <- c(Inf, -Inf)
x[7,1] <- Inf; y[7,1] <- Inf
x[8,1] <- Inf; y[8,1] <- -Inf
x[9,1] <- -Inf; y[9,1] <- -Inf
x[10,] <- Inf; y[10,] <- -Inf
# exact
res1 <- base_kolmogorovsmirnov_twosample(x, y)
res2 <- row_kolmogorovsmirnov_twosample(x, y)
stopifnot(all.equal(res1, res2))
# not exact
res1 <- suppressWarnings(base_kolmogorovsmirnov_twosample(x, y, exact=FALSE))
res2 <- suppressWarnings(row_kolmogorovsmirnov_twosample(x, y, exact=FALSE))
stopifnot(all.equal(res1, res2))


#--- minimal sample size -------------------------------------------------------

# single number in each group
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rnorm(nrow(pars)), ncol=1)
y    <- matrix(rnorm(nrow(pars)), ncol=1)
res1 <- base_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2])
res2 <- row_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2])
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$obs.x, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.y, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.tot, rep(2, nrow(x))))

# 3 observation in both, but only one non NA
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- cbind(matrix(rnorm(nrow(pars)), ncol=1), NA, NaN)
y    <- cbind(NaN, NA, matrix(rnorm(nrow(pars)), ncol=1))
res1 <- base_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2])
res2 <- row_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2])
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$obs.x, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.y, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.tot, rep(2, nrow(x))))

# single equal number in each group
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rep(1, nrow(pars)), ncol=1)
y    <- matrix(rep(1, nrow(pars)), ncol=1)
res1 <- suppressWarnings(base_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2]))
res2 <- suppressWarnings(row_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2]))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$obs.x, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.y, rep(1, nrow(x))))
stopifnot(all.equal(res2$obs.tot, rep(2, nrow(x))))

# three numbers in both groups, all equal
pars <- expand.grid(c("t","g","l"), c(TRUE, FALSE), stringsAsFactors=FALSE)
x    <- matrix(rnorm(nrow(pars)), ncol=1)
x    <- cbind(x,x,x)
y    <- x
res1 <- suppressWarnings(base_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2]))
res2 <- suppressWarnings(row_kolmogorovsmirnov_twosample(x, y, alternative=pars[,1], exact=pars[,2]))
stopifnot(all.equal(res1, res2))


#--- various cases for exact ---------------------------------------------------

# 99 samples with exact = NA should be the same as exact = TRUE
x <- rnorm(99)
y <- rnorm(99)
res1 <- base_kolmogorovsmirnov_twosample(x, y)
res2 <- row_kolmogorovsmirnov_twosample(x, y, exact=NA)
res3 <- row_kolmogorovsmirnov_twosample(x, y, exact=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 100 samples with one NA and exact = NA should be exact = TRUE
x <- c(rnorm(99), NA)
y <- c(NA, rnorm(99))
res1 <- base_kolmogorovsmirnov_twosample(x, y)
res2 <- row_kolmogorovsmirnov_twosample(x, y, exact=NA)
res3 <- row_kolmogorovsmirnov_twosample(x, y, exact=TRUE)
stopifnot(all.equal(res2$exact, TRUE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

# 100 samples with exact = NA should be exact = FALSE
x <- rnorm(100)
y <- rnorm(100)
res1 <- base_kolmogorovsmirnov_twosample(x, y)
res2 <- row_kolmogorovsmirnov_twosample(x, y, exact=NA)
res3 <- row_kolmogorovsmirnov_twosample(x, y, exact=FALSE)
stopifnot(all.equal(res2$exact, FALSE))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))

