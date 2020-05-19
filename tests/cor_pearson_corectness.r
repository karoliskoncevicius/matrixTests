library(matrixTests)

#--- functions -----------------------------------------------------------------

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
    cnf[i] <- ifelse(is.null(res$conf.int), conf[i], attr(res$conf.int, "conf.level"))
    df[i]  <- res$parameter
    mu[i]  <- res$null.value
    al[i]  <- res$alternative
  }

  data.frame(obs.paired=np, cor=cor, df=df, statistic=tst,
             pvalue=p, conf.low=cl, conf.high=ch, alternative=al,
             cor.null=mu, conf.level=cnf, stringsAsFactors=FALSE
             )
}


#--- montecarlo ----------------------------------------------------------------

# 5 observations
x <- matrix(rnorm(5000), ncol=5)
y <- matrix(rnorm(5000), ncol=5)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
cfs  <- sample(seq(0, 1, length.out=nrow(x)))
res1 <- base_cor_pearson(x, y, alts, cfs)
res2 <- row_cor_pearson(x, y, alts, cfs)
stopifnot(all.equal(res1, res2))

# 20 observations
x <- matrix(rnorm(20000), ncol=20)
y <- matrix(rnorm(20000), ncol=20)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
cfs  <- sample(seq(0, 1, length.out=nrow(x)))
res1 <- base_cor_pearson(x, y, alts, cfs)
res2 <- row_cor_pearson(x, y, alts, cfs)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
y <- c(100000000000003, 100000000000002, 100000000000003, 100000000000000)
res1 <- base_cor_pearson(x, y)
res2 <- row_cor_pearson(x, y)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1)
y <- c(1.00000000000003, 1.00000000000002, 1.00000000000003, 1)
res1 <- base_cor_pearson(x, y)
res2 <- row_cor_pearson(x, y)
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# 3 values
x <- matrix(rnorm(9), ncol=3)
y <- matrix(rnorm(9), ncol=3)
alt <- c("two.sided", "greater", "less")
res1 <- base_cor_pearson(x, y)
res2 <- suppressWarnings(row_cor_pearson(x, y))
stopifnot(all.equal(res1, res2))

# three numbers with NAs
x <- matrix(c(rnorm(12), NA, NA, NA), ncol=5)
y <- matrix(c(NA, NA, NA, rnorm(12)), ncol=5)
alt <- c("two.sided", "greater", "less")
res1 <- base_cor_pearson(x, y)
res2 <- suppressWarnings(row_cor_pearson(x, y))
stopifnot(all.equal(res1, res2))

# four numbers (will produce confidence intervals)
x <- matrix(rnorm(12), ncol=4)
y <- matrix(rnorm(12), ncol=4)
alt <- c("two.sided", "greater", "less")
res1 <- base_cor_pearson(x, y)
res2 <- row_cor_pearson(x, y)
stopifnot(all.equal(res1, res2))


#--- parameter edge cases ------------------------------------------------------

# various corner cases with NAs
alt <- c("l", "t", "g")
cfs <- c(0, 0.5, 1)
pars <- expand.grid(alt, cfs, stringsAsFactors=FALSE)
x1 <- matrix(rnorm(10*nrow(pars)), ncol=10)
x2 <- matrix(rnorm(10*nrow(pars)), ncol=10)
x1[sample(length(x1), nrow(pars)*2)] <- NA
x2[sample(length(x2), nrow(pars)*2)] <- NA
res1 <- base_cor_pearson(x1, x2, pars[,1], pars[,2])
res2 <- suppressWarnings(row_cor_pearson(x1, x2, pars[,1], pars[,2]))
stopifnot(all.equal(res1, res2))


#--- constant values -----------------------------------------------------------

# all values are 0
res1 <- suppressWarnings(base_cor_pearson(rep(0,4), rep(0,4)))
res2 <- suppressWarnings(row_cor_pearson(rep(0,4), rep(0,4)))
stopifnot(all.equal(res1, res2))

# all values are constant
res1 <- suppressWarnings(base_cor_pearson(rep(1,4), rep(2,4)))
res2 <- suppressWarnings(row_cor_pearson(rep(1,4), rep(2,4)))
stopifnot(all.equal(res1, res2))

# perfect positive correlation
res1 <- base_cor_pearson(1:4, 1:4)
res2 <- suppressWarnings(row_cor_pearson(1:4, 1:4))
stopifnot(all.equal(res1, res2))

# perfect negative correlation
res1 <- base_cor_pearson(1:4, 4:1)
res2 <- suppressWarnings(row_cor_pearson(1:4, 4:1))
stopifnot(all.equal(res1, res2))

