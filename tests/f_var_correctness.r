library(matrixTests)

#--- functions -----------------------------------------------------------------

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

    # if p-value is NA turn dfs to NA as well
    if(is.na(res$p.value)) res$parameter <- NA

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


#--- montecarlo ----------------------------------------------------------------

# 5 observations
x <- matrix(rnorm(5000), ncol=5)
y <- matrix(rnorm(5000), ncol=5)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
rats <- sample(seq(0.001, 1000, length.out=nrow(x)))
cfs  <- sample(seq(0.001, 0.999, length.out=nrow(x)))  # NOTE: using 0.001 instead of 0 because base doesn't allow 0
res1 <- base_f_var(x, y, rats, alts, cfs)
res2 <- row_f_var(x, y, rats, alts, cfs)
stopifnot(all.equal(res1, res2))

# 20 observations
x <- matrix(rnorm(20000), ncol=20)
y <- matrix(rnorm(20000), ncol=20)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
rats <- sample(seq(0.001, 1000, length.out=nrow(x)))
cfs  <- sample(seq(0.001, 0.999, length.out=nrow(x)))
res1 <- base_f_var(x, y, rats, alts, cfs)
res2 <- row_f_var(x, y, rats, alts, cfs)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
y <- c(100000000000003, 100000000000002, 100000000000003, 100000000000000)
res1 <- base_f_var(x, y)
res2 <- row_f_var(x, y)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0)
y <- c(0.00000000000003, 0.00000000000002, 0.00000000000003, 0)
res1 <- base_f_var(x, y)
res2 <- row_f_var(x, y)
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# two values in both groups
x <- matrix(rnorm(6), ncol=2)
y <- matrix(rnorm(6), ncol=2)
a <- c("two.sided", "greater", "less")
res1 <- base_f_var(x, y, alt=a)
res2 <- row_f_var(x, y, alt=a)
stopifnot(all.equal(res1, res2))

# two values in both groups with NAs
x <- matrix(c(rnorm(6), NA, NA, NA), ncol=3)
y <- matrix(rnorm(6), ncol=2)
a <- c("two.sided", "greater", "less")
res1 <- base_f_var(x, y, alt=a)
res2 <- row_f_var(x, y, alt=a)
stopifnot(all.equal(res1, res2))


#--- parameter edge cases ------------------------------------------------------

# various corner cases with NAs
alt <- c("l", "t", "g")
rat <- c(0.000001, 0.5, 1, 2, 99999)
cfs <- c(0.000001, 0.5, 0.9999)
pars <- expand.grid(rat, alt, cfs, stringsAsFactors=FALSE)
x <- matrix(rnorm(10*nrow(pars)), ncol=10)
y <- matrix(rnorm(10*nrow(pars)), ncol=10)
res1 <- base_f_var(x, y, pars[,1], pars[,2], pars[,3])
res2 <- row_f_var(x, y, pars[,1], pars[,2], pars[,3])
stopifnot(all.equal(res1, res2))


#--- constant values -----------------------------------------------------------

# minimal sample size and both groups have the same values
x <- matrix(0, nrow=3, ncol=2)
y <- matrix(0, nrow=3, ncol=2)
alt <- c("l", "t", "g")
res1 <- base_f_var(x, y, alt=alt)
res1$conf.low[1]  <- NA  # TODO: check why this case returns conf.low = 0 for first row
res1$conf.high[3] <- NA  # TODO: check why this case returns conf.high = Inf for last row
res2 <- suppressWarnings(row_f_var(x, y, alt=alt))
stopifnot(all.equal(res1, res2))

# minimal sample size and both groups have zero variance
x <- matrix(1, nrow=3, ncol=2)
y <- matrix(2, nrow=3, ncol=2)
alt <- c("l", "t", "g")
res1 <- base_f_var(x, y, alt=alt)
res1$conf.low[1]  <- NA  # TODO: check why this case returns conf.low = 0 for first row
res1$conf.high[3] <- NA  # TODO: check why this case returns conf.high = Inf for last row
res2 <- suppressWarnings(row_f_var(x, y, alt=alt))
stopifnot(all.equal(res1, res2))

# variance is equal in both groups
x <- matrix(c(1,1,1,2,2,2), nrow=3, ncol=2)
y <- matrix(c(1,1,1,2,2,2), nrow=3, ncol=2)
alt <- c("l", "t", "g")
res1 <- base_f_var(x, y, alt=alt)
res2 <- suppressWarnings(row_f_var(x, y, alt=alt))
stopifnot(all.equal(res1, res2))

# parameter edge cases when x is constant
alt <- c("l", "t", "g")
rat <- c(0.000001, 0.5, 1, 2, 99999)
cfs <- c(0.000001, 0.5, 0.9999)
pars <- expand.grid(rat, alt, cfs, stringsAsFactors=FALSE)
x <- matrix(1, nrow=nrow(pars), ncol=10)
y <- matrix(rnorm(10*nrow(pars)), ncol=10)
res1 <- base_f_var(x, y, pars[,1], pars[,2], pars[,3])
res2 <- suppressWarnings(row_f_var(x, y, pars[,1], pars[,2], pars[,3]))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$var.x, rep(0, nrow(x))))

# parameter edge cases when y is constant
alt <- c("l", "t", "g")
rat <- c(0.000001, 0.5, 1, 2, 99999)
cfs <- c(0.000001, 0.5, 0.9999)
pars <- expand.grid(rat, alt, cfs, stringsAsFactors=FALSE)
x <- matrix(rnorm(10*nrow(pars)), ncol=10)
y <- matrix(1, nrow=nrow(pars), ncol=10)
res1 <- base_f_var(x, y, pars[,1], pars[,2], pars[,3])
res2 <- suppressWarnings(row_f_var(x, y, pars[,1], pars[,2], pars[,3]))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2$var.y, rep(0, nrow(y))))

