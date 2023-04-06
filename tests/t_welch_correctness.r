library(matrixTests)

#--- functions -----------------------------------------------------------------

base_t_welch <- function(mat1, mat2, null=0, alternative="two.sided", conf=0.95) {
  if(is.vector(mat1)) mat1 <- matrix(mat1, nrow=1)
  if(is.vector(mat2)) mat2 <- matrix(mat2, nrow=1)
  if(length(alternative)==1) alternative <- rep(alternative, nrow(mat1))
  if(length(null)==1) null <- rep(null, nrow(mat1))
  if(length(conf)==1) conf <- rep(conf, nrow(mat1))

  mx <- my <- md <- vx <- vy <- vp <- nx <- ny <- nt <- tst <- p <- cl <- ch <-
    se <- df <- m0 <- cnf <- numeric(nrow(mat1))
  al <- character(nrow(mat1))
  for(i in 1:nrow(mat1)) {
    vec1 <- na.omit(mat1[i,])
    vec2 <- na.omit(mat2[i,])
    res <- t.test(vec1, vec2, alternative=alternative[i], mu=null[i], conf.level=conf[i])

    # if p-value is NA turn stderr and df to NA as well
    if(is.na(res$p.value)) {
      res$stderr <- NA
      res$parameter <- NA
    }

    vx[i]  <- var(vec1)
    vy[i]  <- var(vec2)
    nx[i]  <- length(vec1)
    ny[i]  <- length(vec2)
    nt[i]  <- nx[i] + ny[i]
    mx[i]  <- res$estimate[1]
    my[i]  <- res$estimate[2]
    md[i]  <- mx[i]-my[i]
    tst[i] <- res$statistic
    p[i]   <- res$p.value
    cl[i]  <- res$conf.int[1]
    ch[i]  <- res$conf.int[2]
    df[i]  <- res$parameter
    m0[i]  <- res$null.value
    al[i]  <- res$alternative
    se[i]  <- res$stderr
    cnf[i] <- attr(res$conf.int, "conf.level")
  }

  data.frame(obs.x=nx, obs.y=ny, obs.tot=nt, mean.x=mx, mean.y=my, mean.diff=md,
             var.x=vx, var.y=vy, stderr=se, df=df, statistic=tst, pvalue=p,
             conf.low=cl, conf.high=ch, mean.null=m0, alternative=al,
             conf.level=cnf, stringsAsFactors=FALSE
             )
}


#--- montecarlo ----------------------------------------------------------------

# 3 and 2 observations
x <- matrix(rnorm(3000), ncol=3)
y <- matrix(rnorm(2000), ncol=2)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- sample(seq(-1, 1, length.out=nrow(x)))
cfs  <- sample(seq(0, 1, length.out=nrow(x)))
res1 <- base_t_welch(x, y, mus, alts, cfs)
res2 <- row_t_welch(x, y, mus, alts, cfs)
stopifnot(all.equal(res1, res2))

# 20 observations in each group
x <- matrix(rnorm(20000), ncol=20)
y <- matrix(rnorm(20000), ncol=20)
alts <- sample(c("t", "g", "l"), nrow(x), replace=TRUE)
mus  <- sample(seq(-1, 1, length.out=nrow(x)))
cfs  <- sample(seq(0, 1, length.out=nrow(x)))
res1 <- base_t_welch(x, y, mus, alts, cfs)
res2 <- row_t_welch(x, y, mus, alts, cfs)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000)
y <- c(100000000000003, 100000000000002, 100000000000003)
res1 <- base_t_welch(x, y)
res2 <- row_t_welch(x, y)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(0.00000000000004, 0.00000000000002, 0.00000000000003, 0)
y <- c(0.00000000000003, 0.00000000000002, 0.00000000000003)
res1 <- base_t_welch(x, y)
res2 <- row_t_welch(x, y)
stopifnot(all.equal(res1, res2))

# large sample
x <- rnorm(10^6)
y <- rnorm(10^6)
res1 <- base_t_welch(x, y)
res2 <- row_t_welch(x, y)
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# 2 observations in both groups
x <- matrix(rnorm(6), ncol=2)
y <- matrix(rnorm(6), ncol=2)
alt <- c("two.sided", "greater", "less")
res1 <- base_t_welch(x, y, alternative=alt)
res2 <- row_t_welch(x, y, alternative=alt)
stopifnot(all.equal(res1, res2))

# 2 observations in first group and 4 in second
x <- matrix(rnorm(6), ncol=2)
y <- matrix(rnorm(12), ncol=4)
alt <- c("two.sided", "greater", "less")
res1 <- base_t_welch(x, y, alternative=alt)
res2 <- row_t_welch(x, y, alternative=alt)
stopifnot(all.equal(res1, res2))

# 3 observations in both groups but one is NA
x <- matrix(c(rnorm(6), NA,NA,NA), ncol=3)
y <- matrix(c(NA,NA,NA,rnorm(6)), ncol=3)
alt <- c("two.sided", "greater", "less")
res1 <- base_t_welch(x, y, alternative=alt)
res2 <- row_t_welch(x, y, alternative=alt)
stopifnot(all.equal(res1, res2))


#--- parameter edge cases ------------------------------------------------------

# various corner cases with NAs
mus <- c(-Inf, -1, 0, 1, Inf)
alt <- c("l", "t", "g")
cfs <- c(0, 0.5, 1)
pars <- expand.grid(mus, alt, cfs, stringsAsFactors=FALSE)
x <- matrix(rnorm(10*nrow(pars)), ncol=10)
y <- matrix(rnorm(10*nrow(pars)), ncol=10)
res1 <- base_t_welch(x, y, pars[,1], pars[,2], pars[,3])
res2 <- row_t_welch(x, y, pars[,1], pars[,2], pars[,3])
stopifnot(all.equal(res1, res2))

# NAs in confidence intervals
x <- matrix(rnorm(40), ncol=10)
y <- matrix(rnorm(40), ncol=10)
cnf <- c(0.95, NA, 0.5, NA)
res1 <- base_t_welch(x, y, conf=ifelse(is.na(cnf), 0.95, cnf))
res1[is.na(cnf), c("conf.level", "conf.low", "conf.high")] <- NA
res2 <- row_t_welch(x, y, conf.level=cnf)
stopifnot(all.equal(res1, res2))

# null exactly equal to the mean
res1 <- base_t_welch(c(1,2,3), c(0,0,0), null=2)
res2 <- row_t_welch(c(1,2,3), c(0,0,0), null=2)
stopifnot(all.equal(res2$pvalue, 1))
stopifnot(all.equal(res1, res2))


#--- constant values -----------------------------------------------------------

# first group values are constant
x <- matrix(1, nrow=3, ncol=3)
y <- matrix(rnorm(9), nrow=3, ncol=3)
alt <- c("l", "t", "g")
res1 <- base_t_welch(x, y, alternative=alt)
res2 <- row_t_welch(x, y, alternative=alt)
stopifnot(all.equal(res2$var.x, rep(0, nrow(x))))

# second group values are constant
x <- matrix(rnorm(9), nrow=3, ncol=3)
y <- matrix(1, nrow=3, ncol=3)
alt <- c("l", "t", "g")
res1 <- base_t_welch(x, y, alternative=alt)
res2 <- row_t_welch(x, y, alternative=alt)
stopifnot(all.equal(res2$var.y, rep(0, nrow(y))))

