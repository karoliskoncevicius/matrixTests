library(matrixTests)

#--- functions -----------------------------------------------------------------

# NOTE: we are using two separate versions of base: aov() and oneway.test()
#       aov() has some roudning errors so second version is used for sensitivity
#       oneway.test() does not allow having one observation per group

base_oneway_equalvar <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  st <- sr <- mt <- mr <- ot <- og <- dft <- dfr <- fst <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- mat[i,!bad]
    grp <- factor(groups[!bad])

    res <- summary(aov(vec ~ grp))[[1]]

    # if p-value is NA then turn dfs to NA as well
    if(is.na(res[1,5])) res[,1] <- NA

    st[i]  <- res[1,2]
    sr[i]  <- res[2,2]
    mt[i]  <- res[1,3]
    mr[i]  <- res[2,3]
    dft[i] <- res[1,1]
    dfr[i] <- res[2,1]
    fst[i] <- res[1,4]
    p[i]   <- res[1,5]
    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
  }

  data.frame(obs.tot=ot, obs.groups=og, sumsq.between=st, sumsq.within=sr,
             meansq.between=mt, meansq.within=mr, df.between=dft,
             df.within=dfr, statistic=fst, pvalue=p
             )
}

base_oneway_equalvar2 <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  st <- sr <- mt <- mr <- ot <- og <- dft <- dfr <- fst <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- mat[i,!bad]
    grp <- factor(groups[!bad])

    res <- oneway.test(vec ~ grp, var.equal=TRUE)

    # if p-value is NA then turn dfs to NA as well
    if(is.na(res$p.value)) res$parameter <- NA

    dft[i] <- res$parameter[1]
    dfr[i] <- res$parameter[2]
    fst[i] <- res$statistic
    p[i]   <- res$p.value
    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
    st[i]  <- sum(tapply(vec, grp, length) * (tapply(vec, grp, mean)-mean(vec))^2)
    sr[i]  <- sum((tapply(vec, grp, length)-1) * tapply(vec, grp, var))
    mt[i]  <- st[i] / dft[i]
    mr[i]  <- sr[i] / dfr[i]
  }

  data.frame(obs.tot=ot, obs.groups=og, sumsq.between=st, sumsq.within=sr,
             meansq.between=mt, meansq.within=mr, df.between=dft,
             df.within=dfr, statistic=fst, pvalue=p
             )
}


#--- montecarlo ----------------------------------------------------------------

# two groups
x <- matrix(rnorm(10000), ncol=10)
g <- sample(letters[1:2], 6, replace=TRUE)
g <- sample(c("a", "a", "b", "b", g))  # ensure both groups have at least 2 obs
res1 <- base_oneway_equalvar(x, factor(g))
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))

# lots of groups
x <- matrix(rnorm(100000), ncol=100)
g <- sample(letters[1:15], 100, replace=TRUE)
res1 <- base_oneway_equalvar(x, factor(g))
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
       100000000000003, 100000000000002, 100000000000003, 100000000000000
       )
g <- c(rep("a", 4), rep("b", 4))
res1 <- base_oneway_equalvar2(x, g)
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
       1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
       )
g <- c(rep("a", 4), rep("b", 4))
res1 <- base_oneway_equalvar2(x, g)
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))

# large sample
x <- rnorm(3 * 10^6)
g <- rep(letters[1:3], each=10^6)
res1 <- base_oneway_equalvar2(x, g)
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# two groups one with two values and another with one value
x <- c(1,2,3)
g <- c("a","a","b")
res1 <- base_oneway_equalvar(x, g)
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))

# two groups one with three values and another with one value
x <- c(1,2,3,4)
g <- c("a","a","a","b")
res1 <- base_oneway_equalvar(x, g)
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))

# two groups both with two values
x <- c(1,2,3,4)
g <- c("a","a","b","b")
res1 <- base_oneway_equalvar(x, g)
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))

# multiple groups with one value and one group with two values
x <- rnorm(12)
g <- c("a", letters[1:11])
res1 <- base_oneway_equalvar(x, g)
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))


#--- constant values -----------------------------------------------------------

# all values are zero
x <- c(0,0,0,0)
g <- c("a","a","b","b")
res1 <- base_oneway_equalvar(x, g)
res2 <- suppressWarnings(row_oneway_equalvar(x, g))
stopifnot(all.equal(res1, res2))

# all values are constant
x <- c(1,1,1,1)
g <- c("a","a","b","b")
res1 <- base_oneway_equalvar(x, g)
res2 <- suppressWarnings(row_oneway_equalvar(x, g))
stopifnot(all.equal(res1, res2))

# within group values are constant
x <- c(1,1,2,2)
g <- c("a","a","b","b")
res1 <- base_oneway_equalvar2(x, g)
res2 <- suppressWarnings(row_oneway_equalvar(x, g))
stopifnot(all.equal(res1, res2))

# one group's values are constant
x <- c(1,1,2,3)
g <- c("a","a","b","b")
res1 <- base_oneway_equalvar(x, g)
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))

# all groups have equal variances
x <- c(1:3,2:4,3:5,4:6,5:7)
g <- rep(1:5, each=3)
res1 <- base_oneway_equalvar(x, g)
res2 <- row_oneway_equalvar(x, g)
stopifnot(all.equal(res1, res2))

