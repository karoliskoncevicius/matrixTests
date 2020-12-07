library(matrixTests)

#--- functions -----------------------------------------------------------------

base_bartlett <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ng <- nt <- vt <- ks <- p <- df <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])

    badgr <- names(which(table(grp)==1))
    bad   <- grp %in% badgr
    vec   <- vec[!bad]
    grp   <- grp[!bad]

    res <- bartlett.test(vec, grp)

    # if p-value is NA turn df to NA as well
    if(is.na(res$p.value)) res$parameter <- NA

    ng[i] <- length(unique(grp))
    nt[i] <- length(vec)
    vt[i] <- sum((tapply(vec, grp, length)-1) * tapply(vec, grp, var, na.rm=TRUE)) / (nt[i]-ng[i])
    ks[i] <- res$statistic
    p[i]  <- res$p.value
    df[i]  <- res$parameter
  }

  data.frame(obs.tot=nt, obs.groups=ng, var.pooled=vt, df=df, statistic=ks,
             pvalue=p, stringsAsFactors=FALSE
             )
}


#--- montecarlo ----------------------------------------------------------------

# two groups
x <- matrix(rnorm(10000), ncol=10)
g <- sample(letters[1:2], 10, replace=TRUE)
res1 <- base_bartlett(x, factor(g))
res2 <- row_bartlett(x, g)
stopifnot(all.equal(res1, res2))

# lots of groups
x <- matrix(rnorm(100000), ncol=100)
g <- sample(letters[1:15], 100, replace=TRUE)
res1 <- base_bartlett(x, factor(g))
res2 <- row_bartlett(x, g)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
       100000000000003, 100000000000002, 100000000000003, 100000000000000
       )
g <- c(rep("a", 4), rep("b", 4))
res1 <- base_bartlett(x, g)
res2 <- row_bartlett(x, g)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
       1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
       )
g <- c(rep("a", 4), rep("b", 4))
res1 <- base_bartlett(x, g)
res2 <- row_bartlett(x, g)
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# two groups with two values each
x <- rnorm(4)
g <- c("a","a","b","b")
res1 <- base_bartlett(x, g)
res2 <- row_bartlett(x, g)
stopifnot(all.equal(res1, res2))

# many groups all with two values
x <- rnorm(20)
g <- rep(letters[1:10], each=2)
res1 <- base_bartlett(x, g)
res2 <- row_bartlett(x, g)
stopifnot(all.equal(res1, res2))


#--- constant values -----------------------------------------------------------

# all values are zero
x <- c(0,0,0,0)
g <- c("a","a","b","b")
res1 <- base_bartlett(x, g)
res2 <- suppressWarnings(row_bartlett(x, g))
stopifnot(all.equal(res1, res2))

# all values are constant
x <- c(1,1,1,1)
g <- c("a","a","b","b")
res1 <- base_bartlett(x, g)
res2 <- suppressWarnings(row_bartlett(x, g))
stopifnot(all.equal(res1, res2))

# within group values are constant
x <- c(1,1,2,2)
g <- c("a","a","b","b")
res1 <- base_bartlett(x, g)
res2 <- suppressWarnings(row_bartlett(x, g))
stopifnot(all.equal(res1, res2))

# one group's values are constant
x <- c(1,1,2,3)
g <- c("a","a","b","b")
res1 <- base_bartlett(x, g)
res2 <- suppressWarnings(row_bartlett(x, g))
stopifnot(all.equal(res1, res2))

# all groups have equal variances
x <- c(1:3,2:4,3:5,4:6,5:7)
g <- rep(1:5, each=3)
res1 <- base_bartlett(x, g)
res2 <- row_bartlett(x, g)
stopifnot(all.equal(res1, res2))


#--- groups with one element ---------------------------------------------------

# groups with one remaining values are removed
x <- rnorm(5)
g <- c("a","a","b","b","c")
res1 <- base_bartlett(x, g)
res2 <- suppressWarnings(row_bartlett(x, g))
res3 <- suppressWarnings(row_bartlett(x[-5], g[-5]))
stopifnot(all.equal(res1, res2))
stopifnot(all.equal(res2, res3))
stopifnot(all.equal(res2$obs.groups, 2))

