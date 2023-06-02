library(matrixTests)

#--- functions -----------------------------------------------------------------

base_kruskal <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ot <- og <- df <- chs <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])
    res <- kruskal.test(vec, factor(grp))

    # if p-value is NA turn df to NA as well
    if(is.na(res$p.value)) res$parameter <- NA

    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
    df[i]  <- res$parameter
    chs[i] <- res$statistic
    p[i]   <- res$p.value
  }

  data.frame(obs.tot=ot, obs.groups=og, df, statistic=chs, pvalue=p)
}


#--- montecarlo ----------------------------------------------------------------

# two groups
x <- matrix(rnorm(10000), ncol=10)
g <- sample(letters[1:2], 6, replace=TRUE)
g <- sample(c("a", "a", "b", "b", g))  # ensure both groups have at least 2 obs
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# lots of groups
x <- matrix(rnorm(100000), ncol=100)
g <- sample(letters[1:15], 100, replace=TRUE)
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
       100000000000003, 100000000000002, 100000000000003, 100000000000000
       )
g <- rep(c("a","b"), each=4)
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
       1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
       )
g <- rep(c("a", "b"), each=4)
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# NOTE: turned-off because of precission errors on architectures without long doubles
# large sample
# x <- rnorm(3 * 10^6)
# g <- rep(letters[1:3], each=10^6)
# res1 <- base_kruskal(x, g)
# res2 <- row_kruskalwallis(x, g)
# stopifnot(all.equal(res1, res2))

# infinities in one group
x <- c(Inf, -Inf, 3, 4, 2, 3)
g <- rep(letters[1:2], each=3)
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# infinities in both groups
x <- c(Inf, Inf, 3, 4, Inf, -Inf)
g <- rep(letters[1:2], each=3)
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# infinity is treated as highest rank
x1 <- c(Inf, Inf, 3, 4, 2, -Inf)
x2 <- c(10, 10, 3, 4, 2, -10)
g <- rep(letters[1:2], each=3)
res1 <- row_kruskalwallis(x1, g)
res2 <- row_kruskalwallis(x2, g)
stopifnot(all.equal(res1, res2))


#--- minimal sample size -------------------------------------------------------

# two groups with one value per group
x <- c(1,2)
g <- c("a","b")
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# multiple groups with one value per group
x <- rnorm(12)
g <- letters[1:12]
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# one group has one element
x <- rnorm(10)
g <- rep(letters[1:4], c(1,3,3,3))
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))


#--- ties ----------------------------------------------------------------------

# all numbers are the same
x <- rep(1, 4)
g <- c("a","a","b","b")
res1 <- base_kruskal(x, g)
res2 <- suppressWarnings(row_kruskalwallis(x, g))
stopifnot(all.equal(res1, res2))

# numbers in one group are all the same
x <- c(rep(1, 2), rnorm(4))
g <- c("a","a","b","b","c","c")
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# numbers across multiple groups are all the same
x <- c(rep(1, 4), rnorm(2))
g <- c("a","a","b","b","c","c")
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# groups are identical but have different numbers
x <- rep(rnorm(2), each=3)
g <- c("a","a","b","b","c","c")
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# one number is repeated multiple times
x <- sample(c(rep(1, 5), rnorm(10)))
g <- rep(letters[1:5], each=3)
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))


#--- constant values -----------------------------------------------------------

# all values are zero
x <- c(0,0,0,0)
g <- c("a","a","b","b")
res1 <- base_kruskal(x, g)
res2 <- suppressWarnings(row_kruskalwallis(x, g))
stopifnot(all.equal(res1, res2))

# all values are constant
x <- c(1,1,1,1)
g <- c("a","a","b","b")
res1 <- base_kruskal(x, g)
res2 <- suppressWarnings(row_kruskalwallis(x, g))
stopifnot(all.equal(res1, res2))

# within group values are constant
x <- c(1,1,2,2)
g <- c("a","a","b","b")
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# one group's values are constant
x <- c(1,1,2,3)
g <- c("a","a","b","b")
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

# groups with equal variances give equal results
x <- c(1:3,2:4,3:5,4:6,5:7)
g <- rep(1:5, each=3)
res1 <- base_kruskal(x, g)
res2 <- row_kruskalwallis(x, g)
stopifnot(all.equal(res1, res2))

