library(matrixTests)

#--- special argument cases ----------------------------------------------------

# x can be a vector
x <- rnorm(9)
X <- matrix(x, nrow=1)
g <- rep(letters[1:3], each=3)
stopifnot(all.equal(row_brownforsythe(x, g), row_brownforsythe(X, g)))

# x can be numeric data.frame
x <- iris[,1:4]
X <- as.matrix(x)
g <- c("a","a","a","b")
stopifnot(all.equal(row_brownforsythe(x, g), row_brownforsythe(X, g)))

# x can have 0 rows
x <- matrix(0, nrow=0, ncol=4)
g <- c("a","a","b","b")
stopifnot(all.equal(nrow(row_brownforsythe(x, g)), 0))

# x can have 0 rows and 0 columns
x <- matrix(0, nrow=0, ncol=0)
g <- character()
stopifnot(all.equal(nrow(row_brownforsythe(x, g)), 0))

# groups can be infinities
x <- rnorm(5)
g1 <- c(Inf, Inf, Inf, -Inf, -Inf)
g2 <- c("a","a","a","b","b")
stopifnot(all.equal(row_brownforsythe(x, g1), row_brownforsythe(x, g2)))


#--- missing values ------------------------------------------------------------

# missing values in x are removed
x <- c(NA,NA,rnorm(7),NA)
g <- rep(letters[1:2], each=5)
res1 <- row_brownforsythe(x, g)
res2 <- row_brownforsythe(x[!is.na(x)], g[!is.na(x)])
stopifnot(all.equal(res1, res2))

# missing values in x for whole group are removed correctly
x <- c(rep(NA, 5), rnorm(10))
g <- rep(letters[1:3], each=5)
res1 <- row_brownforsythe(x, g)
res2 <- row_brownforsythe(x[!is.na(x)], g[!is.na(x)])
stopifnot(all.equal(res1, res2))

# no difference between NA and NaN in x
x1 <- c(NA,NA,1:7,NA)
x2 <- c(NaN,NaN,1:7,NaN)
g <- rep(letters[1:2], each=5)
res1 <- row_brownforsythe(x1, g)
res2 <- row_brownforsythe(x2, g)
stopifnot(all.equal(res1, res2))

# missing values in g are removed
x <- rnorm(10)
g <- c(NA, NA, rep(letters[1:2], c(3,4)), NA)
res1 <- suppressWarnings(row_brownforsythe(x, g))
res2 <- row_brownforsythe(x[!is.na(g)], g[!is.na(g)])
stopifnot(all.equal(res1, res2))

# no difference between NA and NaN in g
x <- rnorm(10)
g1 <- c(NA, NA, rep(1:2, c(3,4)), NA)
g2 <- c(NaN, NaN, rep(1:2, c(3,4)), NaN)
res1 <- suppressWarnings(row_brownforsythe(x, g1))
res2 <- suppressWarnings(row_brownforsythe(x, g2))
stopifnot(all.equal(res1, res2))

# correct synergy between NAs in x and NAs in g
x <- c(NA, NA, rnorm(13))
g <- c(NA, "a", NA, rep(letters[2:5], each=3))
res1 <- suppressWarnings(row_brownforsythe(x, g))
res2 <- row_brownforsythe(x[!is.na(x) & !is.na(g)], g[!is.na(x) & !is.na(g)])
stopifnot(all.equal(res1, res2))

# everything can be NA
x <- rep(NA_integer_, 4)
g <- rep(NA, 4)
res <- suppressWarnings(row_brownforsythe(x, g))
stopifnot(all.equal(res$obs.tot, 0))
stopifnot(all.equal(res$obs.groups, 0))


#--- rownames ------------------------------------------------------------------

# when not provided - numbers are used instead.
x <- matrix(rnorm(20), nrow=2)
g <- rep(letters[1:2], each=5)
res <- row_brownforsythe(x, g)
stopifnot(all.equal(rownames(res), c("1", "2")))

# when provided - preserved (matrix)
x <- matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B")))
g <- rep(letters[1:2], each=5)
res <- row_brownforsythe(x, g)
stopifnot(all.equal(rownames(res), rownames(x)))

# when provided - preserved (data.frame)
x <- data.frame(matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B"))))
g <- rep(letters[1:2], each=5)
res <- row_brownforsythe(x, g)
stopifnot(all.equal(rownames(res), rownames(x)))

# when duplicated - made unique
x <- matrix(rnorm(40), nrow=4, dimnames=list(c("A", "A", "B", "B")))
g <- rep(letters[1:2], each=5)
res <- row_brownforsythe(x, g)
stopifnot(all.equal(rownames(res), make.unique(rownames(x))))

