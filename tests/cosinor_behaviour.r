library(matrixTests)

#--- special argument cases ----------------------------------------------------

# x can be a vector
x <- rnorm(9)
X <- matrix(x, nrow=1)
t <- 1:9
stopifnot(all.equal(row_cosinor(x, t), row_cosinor(X, t)))

# x can be numeric data.frame
x <- iris[,1:4]
X <- as.matrix(x)
t <- 1:4
stopifnot(all.equal(row_cosinor(x, t), row_cosinor(X, t)))

# x can have 0 rows
x <- matrix(0, nrow=0, ncol=4)
t <- 1:4
stopifnot(all.equal(nrow(row_cosinor(x, t)), 0))

# x can have 0 rows and 0 columns
x <- matrix(0, nrow=0, ncol=0)
t <- numeric()
stopifnot(all.equal(nrow(row_cosinor(x, t)), 0))


#--- missing values ------------------------------------------------------------

# missing values in x are removed
x <- c(NA,NA,rnorm(7),NA)
t <- 1:10
res1 <- row_cosinor(x, t)
res2 <- row_cosinor(x[!is.na(x)], t[!is.na(x)])
stopifnot(all.equal(res1, res2))

# missing values in x for whole periods are removed correctly
x <- c(rep(NA, 5), rnorm(10))
t <- rep(1:3, each=5)
res1 <- suppressWarnings(row_cosinor(x, t))
res2 <- suppressWarnings(row_cosinor(x[!is.na(x)], t[!is.na(x)]))
stopifnot(all.equal(res1, res2))

# no difference between NA and NaN in x
x1 <- c(NA,NA,1:7,NA)
x2 <- c(NaN,NaN,1:7,NaN)
t  <- 1:10
res1 <- row_cosinor(x1, t)
res2 <- row_cosinor(x2, t)
stopifnot(all.equal(res1, res2))

# missing values in t are removed
x <- rnorm(10)
t <- c(NA, NA, 3:9, NA)
res1 <- suppressWarnings(row_cosinor(x, t))
res2 <- row_cosinor(x[!is.na(t)], t[!is.na(t)])
stopifnot(all.equal(res1, res2))

# no difference between NA and NaN in t
x <- rnorm(10)
t1 <- c(NA, NA, 3:9, NA)
t2 <- c(NaN, NaN, 3:9, NaN)
res1 <- suppressWarnings(row_cosinor(x, t1))
res2 <- suppressWarnings(row_cosinor(x, t2))
stopifnot(all.equal(res1, res2))

# correct synergy between NAs in x and NAs in t
x <- c(NA, NA, rnorm(13))
t <- c(NA, 2, NA, 4:15)
res1 <- suppressWarnings(row_cosinor(x, t))
res2 <- row_cosinor(x[!is.na(x) & !is.na(t)], t[!is.na(x) & !is.na(t)])
stopifnot(all.equal(res1, res2))

# everything can be NA
x <- rep(NA_integer_, 4)
t <- rep(NA_integer_, 4)
res <- suppressWarnings(row_cosinor(x, t))
stopifnot(all.equal(res$obs, 0))


#--- rownames ------------------------------------------------------------------

# when not provided - numbers are used instead.
x <- matrix(rnorm(20), nrow=2)
t <- 1:10
res <- row_cosinor(x, t)
stopifnot(all.equal(rownames(res), c("1", "2")))

# when provided - preserved (matrix)
x <- matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B")))
t <- 1:10
res <- row_cosinor(x, t)
stopifnot(all.equal(rownames(res), rownames(x)))

# when provided - preserved (data.frame)
x <- data.frame(matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B"))))
t <- 1:10
res <- row_cosinor(x, t)
stopifnot(all.equal(rownames(res), rownames(x)))

# when duplicated - made unique
x <- matrix(rnorm(40), nrow=4, dimnames=list(c("A", "A", "B", "B")))
t <- 1:10
res <- row_cosinor(x, t)
stopifnot(all.equal(rownames(res), make.unique(rownames(x))))

