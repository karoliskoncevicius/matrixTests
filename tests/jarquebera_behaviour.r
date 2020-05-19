library(matrixTests)

#--- special argument cases ----------------------------------------------------

# x can be a vector
x <- rnorm(9)
X <- matrix(x, nrow=1)
g <- rep(letters[1:3], each=3)
stopifnot(all.equal(row_jarquebera(x), row_jarquebera(X)))

# x can be numeric data.frame
x <- iris[,1:4]
X <- as.matrix(x)
g <- c("a","a","b","b")
stopifnot(all.equal(row_jarquebera(x), row_jarquebera(X)))

# x can have 0 rows
x <- matrix(0, nrow=0, ncol=4)
g <- c("a","a","b","b")
stopifnot(all.equal(nrow(row_jarquebera(x)), 0))

# x can have 0 rows and 0 columns
x <- matrix(0, nrow=0, ncol=0)
g <- character()
stopifnot(all.equal(nrow(row_jarquebera(x)), 0))


#--- missing values ------------------------------------------------------------

# missing values in x are removed
x <- c(NA,NA,rnorm(7),NA)
res1 <- row_jarquebera(x)
res2 <- row_jarquebera(x[!is.na(x)])
stopifnot(all.equal(res1, res2))

# no difference between NA and NaN in x
x1 <- c(NA,NA,1:7,NA)
x2 <- c(NaN,NaN,1:7,NaN)
res1 <- row_jarquebera(x1)
res2 <- row_jarquebera(x2)
stopifnot(all.equal(res1, res2))

# everything can be NA
x <- rep(NA_integer_, 4)
res <- suppressWarnings(row_jarquebera(x))
stopifnot(all.equal(res$obs, 0))


#--- rownames ------------------------------------------------------------------

# when not provided - numbers are used instead.
x <- matrix(rnorm(20), nrow=2)
res <- row_jarquebera(x)
stopifnot(all.equal(rownames(res), c("1", "2")))

# when provided - preserved (matrix)
x <- matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B")))
res <- row_jarquebera(x)
stopifnot(all.equal(rownames(res), rownames(x)))

# when provided - preserved (data.frame)
x <- data.frame(matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B"))))
res <- row_jarquebera(x)
stopifnot(all.equal(rownames(res), rownames(x)))

# when duplicated - made unique
x <- matrix(rnorm(40), nrow=4, dimnames=list(c("A", "A", "B", "B")))
res <- row_jarquebera(x)
stopifnot(all.equal(rownames(res), make.unique(rownames(x))))

