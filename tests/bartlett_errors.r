library(matrixTests)
source("utils/capture.r")

#--- x argument errors ---------------------------------------------------------

# cannot be missing
err <- 'argument "x" is missing, with no default'
res <- capture(row_bartlett())
stopifnot(all.equal(res$error, err))

# cannot be NULL
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_bartlett(NULL, "a"))
stopifnot(all.equal(res$error, err))

# cannot be character
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_bartlett(c("a", "b"), c("a","b")))
stopifnot(all.equal(res$error, err))

# cannot be logical
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_bartlett(c(TRUE, FALSE), "a"))
stopifnot(all.equal(res$error, err))

# cannot be complex
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_bartlett(complex(c(1,2), c(3,4)), "a"))
stopifnot(all.equal(res$error, err))

# cannot be data.frame containing some non numeric data
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_bartlett(iris, "a"))
stopifnot(all.equal(res$error, err))

# cannot be a list
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_bartlett(as.list(c(1:5)), "a"))
stopifnot(all.equal(res$error, err))

# cannot be in a list
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_bartlett(list(1:5), "a"))
stopifnot(all.equal(res$error, err))


#--- g argument errors ---------------------------------------------------------

# cannot be missing
err <- 'argument "g" is missing, with no default'
res <- capture(row_bartlett(1, ))
stopifnot(all.equal(res$error, err))

# cannot be a matrix with dimensions
err <- '"g" must be a vector with length ncol(x)'
res <- capture(row_bartlett(matrix(1:4, nrow=1), matrix(letters[1:4], nrow=2)))
stopifnot(all.equal(res$error, err))

# cannot be a in a list
err <- '"g" must be a vector with length ncol(x)'
res <- capture(row_bartlett(1:5, list(1:5)))
stopifnot(all.equal(res$error, err))


#--- dimension mismatch errors -------------------------------------------------

# g length must match the observations of x
err <- '"g" must be a vector with length ncol(x)'
res <- capture(row_bartlett(1:5, letters[1:4]))
stopifnot(all.equal(res$error, err))

