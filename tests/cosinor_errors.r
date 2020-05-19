library(matrixTests)
source("utils/capture.r")

#--- x argument errors ---------------------------------------------------------

# cannot be missing
err <- 'argument "x" is missing, with no default'
res <- capture(row_cosinor())
stopifnot(all.equal(res$error, err))

# cannot be NULL
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_cosinor(NULL, "a"))
stopifnot(all.equal(res$error, err))

# cannot be character
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_cosinor(c("a", "b"), c("a","b")))
stopifnot(all.equal(res$error, err))

# cannot be logical
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_cosinor(c(TRUE, FALSE), "a"))
stopifnot(all.equal(res$error, err))

# cannot be complex
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_cosinor(complex(c(1,2), c(3,4)), "a"))
stopifnot(all.equal(res$error, err))

# cannot be data.frame containing some non numeric data
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_cosinor(iris, "a"))
stopifnot(all.equal(res$error, err))

# cannot be a list
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_cosinor(as.list(c(1:5)), "a"))
stopifnot(all.equal(res$error, err))

# cannot be in a list
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_cosinor(list(1:5), "a"))
stopifnot(all.equal(res$error, err))


#--- t argument errors ---------------------------------------------------------

# TODO: decide what error to produce when t is infinity

# cannot be missing
err <- 'argument "t" is missing, with no default'
res <- capture(row_cosinor(1))
stopifnot(all.equal(res$error, err))

# cannot be character
err <- '"t" must be a numeric vector with length ncol(x)'
res <- capture(row_cosinor(matrix(1:4, nrow=1), c("1","2","3","4")))
stopifnot(all.equal(res$error, err))

# cannot be a matrix with dimensions
err <- '"t" must be a numeric vector with length ncol(x)'
res <- capture(row_cosinor(matrix(1:4, nrow=1), matrix(letters[1:4], nrow=2)))
stopifnot(all.equal(res$error, err))

# cannot be a in a list
err <- '"g" must be a numeric vector with length ncol(x)'
res <- capture(row_cosinor(1:5, list(1:5)))
stopifnot(all.equal(res$error, err))


#--- period argument errors ----------------------------------------------------

# cannot be character
err <- '"period" must be a numeric vector with length 1'
res <- capture(row_cosinor(1:5, 1:5, "24"))
stopifnot(all.equal(res$error, err))

# cannot be complex
err <- '"period" must be a numeric vector with length 1'
res <- capture(row_cosinor(1:5, 1:5, complex(24)))
stopifnot(all.equal(res$error, err))

# cannot be in a list
err <- '"period" must be a numeric vector with length 1'
res <- capture(row_cosinor(1:5, 1:5, list(24)))
stopifnot(all.equal(res$error, err))

# cannot be in a data frame
err <- '"period" must be a numeric vector with length 1'
res <- capture(row_cosinor(1:5, 1:5, data.frame(24)))
stopifnot(all.equal(res$error, err))

# cannot be NA
err <- 'all "period" values must be greater than 0 and lower than Inf'
res <- capture(row_cosinor(1:5, 1:5, NA_integer_))
stopifnot(all.equal(res$error, err))

# cannot be NaN
err <- 'all "period" values must be greater than 0 and lower than Inf'
res <- capture(row_cosinor(1:5, 1:5, NaN))
stopifnot(all.equal(res$error, err))

# must be above 0
err <- 'all "period" values must be greater than 0 and lower than Inf'
res <- capture(row_cosinor(1:5, 1:5, 0))
stopifnot(all.equal(res$error, err))

# must be below Infinity
err <- 'all "period" values must be greater than 0 and lower than Inf'
res <- capture(row_cosinor(1:5, 1:5, Inf))
stopifnot(all.equal(res$error, err))


#--- dimension mismatch errors -------------------------------------------------

# t length must match the observations of x
err <- '"t" must be a numeric vector with length ncol(x)'
res <- capture(row_cosinor(1:5, 1:3))
stopifnot(all.equal(res$error, err))

# period must be a single number
err <- '"period" must be a numeric vector with length 1'
res <- capture(row_cosinor(1:5, 1:5, 1:2))
stopifnot(all.equal(res$error, err))

