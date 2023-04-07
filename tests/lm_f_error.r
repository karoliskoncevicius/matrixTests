library(matrixTests)
source("utils/capture.r")

#--- x argument errors ---------------------------------------------------------

# cannot be missing
err <- 'argument "x" is missing, with no default'
res <- capture(row_lm_f())
stopifnot(all.equal(res$error, err))

# cannot be NULL
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(NULL, "a"))
stopifnot(all.equal(res$error, err))

# cannot be character
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(c("a", "b"), c("a","b")))
stopifnot(all.equal(res$error, err))

# cannot be logical
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(c(TRUE, FALSE), "a"))
stopifnot(all.equal(res$error, err))

# cannot be complex
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(complex(c(1,2), c(3,4)), "a"))
stopifnot(all.equal(res$error, err))

# cannot be data.frame containing some non numeric data
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(iris, "a"))
stopifnot(all.equal(res$error, err))

# cannot be a list
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(as.list(c(1:5)), "a"))
stopifnot(all.equal(res$error, err))

# cannot be in a list
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(list(1:5), "a"))
stopifnot(all.equal(res$error, err))


#--- m argument errors ---------------------------------------------------------

#TODO: doublecheck

# cannot be missing
err <- 'argument "m" is missing, with no default'
res <- capture(row_lm_f(1))
stopifnot(all.equal(res$error, err))

# cannot be NULL
err <- '"m" must be a numeric matrix'
res <- capture(row_lm_f(1, NULL))
stopifnot(all.equal(res$error, err))

# cannot be character
err <- '"m" must be a numeric matrix'
res <- capture(row_lm_f(matrix(1:4, nrow=1), c("1","2","3","4")))
stopifnot(all.equal(res$error, err))

# cannot be logical
err <- '"m" must be a numeric matrix'
res <- capture(row_lm_f(1, TRUE))
stopifnot(all.equal(res$error, err))

# cannot be complex
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(complex(c(1,2), c(3,4)), "a"))
stopifnot(all.equal(res$error, err))

# cannot be data.frame containing some non numeric data
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(iris, "a"))
stopifnot(all.equal(res$error, err))

# cannot be a list
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_lm_f(as.list(c(1:5)), "a"))
stopifnot(all.equal(res$error, err))
# cannot be a in a list
err <- '"m" must be a numeric matrix'
res <- capture(row_lm_f(1:5, list(1:5)))
stopifnot(all.equal(res$error, err))


#--- null argument errors ------------------------------------------------------



#--- dimension mismatch errors -------------------------------------------------


#--- nesting errors ------------------------------------------------------------
