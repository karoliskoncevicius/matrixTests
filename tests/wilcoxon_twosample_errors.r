library(matrixTests)
source("utils/capture.r")

#--- x argument errors ---------------------------------------------------------

# cannot be missing
err <- 'argument "x" is missing, with no default'
res <- capture(row_wilcoxon_twosample())
stopifnot(all.equal(res$error, err))

# cannot be NULL
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(NULL, 1:2))
stopifnot(all.equal(res$error, err))

# cannot be character
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(c("1", "2"), 1:2))
stopifnot(all.equal(res$error, err))

# cannot be logical
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(c(TRUE, FALSE), 1:2))
stopifnot(all.equal(res$error, err))

# cannot be complex
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(complex(c(1,2), c(3,4)), 1:2))
stopifnot(all.equal(res$error, err))

# cannot be data.frame containing some non numeric data
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(iris, 1:2))
stopifnot(all.equal(res$error, err))

# cannot be a list
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(as.list(c(1:5)), 1:2))
stopifnot(all.equal(res$error, err))

# cannot be in a list
err <- '"x" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(list(1:5), 1:2))
stopifnot(all.equal(res$error, err))


#--- y argument errors ---------------------------------------------------------

# cannot be missing
err <- 'argument "y" is missing, with no default'
res <- capture(row_wilcoxon_twosample(1))
stopifnot(all.equal(res$error, err))

# cannot be NULL
err <- '"y" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(1, NULL))
stopifnot(all.equal(res$error, err))

# cannot be character
err <- '"y" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(1:2, c("1","2")))
stopifnot(all.equal(res$error, err))

# cannot be logical
err <- '"y" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(1:2, c(TRUE, FALSE)))
stopifnot(all.equal(res$error, err))

# cannot be complex
err <- '"y" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(1:2, complex(c(1,2), c(3,4))))
stopifnot(all.equal(res$error, err))

# cannot be data.frame containing some non numeric data
err <- '"y" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(1:2, iris))
stopifnot(all.equal(res$error, err))

# cannot be a list
err <- '"y" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(1:2, as.list(c(1:5))))
stopifnot(all.equal(res$error, err))

# cannot be in a list
err <- '"y" must be a numeric matrix or vector'
res <- capture(row_wilcoxon_twosample(1:2, list(1:5)))
stopifnot(all.equal(res$error, err))


#--- alternative argument errors -----------------------------------------------

err <- '"alternative" must be a character vector with length 1 or nrow(x)'

# cannot be NA
res <- capture(row_wilcoxon_twosample(x=1:3, y=2:4, alternative=NA))
stopifnot(all.equal(res$error, err))

# cannot be numeric
res <- capture(row_wilcoxon_twosample(x=1:3, y=2:4, alternative=1))
stopifnot(all.equal(res$error, err))

# cannot be complex
res <- capture(row_wilcoxon_twosample(x=1:3, y=2:4, alternative=complex(1)))
stopifnot(all.equal(res$error, err))

# cannot be in a list
res <- capture(row_wilcoxon_twosample(x=1:3, y=2:4, alternative=list("less")))
stopifnot(all.equal(res$error, err))

# cannot be a data frame
res <- capture(row_wilcoxon_twosample(x=1:3, y=2:4, alternative=data.frame("less")))
stopifnot(all.equal(res$error, err))


err <- 'all "alternative" values must be in: two.sided, less, greater'

# must be in correct set
res <- capture(row_wilcoxon_twosample(x=1:3, y=2:4, alternative="ga"))
stopifnot(all.equal(res$error, err))

# error produced even when some are correct
res <- capture(row_wilcoxon_twosample(x=matrix(rnorm(10), nrow=2), y=matrix(rnorm(10), nrow=2), alternative=c("g","c")))
stopifnot(all.equal(res$error, err))


#--- mu argument errors --------------------------------------------------------

err <- '"mu" must be a numeric vector with length 1 or nrow(x)'

# cannot be a character
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, mu="1"))
stopifnot(all.equal(res$error, err))

# cannot be complex
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, mu=complex(1)))
stopifnot(all.equal(res$error, err))

# cannot be in a list
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, mu=list(1)))
stopifnot(all.equal(res$error, err))

# cannot be a data frame
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, mu=data.frame(1)))
stopifnot(all.equal(res$error, err))


err <- 'all "mu" values must be greater than -Inf and lower than Inf'

# cannot be NA
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, mu=NA_integer_))
stopifnot(all.equal(res$error, err))

# cannot be NaN
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, mu=NaN))
stopifnot(all.equal(res$error, err))

# TODO: check if can be made to work with Inf
# cannot be Inf
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, mu=Inf))
stopifnot(all.equal(res$error, err))

# cannot be -Inf
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, mu=-Inf))
stopifnot(all.equal(res$error, err))


#--- exact argument errors -----------------------------------------------------

err <- '"exact" must be a logical vector with length 1 or nrow(x)'

# cannot be non-logical NA
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, exact=NA_integer_))
stopifnot(all.equal(res$error, err))

# cannot be numeric
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, exact=1))
stopifnot(all.equal(res$error, err))

# cannot be character
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, exact="TRUE"))
stopifnot(all.equal(res$error, err))

# cannot be complex
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, exact=complex(1)))
stopifnot(all.equal(res$error, err))

# cannot be in a list
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, exact=list(TRUE)))
stopifnot(all.equal(res$error, err))

# cannot be a data frame
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, exact=data.frame(TRUE)))
stopifnot(all.equal(res$error, err))


#--- correct argument errors ---------------------------------------------------

err <- 'all "correct" values must be in: TRUE, FALSE'

# cannot be NA
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, correct=NA))
stopifnot(all.equal(res$error, err))


err <- '"correct" must be a logical vector with length 1 or nrow(x)'

# cannot be numeric
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, correct=0))
stopifnot(all.equal(res$error, err))

# cannot be character
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, correct="FALSE"))
stopifnot(all.equal(res$error, err))

# cannot be complex
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, correct=complex(1)))
stopifnot(all.equal(res$error, err))

# cannot be in a list
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, correct=list(FALSE)))
stopifnot(all.equal(res$error, err))

# cannot be a data frame
res <- capture(row_wilcoxon_twosample(x=1:3, y=1:3, correct=data.frame(FALSE)))
stopifnot(all.equal(res$error, err))


#--- dimension mismatch errors -------------------------------------------------

# y number of rows must match x number of rows
err <- '"x" and "y" must have the same number of rows'
x <- matrix(rnorm(10), nrow=2)
y <- matrix(rnorm(10), nrow=5)
res <- capture(row_wilcoxon_twosample(x, y))
stopifnot(all.equal(res$error, err))

# mu must match x number of rows
err <- '"mu" must be a numeric vector with length 1 or nrow(x)'
x <- matrix(rnorm(12), nrow=4)
y <- matrix(rnorm(12), nrow=4)
res <- capture(row_wilcoxon_twosample(x, y, mu=c(1,2)))
stopifnot(all.equal(res$error, err))

# alternative must match x number of rows
err <- '"alternative" must be a character vector with length 1 or nrow(x)'
x <- matrix(rnorm(12), nrow=4)
y <- matrix(rnorm(12), nrow=4)
res <- capture(row_wilcoxon_twosample(x, y, alternative=c("g","l")))
stopifnot(all.equal(res$error, err))

# exact must match x number of rows
err <- '"exact" must be a logical vector with length 1 or nrow(x)'
x <- matrix(rnorm(12), nrow=4)
y <- matrix(rnorm(12), nrow=4)
res <- capture(row_wilcoxon_twosample(x, y, exact=c(TRUE, FALSE)))
stopifnot(all.equal(res$error, err))

# correct must match x number of rows
err <- '"correct" must be a logical vector with length 1 or nrow(x)'
x <- matrix(rnorm(12), nrow=4)
y <- matrix(rnorm(12), nrow=4)
res <- capture(row_wilcoxon_twosample(x, y, correct=c(TRUE, FALSE)))
stopifnot(all.equal(res$error, err))

