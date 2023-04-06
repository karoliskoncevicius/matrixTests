library(matrixTests)

#--- special argument cases ----------------------------------------------------

# x can be a vector
x <- rnorm(9)
X <- matrix(x, nrow=1)
m <- model.matrix(~ abs(x))
stopifnot(all.equal(row_lm_f(x, m), row_lm_f(X, m)))

# x can be numeric data.frame
x <- iris[,1:4]
X <- as.matrix(x)
m <- model.matrix(~ abs(1:4))
stopifnot(all.equal(row_lm_f(x, m), row_lm_f(X, m)))

# x can have 0 rows
x <- matrix(0, nrow=0, ncol=4)
m <- model.matrix(~ abs(1:4))
stopifnot(all.equal(nrow(row_lm_f(x, m)), 0))

# x can have 0 rows and 0 columns
x <- matrix(0, nrow=0, ncol=0)
m <- model.matrix(~ numeric())
stopifnot(all.equal(nrow(row_lm_f(x, m)), 0))


#--- missing values ------------------------------------------------------------

# missing values in x are removed
x <- c(NA,NA,rnorm(7),NA)
m <- model.matrix(~ abs(1:length(x)))
res1 <- row_lm_f(x, m)
res2 <- row_lm_f(x[!is.na(x)], m[!is.na(x),])
stopifnot(all.equal(res1, res2))

# no difference between NA and NaN in x
x1 <- c(NA,NA,1:7,NA)
x2 <- c(NaN,NaN,1:7,NaN)
m  <- model.matrix(~ abs(1:length(x1)))
res1 <- row_lm_f(x1, m)
res2 <- row_lm_f(x2, m)
stopifnot(all.equal(res1, res2))

# all x values can be NA
x <- rep(NA_integer_, 4)
m <- model.matrix(~ abs(1:length(x)))
res <- suppressWarnings(row_lm_f(x, m))
stopifnot(all.equal(res$obs, 0))

# TODO: check if model values can be NA


#--- different model order -----------------------------------------------------

# model and null can have different column orders
m1 <- model.matrix(~ mpg + cyl + disp, data=mtcars)
m0 <- model.matrix(~       cyl + disp, data=mtcars)
res1 <- row_lm_f(t(mtcars), m1, m0)
res2 <- row_lm_f(t(mtcars), m1, m0[,c(3,2,1)])
stopifnot(all.equal(res1, res2))


# TODO: beta names


#--- rownames ------------------------------------------------------------------

# when not provided - numbers are used instead.
x <- matrix(rnorm(20), nrow=2)
m <- model.matrix(~ abs(1:ncol(x)))
res <- row_lm_f(x, m)
stopifnot(all.equal(rownames(res), c("1", "2")))

# when provided - preserved (matrix)
x <- matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B")))
m <- model.matrix(~ abs(1:ncol(x)))
res <- row_lm_f(x, m)
stopifnot(all.equal(rownames(res), rownames(x)))

# when provided - preserved (data.frame)
x <- data.frame(matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B"))))
m <- model.matrix(~ abs(1:ncol(x)))
res <- row_lm_f(x, m)
stopifnot(all.equal(rownames(res), rownames(x)))

# when duplicated - made unique
x <- matrix(rnorm(40), nrow=4, dimnames=list(c("A", "A", "B", "B")))
m <- model.matrix(~ abs(1:ncol(x)))
res <- row_lm_f(x, m)
stopifnot(all.equal(rownames(res), make.unique(rownames(x))))

