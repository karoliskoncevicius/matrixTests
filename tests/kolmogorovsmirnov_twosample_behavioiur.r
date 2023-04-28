library(matrixTests)

#--- special argument cases ----------------------------------------------------

# x and y can be vectors
x <- rnorm(9)
y <- rnorm(9)
X <- matrix(x, nrow=1)
Y <- matrix(y, nrow=1)
stopifnot(all.equal(row_kolmogorovsmirnov_twosample(x, y), row_kolmogorovsmirnov_twosample(X, Y)))

# x and y can be numeric data frames
x <- iris[1:50,1:4]
y <- iris[51:100,1:4]
X <- as.matrix(x)
Y <- as.matrix(y)
stopifnot(all.equal(suppressWarnings(row_kolmogorovsmirnov_twosample(x, y)), suppressWarnings(row_kolmogorovsmirnov_twosample(X, Y))))

# x and y can have 0 rows
x <- matrix(0, nrow=0, ncol=4)
y <- matrix(0, nrow=0, ncol=4)
stopifnot(all.equal(nrow(row_kolmogorovsmirnov_twosample(x, y)), 0))

# x and y can have 0 rows and 0 columns
x <- matrix(0, nrow=0, ncol=0)
y <- matrix(0, nrow=0, ncol=0)
stopifnot(all.equal(nrow(row_kolmogorovsmirnov_twosample(x, y)), 0))

# alternative values can be partially completed
x <- rnorm(5)
y <- rnorm(5)
stopifnot(all.equal(row_kolmogorovsmirnov_twosample(x, y, alternative="greater"), row_kolmogorovsmirnov_twosample(x, y, alternative="g")))

# exact can be NA
x <- rnorm(5)
y <- rnorm(5)
stopifnot(all.equal(row_kolmogorovsmirnov_twosample(x, y, exact=NA), row_kolmogorovsmirnov_twosample(x, y, exact=TRUE)))

# TODO: investigate if null can be allowed Inf values.


#--- recycling -----------------------------------------------------------------

# y can be a vector when x is a matrix
x <- matrix(rnorm(10), nrow=2)
y <- 1:5
res1 <- row_kolmogorovsmirnov_twosample(x, y)
res2 <- row_kolmogorovsmirnov_twosample(x, rbind(y,y))
stopifnot(all.equal(res1, res2))

# alternative can be specified for each row
x <- matrix(rnorm(10), nrow=2)
y <- matrix(rnorm(10), nrow=2)
res1 <- row_kolmogorovsmirnov_twosample(x, y, alternative="g")
res2 <- row_kolmogorovsmirnov_twosample(x, y, alternative=c("g","g"))
stopifnot(all.equal(res1, res2))

# exact can be specified for each row
x <- matrix(rnorm(10), nrow=2)
y <- matrix(rnorm(10), nrow=2)
res1 <- row_kolmogorovsmirnov_twosample(x, y, exact=TRUE)
res2 <- row_kolmogorovsmirnov_twosample(x, y, exact=c(TRUE, TRUE))
stopifnot(all.equal(res1, res2))


#--- missing values ------------------------------------------------------------

# missing values in x are removed
x <- c(NA,NA,rnorm(7),NA)
y <- rnorm(10)
res1 <- row_kolmogorovsmirnov_twosample(x, y)
res2 <- row_kolmogorovsmirnov_twosample(x[!is.na(x)], y)
stopifnot(all.equal(res1, res2))

# missing values in y are removed
x <- rnorm(10)
y <- c(NA,NA,rnorm(7),NA)
res1 <- row_kolmogorovsmirnov_twosample(x, y)
res2 <- row_kolmogorovsmirnov_twosample(x, y[!is.na(y)])
stopifnot(all.equal(res1, res2))

# missing values from both x and y are removed
x <- c(NA,NA,rnorm(7),NA)
y <- c(NA,rnorm(7),NA,NA)
res1 <- row_kolmogorovsmirnov_twosample(x, y)
res2 <- row_kolmogorovsmirnov_twosample(x[!is.na(x)], y[!is.na(y)])
stopifnot(all.equal(res1, res2))

# no difference between NA and NaN
x1 <- c(NA,NA,1,2,7:3,NA)
x2 <- c(NaN,NaN,1,2,7:3,NaN)
y1 <- c(NA,NA,18:12,NA)
y2 <- c(NaN,NaN,18:12,NaN)
res1 <- row_kolmogorovsmirnov_twosample(x1, y1)
res2 <- row_kolmogorovsmirnov_twosample(x2, y2)
stopifnot(all.equal(res1, res2))

# everything can be NA
x <- rep(NA_integer_, 4)
y <- rep(NA_integer_, 4)
res <- suppressWarnings(row_kolmogorovsmirnov_twosample(x, y))
stopifnot(all.equal(res$obs.tot, 0))


#--- rownames ------------------------------------------------------------------

# when not provided - numbers are used instead.
x <- matrix(rnorm(20), nrow=2)
y <- matrix(rnorm(20), nrow=2)
res <- row_kolmogorovsmirnov_twosample(x, y)
stopifnot(all.equal(rownames(res), c("1", "2")))

# when not provided - not taken from y
x <- matrix(rnorm(20), nrow=2)
y <- matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B")))
res <- row_kolmogorovsmirnov_twosample(x, y)
stopifnot(all.equal(rownames(res), c("1","2")))

# when provided - preserved (matrix)
x <- matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B")))
y <- matrix(rnorm(20), nrow=2)
res <- row_kolmogorovsmirnov_twosample(x, y)
stopifnot(all.equal(rownames(res), rownames(x)))

# when provided - preserved (data.frame)
x <- data.frame(matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B"))))
y <- matrix(rnorm(20), nrow=2)
res <- row_kolmogorovsmirnov_twosample(x, y)
stopifnot(all.equal(rownames(res), rownames(x)))

# when duplicated - made unique
x <- matrix(rnorm(20), nrow=4, dimnames=list(c("A", "A", "B", "B")))
y <- matrix(rnorm(20), nrow=4)
res <- row_kolmogorovsmirnov_twosample(x, y)
stopifnot(all.equal(rownames(res), make.unique(rownames(x))))

