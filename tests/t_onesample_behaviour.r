library(matrixTests)

#--- special argument cases ----------------------------------------------------

# x and can be vector
x <- rnorm(9)
X <- matrix(x, nrow=1)
stopifnot(all.equal(row_t_onesample(x), row_t_onesample(X)))

# x and be numeric data frame
x <- iris[,1:4]
X <- as.matrix(x)
stopifnot(all.equal(row_t_onesample(x), row_t_onesample(X)))

# x can have 0 rows
x <- matrix(0, nrow=0, ncol=4)
stopifnot(all.equal(nrow(row_t_onesample(x)), 0))

# x can have 0 rows and 0 columns
x <- matrix(0, nrow=0, ncol=0)
stopifnot(all.equal(nrow(row_t_onesample(x)), 0))

# alternative values can be partially completed
x <- rnorm(5)
stopifnot(all.equal(row_t_onesample(x, alternative="greater"), row_t_onesample(x, alternative="g")))

# null can be +Inf
x <- rnorm(5)
res <- row_t_onesample(x, null=Inf)
stopifnot(all.equal(res$pvalue, 0))
stopifnot(all.equal(res$statistic, -Inf))

# null can be -Inf
x <- rnorm(5)
res <- row_t_onesample(x, null=-Inf)
stopifnot(all.equal(res$pvalue, 0))
stopifnot(all.equal(res$statistic, Inf))

#conf.level can be 0
x <- rnorm(5)
res <- row_t_onesample(x, conf.level=0)
stopifnot(all.equal(res$conf.low, res$conf.high))

# conf.level can be 1
x <- rnorm(5)
res <- row_t_onesample(x, conf.level=1)
stopifnot(all.equal(res$conf.low, -Inf))
stopifnot(all.equal(res$conf.high, Inf))


#--- recycling -----------------------------------------------------------------

# null can be specified for each row
x <- matrix(rnorm(10), nrow=2)
res1 <- row_t_onesample(x, null=2)
res2 <- row_t_onesample(x, null=c(2,2))
stopifnot(all.equal(res1, res2))

# alternative can be specified for each row
x <- matrix(rnorm(10), nrow=2)
res1 <- row_t_onesample(x, alternative="g")
res2 <- row_t_onesample(x, alternative=c("g","g"))
stopifnot(all.equal(res1, res2))

# conf.level can be specified for each row
x <- matrix(rnorm(10), nrow=2)
res1 <- row_t_onesample(x, conf.level=0.98)
res2 <- row_t_onesample(x, conf.level=c(0.98,0.98))
stopifnot(all.equal(res1, res2))


#--- missing values ------------------------------------------------------------

# missing values in x are removed
x <- c(NA,NA,rnorm(7),NA)
res1 <- row_t_onesample(x)
res2 <- row_t_onesample(x[!is.na(x)])
stopifnot(all.equal(res1, res2))

# no difference between NA and NaN
x1 <- c(NA,NA,1:7,NA)
x2 <- c(NaN,NaN,1:7,NaN)
res1 <- row_t_onesample(x1)
res2 <- row_t_onesample(x2)
stopifnot(all.equal(res1, res2))

# everything can be NA
x <- rep(NA_integer_, 4)
res <- suppressWarnings(row_t_onesample(x))
stopifnot(all.equal(res$obs, 0))


#--- rownames ------------------------------------------------------------------

# when not provided - numbers are used instead.
x <- matrix(rnorm(20), nrow=2)
res <- row_t_onesample(x)
stopifnot(all.equal(rownames(res), c("1", "2")))

# when provided - preserved (matrix)
x <- matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B")))
res <- row_t_onesample(x)
stopifnot(all.equal(rownames(res), rownames(x)))

# when provided - preserved (data.frame)
x <- data.frame(matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B"))))
res <- row_t_onesample(x)
stopifnot(all.equal(rownames(res), rownames(x)))

# when duplicated - made unique
x <- matrix(rnorm(40), nrow=4, dimnames=list(c("A", "A", "B", "B")))
res <- row_t_onesample(x)
stopifnot(all.equal(rownames(res), make.unique(rownames(x))))

