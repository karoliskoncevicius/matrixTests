library(matrixTests)

#--- special argument cases ----------------------------------------------------

# x and can be vector
x <- rnorm(9)
X <- matrix(x, nrow=1)
stopifnot(all.equal(row_wilcoxon_onesample(x), row_wilcoxon_onesample(X)))

# x and be numeric data frame
x <- iris[,1:4]
X <- as.matrix(x)
stopifnot(all.equal(row_wilcoxon_onesample(x), row_wilcoxon_onesample(X)))

# x can have 0 rows
x <- matrix(0, nrow=0, ncol=4)
stopifnot(all.equal(nrow(row_wilcoxon_onesample(x)), 0))

# x can have 0 rows and 0 columns
x <- matrix(0, nrow=0, ncol=0)
stopifnot(all.equal(nrow(row_wilcoxon_onesample(x)), 0))

# alternative values can be partially completed
x <- rnorm(5)
stopifnot(all.equal(row_wilcoxon_onesample(x, alternative="greater"), row_wilcoxon_onesample(x, alternative="g")))

# exact can be NA
x <- rnorm(5)
stopifnot(all.equal(row_wilcoxon_onesample(x, exact=NA), row_wilcoxon_onesample(x, exact=TRUE)))

# TODO: investigate if null can be allowed Inf values.


#--- recycling -----------------------------------------------------------------

# null can be specified for each row
x <- matrix(rnorm(10), nrow=2)
res1 <- row_wilcoxon_onesample(x, null=2)
res2 <- row_wilcoxon_onesample(x, null=c(2,2))
stopifnot(all.equal(res1, res2))

# alternative can be specified for each row
x <- matrix(rnorm(10), nrow=2)
res1 <- row_wilcoxon_onesample(x, alternative="g")
res2 <- row_wilcoxon_onesample(x, alternative=c("g","g"))
stopifnot(all.equal(res1, res2))

# exact can be specified for each row
x <- matrix(rnorm(10), nrow=2)
res1 <- row_wilcoxon_onesample(x, exact=TRUE)
res2 <- row_wilcoxon_onesample(x, exact=c(TRUE, TRUE))
stopifnot(all.equal(res1, res2))

# correct can be specified for each row
x <- matrix(rnorm(10), nrow=2)
res1 <- row_wilcoxon_onesample(x, exact=FALSE, correct=FALSE)
res2 <- row_wilcoxon_onesample(x, exact=FALSE, correct=c(FALSE, FALSE))
stopifnot(all.equal(res1, res2))


#--- missing values ------------------------------------------------------------

# missing values in x are removed
x <- c(NA,NA,rnorm(7),NA)
res1 <- row_wilcoxon_onesample(x)
res2 <- row_wilcoxon_onesample(x[!is.na(x)])
stopifnot(all.equal(res1, res2))

# no difference between NA and NaN
x1 <- c(NA,NA,1:7,NA)
x2 <- c(NaN,NaN,1:7,NaN)
res1 <- row_wilcoxon_onesample(x1)
res2 <- row_wilcoxon_onesample(x2)
stopifnot(all.equal(res1, res2))

# everything can be NA
x <- rep(NA_integer_, 4)
res <- suppressWarnings(row_wilcoxon_onesample(x))
stopifnot(all.equal(res$obs, 0))


#--- rownames ------------------------------------------------------------------

# when not provided - numbers are used instead.
x <- matrix(rnorm(20), nrow=2)
res <- row_wilcoxon_onesample(x)
stopifnot(all.equal(rownames(res), c("1", "2")))

# when provided - preserved (matrix)
x <- matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B")))
res <- row_wilcoxon_onesample(x)
stopifnot(all.equal(rownames(res), rownames(x)))

# when provided - preserved (data.frame)
x <- data.frame(matrix(rnorm(20), nrow=2, dimnames=list(c("A", "B"))))
res <- row_wilcoxon_onesample(x)
stopifnot(all.equal(rownames(res), rownames(x)))

# when duplicated - made unique
x <- matrix(rnorm(40), nrow=4, dimnames=list(c("A", "A", "B", "B")))
res <- row_wilcoxon_onesample(x)
stopifnot(all.equal(rownames(res), make.unique(rownames(x))))

