library(matrixTests)
source("utils/capture.r")

#--- less than 2 observations --------------------------------------------------

wrn <- 'row_jarquebera: 1 of the rows had less than 2 total observations.\nFirst occurrence at row 1'
nacolumns <- c("skewness", "kurtosis", "statistic", "pvalue")

# single observation
res <- capture(row_jarquebera(1))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 1))
stopifnot(all.equal(res$value$df, 2))

# single observation with some NA values
res <- capture(row_jarquebera(c(1,NA,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 1))
stopifnot(all.equal(res$value$df, 2))

# zero observations
res <- capture(row_jarquebera(NA_integer_))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 0))
stopifnot(all.equal(res$value$df, 2))


#--- all values are constant ---------------------------------------------------

wrn <- 'row_jarquebera: 1 of the rows had essentially constant values.\nFirst occurrence at row 1'
nacolumns <- c("skewness", "kurtosis", "statistic", "pvalue")

# two equal observations
res <- capture(row_jarquebera(c(1,1)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 2))
stopifnot(all.equal(res$value$df, 2))

# three equal observations with some NA values
res <- capture(row_jarquebera(c(0,0,0,NA,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 3))
stopifnot(all.equal(res$value$df, 2))

