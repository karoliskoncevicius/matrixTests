library(matrixTests)
source("utils/capture.r")

#--- less than 2 observations --------------------------------------------------

wrn <- 'row_andersondarling: 1 of the rows had less than 8 total observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# single observation
res <- capture(row_andersondarling(1))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 1))

# seven observations with some NA values
res <- capture(row_andersondarling(c(rnorm(7),NA,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 7))

# zero observations
res <- capture(row_andersondarling(NA_integer_))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 0))


#--- all values are constant ---------------------------------------------------

wrn <- 'row_andersondarling: 1 of the rows had essentially constant values.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# 8 equal observations
res <- capture(row_andersondarling(rep(1,8)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 8))

# 8 equal observations with some NA values
res <- capture(row_andersondarling(c(rep(0.5,8),NA,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 8))

