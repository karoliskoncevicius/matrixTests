library(matrixTests)
source("utils/capture.r")

#--- x has less than 1 observation ---------------------------------------------

wrn <- 'row_kolmogorovsmirnov_twosample: 1 of the rows had less than 1 "x" observation.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# 0 observations
res <- capture(row_kolmogorovsmirnov_twosample(numeric(), 1))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.tot, 1))

# only NA observations
res <- capture(row_kolmogorovsmirnov_twosample(c(NA_integer_, NA_integer_), 1))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.tot, 1))


#--- y has less than 1 observations --------------------------------------------

wrn <- 'row_kolmogorovsmirnov_twosample: 1 of the rows had less than 1 "y" observation.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# 0 observations
res <- capture(row_kolmogorovsmirnov_twosample(1, numeric()))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 1))

# only NA observations
res <- capture(row_kolmogorovsmirnov_twosample(1, c(NA_integer_, NA_integer_)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 1))


#--- both x and y have less than 1 observation ---------------------------------

wrnx <- 'row_kolmogorovsmirnov_twosample: 1 of the rows had less than 1 "x" observation.\nFirst occurrence at row 1'
wrny <- 'row_kolmogorovsmirnov_twosample: 1 of the rows had less than 1 "y" observation.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# 0 observations
res <- capture(row_kolmogorovsmirnov_twosample(numeric(), numeric()))
stopifnot(all.equal(res$warning[1:2], c(wrnx, wrny)))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 0))

# only NA observations
res <- capture(row_kolmogorovsmirnov_twosample(c(NA_integer_, NA_integer_), c(NA_integer_, NA_integer_)))
stopifnot(all.equal(res$warning[1:2], c(wrnx, wrny)))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 0))


#--- ties are present ----------------------------------------------------------

wrn <- 'row_kolmogorovsmirnov_twosample: 1 of the rows had ties: asymptotic p-values will be approximate.\nFirst occurrence at row 1'

# warning when exact=TRUE
res <- capture(row_kolmogorovsmirnov_twosample(c(3,3,1,2), c(2,0,5,-2), exact=FALSE))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$exact, FALSE))
stopifnot(all.equal(res$value$obs.x, 4))
stopifnot(all.equal(res$value$obs.y, 4))
stopifnot(all.equal(res$value$obs.tot, 8))

# no warning when exact=FALSE
res <- capture(row_kolmogorovsmirnov_twosample(c(3,3,1,2), c(2,0,5,-2), exact=TRUE))
stopifnot(all.equal(res$warning, NULL))
stopifnot(all.equal(res$value$exact, TRUE))
stopifnot(all.equal(res$value$obs.x, 4))
stopifnot(all.equal(res$value$obs.y, 4))
stopifnot(all.equal(res$value$obs.tot, 8))

