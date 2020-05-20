library(matrixTests)
source("utils/capture.r")

#--- x has less than 2 observations --------------------------------------------

wrn <- 'row_f_var: 1 of the rows had less than 1 "x" observation.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# 0 observations
res <- capture(row_f_var(numeric(), 1:2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 2))
stopifnot(all.equal(res$value$obs.tot, 2))

# 1 observations
res <- capture(row_f_var(1, 1:2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 1))
stopifnot(all.equal(res$value$obs.y, 2))
stopifnot(all.equal(res$value$obs.tot, 3))

# 1 observation + NA
res <- capture(row_f_var(c(1, NA), 1:2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 1))
stopifnot(all.equal(res$value$obs.y, 2))
stopifnot(all.equal(res$value$obs.tot, 3))

# only NAs
res <- capture(row_f_var(c(NA_integer_, NA_integer_), 1:2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 2))
stopifnot(all.equal(res$value$obs.tot, 2))


#--- y has less than 2 observations --------------------------------------------

wrn <- 'row_f_var: 1 of the rows had less than 1 "y" observation.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# 0 observations
res <- capture(row_f_var(1:2, numeric()))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 2))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 2))

# 1 observations
res <- capture(row_f_var(1:2, 1))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 2))
stopifnot(all.equal(res$value$obs.y, 1))
stopifnot(all.equal(res$value$obs.tot, 3))

# 1 observation + NA
res <- capture(row_f_var(1:2, c(1, NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 2))
stopifnot(all.equal(res$value$obs.y, 1))
stopifnot(all.equal(res$value$obs.tot, 3))

# only NAs
res <- capture(row_f_var(1:2, c(NA_integer_, NA_integer_)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 2))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 2))


#--- x has 0 variance ----------------------------------------------------------

wrn <- 'row_f_var: 1 of the rows had zero variance in "x".\nFirst occurrence at row 1'

# constant values
res <- capture(row_f_var(rep(1,5), 1:5))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.ratio, 0))

# constant values with NAs
res <- capture(row_f_var(c(NA,rep(1,5),NA), 1:5))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.ratio, 0))


#--- y has 0 variance ----------------------------------------------------------

wrn <- 'row_f_var: 1 of the rows had zero variance in "y".\nFirst occurrence at row 1'

# constant values
res <- capture(row_f_var(1:5, rep(1,5)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$var.y, 0))
stopifnot(all.equal(res$value$var.ratio, Inf))

# constant values with NAs
res <- capture(row_f_var(1:5, c(NA,rep(1,5),NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$var.y, 0))
stopifnot(all.equal(res$value$var.ratio, Inf))


#--- both have 0 variance ------------------------------------------------------

wrnx <- 'row_f_var: 1 of the rows had zero variance in "x".\nFirst occurrence at row 1'
wrny <- 'row_f_var: 1 of the rows had zero variance in "y".\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# x and y have constant values
res <- capture(row_f_var(rep(1,5), rep(2,5)))
stopifnot(all.equal(res$warning, c(wrnx, wrny)))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.y, 0))
stopifnot(all.equal(res$value$var.ratio, NaN))

