library(matrixTests)
source("utils/capture.r")

#--- less than 3 observations --------------------------------------------------

wrn <- 'row_t_equalvar: 1 of the rows had less than 3 total observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# 0 observations
res <- capture(row_t_equalvar(NA_integer_, NA_integer_))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 0))

# 1 observation in each group
res <- capture(row_t_equalvar(1, 2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 1))
stopifnot(all.equal(res$value$obs.y, 1))
stopifnot(all.equal(res$value$obs.tot, 2))

# 2 observations in x and 0 in y
res <- capture(row_t_equalvar(1:2, numeric()))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 2))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 2))

# 0 observations in x and 2 in y
res <- capture(row_t_equalvar(numeric(), 1:2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 2))
stopifnot(all.equal(res$value$obs.tot, 2))

# both have 2 observations but one is NA in each
res <- capture(row_t_equalvar(c(1,NA), c(NA,2)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 1))
stopifnot(all.equal(res$value$obs.y, 1))
stopifnot(all.equal(res$value$obs.tot, 2))


#--- one group has no observations ---------------------------------------------

wrn <- 'row_t_equalvar: 1 of the rows had zero "x" observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# no observations in x
res <- capture(row_t_equalvar(numeric(), 1:3))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 3))
stopifnot(all.equal(res$value$obs.tot, 3))

# only NA observations in x
res <- capture(row_t_equalvar(c(NA_integer_, NA_integer_), 1:3))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 3))
stopifnot(all.equal(res$value$obs.tot, 3))


wrn <- 'row_t_equalvar: 1 of the rows had zero "y" observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# no observations in y
res <- capture(row_t_equalvar(1:3, numeric()))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 3))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 3))

# only NA observations in y
res <- capture(row_t_equalvar(1:3, c(NA_integer_, NA_integer_)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 3))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 3))


#--- constant values -----------------------------------------------------------

wrn <- 'row_t_equalvar: 1 of the rows were essentially constant.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# all values are the same
res <- capture(row_t_equalvar(c(1,1,1,1), c(1,1,1,1)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.y, 0))
stopifnot(all.equal(res$value$var.pooled, 0))

# all values are the equal within each group
res <- capture(row_t_equalvar(c(1,1,1,1), c(2,2,2,2)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.y, 0))
stopifnot(all.equal(res$value$var.pooled, 0))

# all values are the equal within each group after NAs
res <- capture(row_t_equalvar(c(NA,1,1,1,1), c(2,2,2,2,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.y, 0))
stopifnot(all.equal(res$value$var.pooled, 0))

