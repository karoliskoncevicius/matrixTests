library(matrixTests)
source("utils/capture.r")

#--- less than 2 x observations ------------------------------------------------

wrn <- 'row_t_welch: 1 of the rows had less than 2 "x" observations.\nFirst occurrence at row 1'
nacolumns <- c("stderr", "df", "statistic", "pvalue", "conf.low", "conf.high")

# 0 observations in both
res <- capture(row_t_welch(numeric(), numeric()))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 0))

# 0 observations in x
res <- capture(row_t_welch(numeric(), 1:3))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 3))
stopifnot(all.equal(res$value$obs.tot, 3))

# 1 observation in x
res <- capture(row_t_welch(1, 1:3))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 1))
stopifnot(all.equal(res$value$obs.y, 3))
stopifnot(all.equal(res$value$obs.tot, 4))

# 2 observation in x but one is NA
res <- capture(row_t_welch(c(1,NA), 1:3))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 1))
stopifnot(all.equal(res$value$obs.y, 3))
stopifnot(all.equal(res$value$obs.tot, 4))


#--- less than 2 y observations ------------------------------------------------

wrn <- 'row_t_welch: 1 of the rows had less than 2 "y" observations.\nFirst occurrence at row 1'
nacolumns <- c("stderr", "df", "statistic", "pvalue", "conf.low", "conf.high")

# 0 observations in y
res <- capture(row_t_welch(1:3, numeric()))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 3))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.tot, 3))

# 1 observation in y
res <- capture(row_t_welch(1:3, 1))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 3))
stopifnot(all.equal(res$value$obs.y, 1))
stopifnot(all.equal(res$value$obs.tot, 4))

# 2 observation in y but one is NA
res <- capture(row_t_welch(1:3, c(1,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 3))
stopifnot(all.equal(res$value$obs.y, 1))
stopifnot(all.equal(res$value$obs.tot, 4))


#--- constant values -----------------------------------------------------------

wrn <- 'row_t_welch: 1 of the rows had essentially constant values.\nFirst occurrence at row 1'
nacolumns <- c("stderr", "df", "statistic", "pvalue", "conf.low", "conf.high")

# all values are the same
res <- capture(row_t_welch(c(1,1,1,1), c(1,1,1,1)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.y, 0))

# all values are the equal within each group
res <- capture(row_t_welch(c(1,1,1,1), c(2,2,2,2)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.y, 0))

# all values are the equal within each group after NAs
res <- capture(row_t_welch(c(NA,1,1,1,1), c(2,2,2,2,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.y, 0))

