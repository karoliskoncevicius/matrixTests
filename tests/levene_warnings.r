library(matrixTests)
source("utils/capture.r")

#--- removing na groups --------------------------------------------------------

wrn <- '2 columns dropped due to missing group information'

# 2 NAs
x <- 1:10
g <- c(1,1,1,1,NA,NA,2,2,2,2)
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 8))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- less than 2 groups --------------------------------------------------------

wrn <- 'row_levene: 1 of the rows had less than 2 groups with enough observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# all values in one group
x <- rnorm(10)
g <- rep("a", 10)
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 10))
stopifnot(all.equal(res$value$obs.groups, 1))

# two groups but one has only NAs
x <- c(rnorm(3), NA, NA, NA)
g <- rep(c("a","b"), each=3)
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 3))
stopifnot(all.equal(res$value$obs.groups, 1))


#--- no groups with 3 observations ---------------------------------------------

wrn <- 'row_levene: 1 of the rows had no groups with at least 3 observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# 10 groups all with a single value
x <- 1:10
g <- letters[1:10]
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 10))
stopifnot(all.equal(res$value$obs.groups, 10))

# 5 groups 2 observations each
x <- 1:10
g <- rep(letters[1:5], each=2)
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 10))
stopifnot(all.equal(res$value$obs.groups, 5))

# two groups one with 3 observations but one of them is NA
x <- c(rnorm(3), NA)
g <- c("a", "b", "b","b")
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 3))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- close to constant variance ------------------------------------------------

wrn <- 'row_levene: 1 of the rows had essentially constant absolute residuals from the mean: results might be unreliable.\nFirst occurrence at row 1'

# two groups - close to constant values within group
x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
       1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
       )
g <- rep(c("a","b"), each=4)
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 8))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- constant variance after residuals -----------------------------------------

wrn <- 'row_levene: 1 of the rows had zero within group variance of absolute residuals from the mean.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# two groups - all zero values
x <- rep(0,6)
g <- rep(letters[1:2], each=3)
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 6))
stopifnot(all.equal(res$value$obs.groups, 2))

# two groups - all constant values
x <- rep(1,6)
g <- rep(letters[1:2], each=3)
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 6))
stopifnot(all.equal(res$value$obs.groups, 2))

# three groups - constant values plus NAs
x <- c(1,1,1,NA,2,2,2,NA,3,3,3,NA)
g <- rep(1:3, each=4)
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 9))
stopifnot(all.equal(res$value$obs.groups, 3))

# values become constant after residuals
x <- c(1,1,2,2,3,3,4,4,5,5,6,6)
g <- rep(letters[1:3], each=4)
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 12))
stopifnot(all.equal(res$value$obs.groups, 3))

# values become constant after residuals with different groups sizes
x <- c(1,1,1,2,2,2,3,3,4,4,5,6,7)
g <- rep(letters[1:4], c(6,4,2,1))
res <- capture(row_levene(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 13))
stopifnot(all.equal(res$value$obs.groups, 4))

