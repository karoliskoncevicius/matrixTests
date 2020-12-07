library(matrixTests)
source("utils/capture.r")

#--- removing na groups --------------------------------------------------------

wrn <- '2 columns dropped due to missing group information'

# 2 NAs
x <- 1:10
g <- c(1,1,1,1,NA,NA,2,2,2,2)
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 8))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- groups with one observations ----------------------------------------------

# TODO: In order to be consistent with tests that remove empty groups either:
#       1) Only show this warning in cases where group with 1 observation is removed.
#       2) Make other tests show similar warnings where empty groups are removed.
wrn <- 'row_oneway_welch: 1 of the rows had groups with less than 2 observations: those groups were removed.\nFirst occurrence at row 1'

# 3 groups with one having only one observation
x <- 1:5
g <- c("a","a","b","b","c")
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 4))
stopifnot(all.equal(res$value$obs.groups, 2))

# 3 groups with 1 having one observation due to NA values
x <- c(1:4, NA, 1)
g <- c("a","a","b","b","c","c")
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 4))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- less than 2 groups --------------------------------------------------------

wrn <- 'row_oneway_welch: 1 of the rows had less than 2 groups with enough observations.\nFirst occurrence at row 1'
nacolumns <- c("df.between", "df.within", "statistic", "pvalue")

# all values in one group
x <- 1:10
g <- rep("a", 10)
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 10))
stopifnot(all.equal(res$value$obs.groups, 1))

# many groups but all have one observation
x <- 1:10
g <- letters[1:10]
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning[2], wrn)) # TODO: fix the extra warning produced in this situation
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 0))
stopifnot(all.equal(res$value$obs.groups, 0))

# two groups but one has only a single observation and is removed
x <- 1:4
g <- c("a","a","a","b")
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 3))
stopifnot(all.equal(res$value$obs.groups, 1))

# two groups but one has only NAs
x <- c(1:3, NA, NA, NA)
g <- rep(c("a","b"), each=3)
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 3))
stopifnot(all.equal(res$value$obs.groups, 1))


#--- all groups have zero variance ---------------------------------------------

wrn <- 'row_oneway_welch: 1 of the rows had zero variance in all of the groups.\nFirst occurrence at row 1'
nacolumns <- c("df.between", "df.within", "statistic", "pvalue")

# two groups - all values are constant
x <- c(1,1,1,1)
g <- c("a","a","b","b")
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 4))
stopifnot(all.equal(res$value$obs.groups, 2))

# two groups - constant values within each group
x <- c(3,3,0,0)
g <- c("a","a","b","b")
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 4))
stopifnot(all.equal(res$value$obs.groups, 2))

# two groups - constant values plus NA within each group
x <- c(3,3,NA,0,0,NA)
g <- c("a","a","a","b","b","b")
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 4))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- one of the groups has zero variance ---------------------------------------

wrn <- 'row_oneway_welch: 1 of the rows had groups with zero variance: result might be unreliable.\nFirst occurrence at row 1'

# two groups - one with constant values
x <- c(1,2,1,1)
g <- c("a","a","b","b")
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$statistic, Inf))

# three groups - constant values within two of them
x <- c(1,2,3,3,5,5)
g <- c("a","a","b","b","c","c")
res <- capture(row_oneway_welch(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$statistic, Inf))

