library(matrixTests)
source("utils/capture.r")

#--- hes less than 1 paired observation ----------------------------------------

wrn <- 'row_wilcoxon_paired: 1 of the rows had less than 1 remaining paired "x-y" observation.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# no observations
res <- capture(row_wilcoxon_paired(numeric(), numeric()))
stopifnot(all.equal(res$warning[2], wrn))  # TODO: check why two warnings appear
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.paired, 0))

# x only has NAs
res <- capture(row_wilcoxon_paired(c(NA, NaN, NA), 1:3))
stopifnot(all.equal(res$warning[2], wrn))  # TODO: check why two warnings appear
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 3))
stopifnot(all.equal(res$value$obs.paired, 0))

# y only has NAs
res <- capture(row_wilcoxon_paired(1:3, c(NA, NaN, NA)))
stopifnot(all.equal(res$warning[2], wrn))  # TODO: check why two warnings appear
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 3))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.paired, 0))

# both x and y only NA observations
res <- capture(row_wilcoxon_paired(c(NA,NaN,NA), c(NA, NA, NaN)))
stopifnot(all.equal(res$warning[2], wrn))  # TODO: check why two warnings appear
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 0))
stopifnot(all.equal(res$value$obs.paired, 0))

# no common complete observations
res <- capture(row_wilcoxon_paired(c(1, NA, 3, NA), c(NA, 2, NA, 4)))
stopifnot(all.equal(res$warning[2], wrn))  # TODO: check why two warnings appear
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.x, 2))
stopifnot(all.equal(res$value$obs.y, 2))
stopifnot(all.equal(res$value$obs.paired, 0))


#--- ties are present ----------------------------------------------------------

wrn <- 'row_wilcoxon_paired: 1 of the rows had ties: cannot compute exact p-values with ties.\nFirst occurrence at row 1'

# warning when exact=TRUE
res <- capture(row_wilcoxon_paired(c(3,3,1,2), c(2,2,5,-2), exact=TRUE))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.x, 4))
stopifnot(all.equal(res$value$obs.y, 4))
stopifnot(all.equal(res$value$obs.paired, 4))

# no warning when exact=FALSE
res <- capture(row_wilcoxon_paired(c(3,3,1,2), c(2,2,5,-2), exact=FALSE))
stopifnot(all.equal(res$warning, NULL))
stopifnot(all.equal(res$value$obs.x, 4))
stopifnot(all.equal(res$value$obs.y, 4))
stopifnot(all.equal(res$value$obs.paired, 4))


#--- warning about zero values -------------------------------------------------

wrn1 <- 'row_wilcoxon_paired: 1 of the rows had observations with "x-y" equal "null" that were removed.\nFirst occurrence at row 1'
wrn2 <- 'row_wilcoxon_paired: 1 of the rows had zeroes: cannot compute exact p-values with zeroes.\nFirst occurrence at row 1'

# warning when exact = TRUE
res <- capture(row_wilcoxon_paired(1:4, c(0,2,4,3), null=1, exact=TRUE))
stopifnot(all.equal(res$warning, c(wrn1, wrn2)))
stopifnot(all.equal(res$value$obs.x, 4))
stopifnot(all.equal(res$value$obs.y, 4))
stopifnot(all.equal(res$value$obs.paired, 2))

# no warning when exact = FALSE
res <- capture(row_wilcoxon_paired(1:4, c(0,2,4,3), null=1, exact=FALSE))
stopifnot(all.equal(res$warning, wrn1))
stopifnot(all.equal(res$value$obs.x, 4))
stopifnot(all.equal(res$value$obs.y, 4))
stopifnot(all.equal(res$value$obs.paired, 2))

# if both ties and zeroes are presnet - only warn about zero
res <- capture(row_wilcoxon_paired(1:5, c(1,1,2,9,9), exact=TRUE))
stopifnot(all.equal(res$warning, c(wrn1, wrn2)))
stopifnot(all.equal(res$value$exact, FALSE))

