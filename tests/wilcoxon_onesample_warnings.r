library(matrixTests)
source("utils/capture.r")

#--- x has less than 1 observation ---------------------------------------------

wrn <- 'row_wilcoxon_onesample: 1 of the rows had less than 1 remaining "x" observation.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# 0 observations
res <- capture(row_wilcoxon_onesample(numeric()))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 0))

# only NA observations
res <- capture(row_wilcoxon_onesample(c(NA_integer_,NA_integer_)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 0))


#--- values equal to null are removed ------------------------------------------

wrn <- 'row_wilcoxon_onesample: 1 of the rows had observations with "x" equal "null" that were removed.\nFirst occurrence at row 1'

# a few x values equal to null
res <- capture(row_wilcoxon_onesample(c(1,2,3,4,1), null=1, exact=FALSE))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 3))

# one x value equal to null plus NAs
res <- capture(row_wilcoxon_onesample(c(NA,3,4,1), null=1, exact=FALSE))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 2))

wrn1 <- 'row_wilcoxon_onesample: 1 of the rows had observations with "x" equal "null" that were removed.\nFirst occurrence at row 1'
wrn2 <- 'row_wilcoxon_onesample: 1 of the rows had zeroes: cannot compute exact p-values with zeroes.\nFirst occurrence at row 1'

# a few x values equal to null and exact is changed to FALSE
res <- capture(row_wilcoxon_onesample(c(1,2,3,4,1), null=1, exact=TRUE))
stopifnot(all.equal(res$warning, c(wrn1,wrn2)))
stopifnot(all.equal(res$value$obs, 3))
stopifnot(all.equal(res$value$exact, FALSE))


#--- ties are present ----------------------------------------------------------

wrn <- 'row_wilcoxon_onesample: 1 of the rows had ties: cannot compute exact p-values with ties.\nFirst occurrence at row 1'

# warning when exact=TRUE
res <- capture(row_wilcoxon_onesample(c(1,2,3,4,1), exact=TRUE))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 5))
stopifnot(all.equal(res$value$exact, FALSE))

# no warning when exact=FALSE
res <- capture(row_wilcoxon_onesample(c(1,2,3,4,1), exact=FALSE))
stopifnot(all.equal(res$warning, NULL))
stopifnot(all.equal(res$value$obs, 5))
stopifnot(all.equal(res$value$exact, FALSE))


#--- warning about zero values -------------------------------------------------

wrn1 <- 'row_wilcoxon_onesample: 1 of the rows had observations with "x" equal "null" that were removed.\nFirst occurrence at row 1'
wrn2 <- 'row_wilcoxon_onesample: 1 of the rows had zeroes: cannot compute exact p-values with zeroes.\nFirst occurrence at row 1'

# warning when exact=TRUE
res <- capture(row_wilcoxon_onesample(c(1,2,3,4,1), null=2, exact=TRUE))
stopifnot(all.equal(res$warning, c(wrn1,wrn2)))
stopifnot(all.equal(res$value$obs, 4))
stopifnot(all.equal(res$value$exact, FALSE))

# no warning about ties when exact=FALSE
res <- capture(row_wilcoxon_onesample(c(1,2,3,4,1), null=2, exact=FALSE))
stopifnot(all.equal(res$warning, wrn1))
stopifnot(all.equal(res$value$obs, 4))
stopifnot(all.equal(res$value$exact, FALSE))

# if both ties and zeroes are presnet - only warn about zero
res <- capture(row_wilcoxon_onesample(c(0,1,1,2,3), exact=TRUE))
stopifnot(all.equal(res$warning, c(wrn1, wrn2)))
stopifnot(all.equal(res$value$exact, FALSE))

