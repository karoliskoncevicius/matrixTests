#--- libraries -----------------------------------------------------------------

library(matrixTests)

source("R/plot.R", chdir=TRUE)

#--- parameters ----------------------------------------------------------------

ncols <- c(10,100,1000)
nrows <- 10^c(1:5)
nreps <- 5

#--- benchmark -----------------------------------------------------------------

resM <- matrix(nrow=length(nrows), ncol=length(ncols))
rownames(resM) <- nrows
colnames(resM) <- ncols
resB <- resM

for(nr in nrows) {
  for(nc in ncols) {
    X <- matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
    resM[as.character(nr), as.character(nc)] <- median(replicate(nreps, system.time(row_andersondarling(X))[3]))
    resB[as.character(nr), as.character(nc)] <- median(replicate(nreps, system.time(for(i in seq.int(nr)) nortest::ad.test(X[i,]))[3]))
    cat("rows:", nr, "cols:", nc, "base:", round(resB[as.character(nr),as.character(nc)],3), "(s)", "matrixTests:", round(resM[as.character(nr),as.character(nc)],3), "(s)", "\n")
  }
}

#--- plot ----------------------------------------------------------------------

png("row_andersondarling", width=1000, height=800)
plotBenchResults(resB, resM, "nortest::ad.test(x)", "row_andersondarling(x)")
invisible(dev.off())

