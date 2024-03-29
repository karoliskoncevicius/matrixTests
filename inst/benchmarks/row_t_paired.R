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
    Y <- matrix(rnorm(nr*nc), nrow=nr, ncol=nc)
    resM[as.character(nr), as.character(nc)] <- median(replicate(nreps, system.time(row_t_paired(X, Y))[3]))
    resB[as.character(nr), as.character(nc)] <- median(replicate(nreps, system.time(for(i in seq.int(nr)) t.test(X[i,], Y[i,], paired=TRUE))[3]))
    cat("rows:", nr, "cols:", nc, "base:", round(resB[as.character(nr),as.character(nc)],3), "(s)", "matrixTests:", round(resM[as.character(nr),as.character(nc)],3), "(s)", "\n")
  }
}

#--- plot ----------------------------------------------------------------------

png("row_t_paired.png", width=1000, height=800)
plotBenchResults(resB, resM, "t.test(x, y, paired=TRUE)", "row_t_paired(x, y)")
invisible(dev.off())

