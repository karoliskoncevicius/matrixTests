################################################################################
################################## FUNCTIONS ###################################
################################################################################

plotBenchResults <- function(timesB, timesM, nameB, nameM) {

  nrows <- as.numeric(rownames(timesB))
  ncols <- as.numeric(colnames(timesB))

  timesD <- timesB/timesM
  timesD[timesD<1] <- -1/timesD[timesD<1]
  timesD[round(abs(timesD),6)==1] <- 0

  colsB <- colorRampPalette(c("lightgreen", "limegreen", "darkgreen"))(ncol(timesB))
  colsM <- colorRampPalette(c("lightblue", "cornflowerblue", "darkblue"))(ncol(timesM))
  colsD <- colorRampPalette(c("pink", "red", "darkred"))(ncol(timesD))

  layout(rbind(1:2, c(3,5), c(4,6)), heights=rbind(c(0.6, 0.6), c(0.2, 0.2), c(0.2, 0.2)))

  par(cex=1.5, cex.axis=0.7, cex.lab=0.8, cex.main=1, cex.sub=0.8, lab=c(5,5,7),
      las=1, mgp=c(2, 0.5, 0), tck=-0.015
      )

  plot.new()
  plot.window(xlim=range(log10(nrows)), ylim=range(pretty(c(timesB, timesM))))
  grid(nx=0, ny=NULL)
  Map(lines, asplit(timesB, 2), col=colsB, lwd=3)
  Map(lines, asplit(timesM, 2), col=colsM, lwd=3)
  axis(1, at=seq.int(nrows), labels=nrows)
  axis(2)
  title(xlab="number of rows")
  title(ylab="elapsed time (seconds)")
  title(main="raw speed")
  l1 <- legend("topleft", legend=paste(ncols, "columns"), lwd=2, col=c("grey90", "grey70", "grey50"), cex=0.7, bty="n")
  legend(l1$rect$left, l1$rect$top-l1$rect$h, legend=c(nameB, nameM), lwd=2, col=c(colsB[2], colsM[2]), cex=0.7, bty="n")

  plot.new()
  plot.window(xlim=range(log10(nrows)), ylim=range(pretty(timesD)))
  grid(nx=0, ny=NULL)
  Map(lines, asplit(timesD, 2), col=c("pink", "red", "darkred"), lwd=3)
  axis(1, at=seq.int(nrows), labels=nrows)
  axis(2)
  title(xlab="number of rows")
  title(ylab="fold speed change")
  title(main="fold speed difference")
  legend("topleft", legend=paste(ncols, "columns"), lwd=2, col=colsD, cex=0.7, bty="n")

  plotTable <- function(times, cols) {
    par(mar=c(2.1,4.1,2.1,2.1))
    plot.new()
    plot.window(xlim=c(1,nrow(times)), ylim=c(ncol(times), 0))
    text(row(times), col(times), format(times, digits=3, nsmall=3), xpd=TRUE)
    par(xpd=TRUE)
    segments(rep(0.3, ncol(times)), 1:ncol(times), rep(0.6, ncol(times)), 1:ncol(times), col=cols, lwd=3)
  }

  plotTable(timesB, colsB)
  title("base times (s)")

  plotTable(timesM, colsM)
  title("matrixTests times (s)")

  plotTable(timesD, colsD)
  title("fold speed increase")
  title(sub="(number of times matrixTests was faster than base)", line=1)
}

