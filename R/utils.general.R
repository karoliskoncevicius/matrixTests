rowVars <- function(x, n=NULL, m=NULL, na.rm=FALSE) {
  if(is.null(n))
    n <- rep.int(ncol(x), nrow(x)) - rowSums(is.na(x))
  if(is.null(m))
    m <- rowMeans(x, na.rm=na.rm)
  res <- rowSums((x-m)^2, na.rm=na.rm) / (n-1)
  res[n <= 1] <- NA
  res[!is.finite(m)] <- NaN
  res
}

rowTies <- function(x) {
  dupRows <- apply(x, 1, anyDuplicated, incomparables=NA) != 0
  if (any(dupRows)) {
    # convert the values of each row to ranks
    x <- matrixStats::rowRanks(x[dupRows,, drop=FALSE], ties.method="dense")

    # construct a matrix where ranks are treated as columns
    # in this way column counts will track the number of repeated values[
    inds  <- nrow(x) * (x - 1) + row(x)
    msize <- nrow(x) * (ncol(x) - 1) + nrow(x)
    t <- tabulate(inds, msize)

    res <- matrix(1L, nrow=length(dupRows), ncol=ncol(x))
    res[dupRows,][1:msize] <- t
  } else {
    res <- matrix(1L, nrow=length(dupRows), ncol=1)
  }
  res
}

showWarning <- function(isWarning, err) {
  if(any(isWarning, na.rm=TRUE)) {
    parentFun <- deparse(as.list(sys.call(-1))[[1]])
    grandFun  <- as.list(sys.call(-2))
    if(length(grandFun) > 0) {
      grandFun <- deparse(grandFun[[1]])
      if(grandFun %in% getNamespaceExports("matrixTests")) {
        parentFun <- grandFun
      }
    }
    pref <- "row"
    if(grepl("^col_", parentFun)) pref <- "column"
    n <- sum(isWarning, na.rm=TRUE)
    i <- match(TRUE, isWarning)
    err <- paste0(parentFun, ": ", n, ' of the ', pref, 's ', err, ".",
                  '\nFirst occurrence at ', pref, ' ', i
                  )
    warning(err, call.=FALSE)
  }
}
