rowVars <- function(x, n=NULL, m=NULL, na.rm=FALSE) {
  if(is.null(n))
    n <- ncol(x) - matrixStats::rowCounts(x, value=NA)
  if(is.null(m))
    m <- rowMeans(x, na.rm=na.rm)
  res <- rowSums((x-m)^2, na.rm=na.rm) / (n-1)
  res[n <= 1] <- NA
  res[!is.finite(m)] <- NaN
  res
}

# NOTE: assumes that the input is a matrix with ranks in each row
#       the behaviour of this function is quite trikcy
#       it converts the matrix of ranks into an index-form representation of a different matrix,
#       where rows represent the original rows but columns correspond to ranks (values of r).
#       then it counts the occurance of each such index to get table of ranks for each row
rowRankTies <- function(r) {
  dupRows <- apply(r, 1, anyDuplicated, incomparables=NA) != 0
  if (any(dupRows)) {
    r <- r[dupRows,, drop=FALSE]
    storage.mode(r) <- "integer"

    inds  <- nrow(r) * (r - 1) + row(r)
    msize <- nrow(r) * (ncol(r) - 1) + nrow(r)
    t <- tabulate(inds, msize)

    res <- matrix(1L, nrow=length(dupRows), ncol=ncol(r))
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
