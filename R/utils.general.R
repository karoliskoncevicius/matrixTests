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
  if(any(dupRows)) {
    dups <- matrix(FALSE, nrow=nrow(x), ncol=ncol(x))
    dups[dupRows,] <- t(apply(x[dupRows,,drop=FALSE], 1, duplicated, incomparables=NA))
    dups <- cbind(which(dups, arr.ind=TRUE), val=x[dups])
    dups <- dups[order(dups[,1], dups[,3]),,drop=FALSE]

    sp <- split(dups[,3], dups[,1])
    sp <- lapply(sp, function(x) rle(x)$length+1)
    cl <- lapply(sp, seq_along)

    dups <- cbind(rep(unique(dups[,1]), lengths(sp)),
                  unlist(cl),
                  unlist(sp)
                  )

    res <- matrix(0L, nrow=nrow(x), ncol=max(dups[,2]))
    res[dups[,1:2,drop=FALSE]] <- dups[,3]
  } else {
    res <- matrix(0L, nrow=nrow(x), ncol=1)
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
