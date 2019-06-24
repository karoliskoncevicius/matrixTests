rowVars <- function(x, na.rm=FALSE) {
  M <- rowMeans(x, na.rm=na.rm)
  n <- rowSums(!is.na(x))-1
  ifelse(n > 0, rowSums((x-M)^2, na.rm=na.rm) / n, NA)
}

rowTies <- function(x) {
  dups <- t(apply(x, 1, duplicated, incomparables=NA))
  dups <- cbind(which(dups, arr.ind=TRUE), val=x[dups])
  inds <- !duplicated(dups[,c(1,3),drop=FALSE])
  dups[,3] <- stats::ave(dups[,1], dups[,1], dups[,3], FUN=length)+1
  dups <- dups[inds,,drop=FALSE]
  dups[,2] <- stats::ave(dups[,2], dups[,1], FUN=seq_along)

  res <- matrix(0L, nrow=nrow(x), ncol=max(1, dups[,2]))
  res[dups[,1:2,drop=FALSE]] <- dups[,3]
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
