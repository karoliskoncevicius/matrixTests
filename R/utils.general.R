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

rowRankTies <- function(r) {
  if(storage.mode(r) != "integer") {
    storage.mode(r) <- "integer"
  }
  r <- t(r)
  res <- vector("list", ncol(r))
  for(i in seq_len(ncol(r))) {
    res[[i]] <- tabulate(r[,i], nrow(r))
  }
  if(length(res) == 0) {  # needed to work with zero-row inputs
    res <- list(integer())
  }
  do.call(rbind, res)
}

showWarning <- function(w, fname, msg) {
  if(any(w, na.rm=TRUE)) {
    callstack <- paste(deparse(sys.calls()), collapse="")
    prefix    <- ifelse(grepl(paste0("col_", fname), callstack), "col", "row")
    fname     <- paste0(prefix, "_", fname)
    prefix    <- ifelse(prefix == "col", "column", "row")
    n <- sum(w, na.rm=TRUE)
    i <- match(TRUE, w)
    msg <- paste0(fname, ": ", n, ' of the ', prefix, 's ', msg, ".",
                  '\nFirst occurrence at ', prefix, ' ', i
                  )
    warning(msg, call.=FALSE)
  }
}
