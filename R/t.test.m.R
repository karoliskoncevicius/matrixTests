
t.test.m.default <- function(x, y=NULL,
                             alternative=c("two.sided", "less", "greater"),
                             mu=0, paired=FALSE, var.equal=FALSE,
                             conf.level=0.95, ...
                             ) {

  alternative <- match.arg(alternative)

  if(is.vector(x))
    x <- matrix(x, nrow=1)
  if(is.vector(y))
    y <- matrix(y, nrow=1)

  if(!is.matrix(x))
    stop("'x' should be a matrix or a vector")

  if(!missing(mu) && ((length(mu) != 1 & length(mu) != nrow(x)) || is.na(mu)))
    stop("'mu' must be a single number or a vector of length nrow(x)")

  if(!missing(conf.level) &&
     (length(conf.level)!=1 & length(conf.level)!=nrow(x)))
    stop("'conf.level' must be a single number or a vector of length nrow(x)")

  if(any(!is.finite(conf.level) | conf.level<0 | conf.level>1))
    stop("'conf.level' must be between 0 and 1")

  if(!is.null(y)) {
    if(!is.matrix(y))
      stop("'y' should be a matrix or a vector")
    if(!all.equal(dim(x), dim(y)))
      stop("'x' and 'y' must have same dimensions")
    xnames <- paste0(deparse(substitute(x)), "[", seq_along(x), " ,]")
    ynames <- paste0(deparse(substitute(y)), "[", seq_along(y), " ,]")
    dnames <- paste(xnames, "and", ynames)
    if(paired) {
      x[is.na(y)] <- NA
      y[is.na(x)] <- NA
    }
  } else {
    dnames <- paste0(deparse(substitute(x)), "[", seq_along(x), " ,]")
    if(paired)
      stop("'y' is missing for paired test")
  }

  if(paired) {
    x <- x - y
    y <- NULL
  }

  nxs <- rowSums(!is.na(x))
  mxs <- rowMeans(x, na.rm=TRUE)
  vxs <- rowVars(x, na.rm=TRUE)
  messages <- rep("OK", length(mxs))
  if(is.null(y)) {
    bad           <- nxs<2
    mxs[bad]      <- NA
    messages[bad] <- "not enough 'x' observations"
    dfs     <- nxs-1
    stderrs <- sqrt(vxs/nxs)
    bad           <- stderrs < 10 * .Machine$double.eps * abs(mxs)
    mxs[bad]      <- NA
    messages[bad] <- "data are essentially constant"
    tstats <- (mxs - mu)/stderrs
    method <- ifelse(paired, "Paired t-test", "One Sample t-test")
    # estimate <- setNames(mxs, ifelse(paired, "mean of the differences", "mean of x"))
    estimates <- matrix(mxs, ncol=1)
    colnames(estimates) <- ifelse(paired, "Paired t-test", "One Sample t-test")
  } else {
    nys <- rowSums(!is.na(y))
    mys <- rowMeans(y, na.rm=TRUE)
    vys <- rowVars(y, na.rm=TRUE)
    bad           <- nxs<1 | (!var.equal & nxs<2)
    mxs[bad]      <- NA
    messages[bad] <- "not enough 'x' observations"
    bad           <- nys<1 | (!var.equal & nys<2)
    mys[bad]      <- NA
    messages[bad] <- "not enough 'y' observations"
    bad           <- var.equal & (nxs+nys)<3
    mxs[bad]      <- NA
    mys[bad]      <- NA
    messages[bad] <- "not enough observations"

    method    <- paste(ifelse(!var.equal, "Welch", ""), "Two Sample t-test")
    estimates <- cbind(mxs, mys)
    colnames(estimates) <- c("mean of x", "mean of y")
    if(var.equal) {
      dfs <- nxs + nys - 2
      vs <- rep(0, length(nxs))
      vs[nxs>1] <- vs + (nxs-1) * vxs
      vs[nys>1] <- vs + (nys-1) * vys
      vs <- vs/dfs
      stderrs <- sqrt(vs * (1/nxs + 1/nys))
    } else {
      stderrxs <- sqrt(vxs/nxs)
      stderrys <- sqrt(vys/nys)
      stderrs  <- sqrt(stderrxs^2 + stderrys^2)
      dfs <- stderrs^4/(stderrxs^4/(nxs - 1) + stderrys^4/(nys - 1))
    }
    bad <- stderrs < 10 * .Machine$double.eps * rowMaxs(abs(cbind(mxs, mys)))
    mxs[bad] <- NA
    mys[bad] <- NA
    messages[bad] <- "data are essentially constant"
    tstats <- (mxs - mys - mu)/stderrs
  }
  if(alternative=="less") {
    pvals <- pt(tstats, dfs)
    cints <- cbind(-Inf, tstats + qt(conf.level, dfs))
  } else if (alternative=="greater") {
    pval <- pt(tstats, dfs, lower.tail=FALSE)
    cints <- cbind(tstats - qt(conf.level, dfs), Inf)
  } else {
    pvals  <- 2 * pt(-abs(tstats), dfs)
    alpha <- 1 - conf.level
    cints <- qt(1 - alpha/2, dfs)
    cints <- tstats + cbind(-cints, cints)
  }
  cints <- mu + cints * stderrs

  names(tstats) <- "t"
  names(dfs) <- "df"
  names(mu) <- ifelse((paired || !is.null(y)), "difference in means", "mean")
  attr(cints, "conf.level") <- conf.level
  rval <- Map(list, statistic=as.list(tstats), parameter=as.list(dfs),
              p.value=as.list(pvals), conf.int=split(cints, row(cints)),
              estimate=split(estimates, row(estimates)), null.value=as.list(mu),
              alternative=as.list(alternative), method=as.list(method),
              data.name=as.list(dnames)
              )
  # rval <- list(statistic=tstats, parameter=dfs, p.value=pvals,
  #              conf.int=cints, estimate=estimates, null.value=mu,
  #              alternative=alternative, method=method, data.name=dnames
  #              )
  rval <- Map(`class<-`, rval, "htest")
  rval
}


t.test.m.formula <- function (formula, data, subset, na.action, ...) {
    if (missing(formula) || (length(formula) != 3L) || (length(attr(terms(formula[-2L]),
        "term.labels")) != 1L))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (nlevels(g) != 2L)
        stop("grouping factor must have exactly 2 levels")
    DATA <- setNames(split(mf[[response]], g), c("x", "y"))
    y <- do.call("t.test", c(DATA, list(...)))
    y$data.name <- DNAME
    if (length(y$estimate) == 2L)
        names(y$estimate) <- paste("mean in group", levels(g))
    y
}
