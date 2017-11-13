context("ievora")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_ievora <- function(mat, groups, tCut=0.05, bCut=0.001) {
  doDV <- function(tmp.v,pheno.v) {
    co.idx <- which(pheno.v==0)
    ca.idx <- which(pheno.v==1)
    bt.o <- bartlett.test(x=tmp.v,g=pheno.v)
    pv <- bt.o$p.value
    logR <- log2(var(tmp.v[ca.idx])/var(tmp.v[co.idx]))
    avCA <- mean(tmp.v[ca.idx])
    avCO <- mean(tmp.v[co.idx])
    out.v <- c(logR,pv,avCA,avCO)
    names(out.v) <- c("log(V1/V0)","P(BT)","Av1","Av0")
    return(out.v)
  }

  doTT <- function(tmp.v,pheno.v){
    tt.o <- t.test(tmp.v ~ pheno.v)
    out.v <- c(-tt.o$stat,tt.o$p.val)
    names(out.v) <- c("t","P")
    return(out.v)
  }

  iEVORA <- function(data.m,pheno.v,thDV=0.001,thDM=0.05) {

    statDVC.m <- t(apply(data.m,1,doDV,pheno.v))
    qvDVC.v <- p.adjust(statDVC.m[,2], "fdr")
    dvc.idx <- which(qvDVC.v < thDV)
    nDVC <- length(dvc.idx)
    if( nDVC > 0 ){
      statDMC.m <- t(apply(data.m[dvc.idx,],1,doTT,pheno.v))
      tmp.s <- sort(statDMC.m[,2],decreasing=FALSE,index.return=TRUE)
      pvDMC.v <- tmp.s$x
      ntop <- length(which(pvDMC.v < thDM))
      if(ntop > 0){
        topDVMC.m <- cbind(statDMC.m[tmp.s$ix[1:ntop],],
                           statDVC.m[dvc.idx[tmp.s$ix[1:ntop]],c(3:4,1:2)],
                           qvDVC.v[dvc.idx[tmp.s$ix[1:ntop]]]
                           )
        colnames(topDVMC.m) <- c("t","P(TT)","Av1","Av0","log[V1/V0]","P(BT)","q(BT)")
        rownames(topDVMC.m) <- rownames(statDMC.m)[tmp.s$ix[1:ntop]]
      } else {
        topDVMC.m <- NULL
      }
    } else {
      topDVMC.m <- NULL
    }
    topDVMC.m
  }

  res <- as.data.frame(iEVORA(mat, groups, bCut, tCut))
  res$Av1 <- rowMeans(mat[rownames(res),groups==1], na.rm=TRUE)
  res$Av0 <- rowMeans(mat[rownames(res),groups==0], na.rm=TRUE)
  res$var1 <- rowVars(mat[rownames(res),groups==1], na.rm=TRUE)
  res$var0 <- rowVars(mat[rownames(res),groups==0], na.rm=TRUE)
  res$logR <- log2(res$var1/res$var0)
  res$obs1 <- apply(mat[rownames(res),groups==1], 1, function(x) sum(!is.na(x)))
  res$obs0 <- apply(mat[rownames(res),groups==0], 1, function(x) sum(!is.na(x)))
  res$significant <- TRUE
  res$rank <- c(1:nrow(res))
  res <- res[,c(4,3,9,8,12,11,10,1,2,6,7,13,14)]
  colnames(res) <- c("mean.0", "mean.1", "var.0", "var.1", "obs.0", "obs.1",
                     "logR", "t.statistic", "tt.p.value", "bt.p.value",
                     "bt.q.value", "significant", "rank"
                     )
  res
}

################################################################################
#################### TEST CONSISTENCY WITH ORIGINAL VERSION ####################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X1 <- matrix(rnorm(10000,0,0.5), ncol=50)
  X2 <- matrix(rnorm(10000,0,1), ncol=50)
  X <- cbind(X1, X2)
  X[sample(length(X),500)] <- NA
  rownames(X) <- 1:nrow(X)
  groups <- rep(c(0,1), each=50)

  t1 <- base_ievora(X, groups, 0.05, 0.001)
  t2 <- ievora(X, groups, 0.05, 0.001)
  t2 <- t2[t2$significant,]
  t2 <- t2[order(t2$rank),]

  expect_equal(t1, t2)
})

################################################################################
############################### TEST PARAMETERS ################################
################################################################################

###################################### X #######################################

test_that("x cannot be missing", {
  er <- 'argument "x" is missing, with no default'
  expect_error(ievora(), er)
})

test_that("invalid x values produce errors", {
  er <- "'x' must be a numeric matrix or a vector"
  expect_error(ievora(iris), er)
  expect_error(ievora(matrix(LETTERS[1:10], ncol=2)), er)
  expect_error(ievora(list(1:10)), er)
  expect_error(ievora(NA), er)
  expect_error(ievora(NULL), er)
  expect_error(ievora(matrix(nrow=0, ncol=5)), er)
  expect_error(ievora(matrix(nrow=5, ncol=0)), er)
})

test_that("data.frame values pass when they are all numeric", {
  expect_error(ievora(t(iris[,-5]), iris$Species=="setosa"), NA)
})

test_that("zero rows matrix pass", {
  expect_error(ievora(matrix(NA_integer_, nrow=0, ncol=2), c(0,1)), NA)
})

#################################### GROUPS ####################################

test_that("groups cannot be missing", {
  er <- 'argument "groups" is missing, with no default'
  expect_error(ievora(1:5), er)
})
