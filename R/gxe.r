genexE.association.test <- function(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), E,
                             method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                             test = c("score", "wald", "lrt"), df=c(1,2,3), 
                             K, eigenK, beg = 1, end = ncol(x), p = 0, 
                             tol = .Machine$double.eps^0.25, ...) {

  if(beg < 1 || end > ncol(x)) stop("range too wide")
  if(is.null(x@mu)) stop("Need mu to be set in x")
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  if(length(E) != nrow(x)) stop("Dimensions of E and x mismatch")
  
  X <- as.matrix(X)
  if(nrow(X) != nrow(x)) stop("Dimensions of Y and x mismatch")

  response <- match.arg(response)
  test <- match.arg(test)
  
  # preparation de X 
  if(p > 0) {
    if((method == "lmm" & response == "quantitative" & test == "score") | 
       (method == "lmm" & response == "binary") |
       (method == "lm")) { # il faut ajouter les PCs à X
      X <- cbind(X, eigenK$vectors[,seq_len(p)])
      X <- trans.X(X, mean.y = mean(Y))
    } else { 
      X <- trans.X(X, eigenK$vectors[,seq_len(p)], mean(Y))
    }
  } else {
    X <- trans.X(X, mean.y = mean(Y))
  }
  
  # check dimensions before anything
  n <- nrow(x)
  if(!missing(K)) {
    if(n != nrow(K) | n != ncol(K)) stop("K and x dimensions don't match")
  }
  if(!missing(eigenK)) {
    if(n != nrow(eigenK$vectors) | n != ncol(eigenK$vectors) | n != length(eigenK$values)) 
      stop("eigenK and x dimensions don't match")
  }

  # random effect
  if(match.arg(method) == "lmm") { 
     # if(response == "binary" & test == "score") {
      # warning('Binary phenotype and method = "lmm" force test = "score"')
      # test <- "score"
     # }

    if(test == "score" | response == "binary") {
      if(missing(K)) stop("For a score test and for binary traits, argument K is mandatory")
      # avec le score test on peut gérer les données manquantes dans Y
      if( any(is.na(Y)) ) {
        w <- !is.na(Y)
        X <- as.matrix(X[w,])
        Y <- Y[w]
        if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
        warning(sum(!w), 'individuals with missing phenotype are ignored.\n')
      } 
    } else {
      if(missing(eigenK)) 
        stop("For quantitative Wald and LRT tests, argument eigenK is mandatory")
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
      if( any(is.na(E)) )
        stop("Can't handle missing data in E, please recompute eigenK for the individuals with non-missing phenotype")
    }

    if(response == "quantitative") { # score (argument K), wald ou lrt (eigen K) possibles
      if(test == "score") {
        if (df==1) {
          stop("For 1 df interaction test, score is not judicious (computationnaly heavy). Try LRT or Wald tesst.")
        } else if (df %in% 2:3) {
          if (df==2) model <- lmm.aireml(Y, X = cbind(X, E), K, get.P = TRUE, ... )
          if (df==3) model <- lmm.aireml(Y, X = X, K, get.P = TRUE, ... )
          t <- .Call("gg_GxE_lmm_score_bed", PACKAGE = "GEnX", x@bed, x@p, model$Py, model$P, E, df, beg-1, end-1)
          t$p <- pchisq( t$score, df = df, lower.tail=FALSE)
        } else stop("df must be equal to 1, 2, or 3.")		
      } else if(test == "wald") {
        if (df %in% 1:3) {
          X <- cbind(X, E, 0, 0) # space for the SNP, E and SNPxE
          t <- .Call("gg_GxE_lmm_wald_bed", PACKAGE = "GEnX", x@bed, x@p, Y, X, p, eigenK$values, eigenK$vectors, df, beg-1, end-1, tol)
          t$p <- pchisq( t$Wald, df = df, lower.tail=FALSE)
        } else stop("df must be equal to 1, 2, or 3.")
      } else { # test == "lrt"
        if (df %in% 1:3) {
          X <- cbind(X, E, 0, 0) # space for the SNP and interaction
          t <- .Call("gg_GxE_lmm_lrt_bed", PACKAGE = "GEnX", x@bed, x@p, Y, X, p, eigenK$values, eigenK$vectors, df, beg-1, end-1, tol)
          t$p <- pchisq( t$LRT, df = df, lower.tail=FALSE)
          } else stop("df must be equal to 1, 2, or 3.")
      }
    } else { # response == "binary", seulement le score test, avec argument K
      if(test == "score") {
        if (df==1) {
          stop("For 1 df interaction test, score is not judicious (computationnaly heavy). Try 2 ou 3 df test or Wald test.")
        } else if (df %in% 2:3) {
          if (df==2) {
            X <- cbind(X, E)
            model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
          }
          if (df==3) model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
          omega <- model$BLUP_omega + X%*%model$BLUP_beta
          pi <- 1/(1+exp(-omega))
          t <- .Call("gg_GxE_lmm_score_bed", PACKAGE = "GEnX", x@bed, x@p, Y-pi, model$P, E, df, beg-1, end-1)
          t$p <- pchisq( t$score, df = df, lower.tail=FALSE)
        } else stop("df must be equal to 1, 2, or 3.")		
      } else if(test == "wald") {
        #stop("Wald tests for binary trait not available")
        X <- cbind(X, E, 0, 0) # E and space for the SNP and SNPxE
        if (df %in% 1:3) {
          t <- .Call("gg_GxE_logitmm_wald_bed", PACKAGE = "GEnX", x@bed, x@p, Y, X, K, df, beg-1, end-1, tol)
          t$p <- pchisq( t$Wald, df = df, lower.tail=FALSE)
        } else stop("df must be equal to 1, 2, or 3.")
      } else stop("LRT tests for binary trait not available")
    }
  }

  # only fixed effects
  if(match.arg(method) == "lm") {
    stop("'lm' methods are not implemented yet.")
    if(test != "wald") warning('Method = "lm" force test = "wald"')
    if(response == "quantitative") {
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
      X <- cbind(X, E, 0, 0)
      if (df==1) {
        t <- .Call("gg_GxE_lm_quanti_1df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1);
        t$Wald <- (t$beta_ExSNP/t$sd_ExSNP)**2
        t$p <- pchisq( t$Wald, df = 1, lower.tail=FALSE)
      } else if (df==2) {
        t <- .Call("gg_GxE_lm_quanti_2df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1);
        t$p <- pchisq( t$Wald, df = 2, lower.tail=FALSE)
      } else if (df==3) {
        t <- .Call("gg_GxE_lm_quanti_3df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1);
        t$p <- pchisq( t$Wald, df = 3, lower.tail=FALSE)
      }
      #if( any(is.na(Y)) ) stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
      #X <- cbind(X, E, 0, 0)
      #if (df %in% 1:3) {
      #  t <- .Call("gg_GxE_lm_quanti_bed", PACKAGE = "GEnX", x@bed, x@p, Y, X, df, beg-1, end-1);
      #  t$p <- pchisq( t$Wald, df = df, lower.tail=FALSE)
      #} else stop("df must be equal to 1, 2, or 3.")
    }
    if(response == "binary") {
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
      X <- cbind(X, E, 0, 0)
      if (df==1) {
        t <- .Call("gg_GxE_logit_wald_1df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1, tol);
        t$Wald <- (t$beta_ExSNP/t$sd_ExSNP)**2
        t$p <- pchisq( t$Wald, df = 1, lower.tail=FALSE)
      } else if (df==2) {
        t <- .Call("gg_GxE_logit_wald_2df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1, tol);
        t$p <- pchisq( t$Wald, df = 2, lower.tail=FALSE)
      } else if (df==3) {
        t <- .Call("gg_GxE_logit_wald_3df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1, tol);
        t$p <- pchisq( t$Wald, df = 3, lower.tail=FALSE)
      }
      #if( any(is.na(Y)) ) stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
      #X <- cbind(X, E, 0, 0)
      #if (df %in% 1:3) {
      #  t <- .Call("gg_GxE_logit_wald_bed", PACKAGE = "GEnX", x@bed, x@p, Y, X, df, beg-1, end-1, tol);
      #  t$p <- pchisq( t$Wald, df = df, lower.tail=FALSE)
      #} else stop("df must be equal to 1, 2, or 3.")
    }
  }
  L <- list(chr = x@snps$chr, pos = x@snps$pos, id  = x@snps$id)
  if(beg > 1 | end < ncol(x))  # avoid copy
  L <- lapply(L, function(l) l[beg:end])

  data.frame( c( L, t) )
}



# Check AND replaces by QR decomposition...
trans.X <- function(X, PCs = matrix(0, nrow=nrow(X), ncol=0), mean.y = 1) {
  if(any(is.na(X)))
    stop("Covariates can't be NA")

  PCs <- as.matrix(PCs) # cas où p = 1
  n.X  <- ncol(X)
  n.pc <- ncol(PCs)
  n <- n.X + n.pc

  qr.X <- qr( cbind(PCs, X) );
  if(qr.X$rank < n) {
    warning("Covariate matrix X is not full rank, removing col(s) ", paste(qr.X$pivot[ seq(qr.X$rank+1,n) ] - n.pc , collapse = ", "))
    X <- X[ , qr.X$pivot[seq(n.pc+1, qr.X$rank)] - n.pc]
    qr.X <- qr(X)
  }
  if(mean.y > 1e-4) {
    X1 <- cbind(1,X);
    qr.X1 <- qr(X1);
    if(qr.X1$rank == ncol(X1)) {
      warning("An intercept column was added to the covariate matrix X")
      X <- X1;
      qr.X <- qr.X1
    }
  }
  if( qr.X$rank == ncol(X) )
    qr.Q(qr.X)
  else
    qr.Q(qr(X))
}


genexE.association.test_old <- function(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), E,
                             method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                             test = c("score", "wald", "lrt"), df=c(1,2,3), 
                             K, eigenK, beg = 1, end = ncol(x), p = 0, 
                             tol = .Machine$double.eps^0.25, ...) {

  if(beg < 1 || end > ncol(x)) stop("range too wide")
  if(is.null(x@mu)) stop("Need mu to be set in x")
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  if(length(E) != nrow(x)) stop("Dimensions of E and x mismatch")
  
  X <- as.matrix(X)
  if(nrow(X) != nrow(x)) stop("Dimensions of Y and x mismatch")

  response <- match.arg(response)
  test <- match.arg(test)
  
  # preparation de X 
  if(p > 0) {
    if((method == "lmm" & response == "quantitative" & test == "score") | 
       (method == "lmm" & response == "binary") |
       (method == "lm")) { # il faut ajouter les PCs à X
      X <- cbind(X, eigenK$vectors[,seq_len(p)])
      X <- trans.X(X, mean.y = mean(Y))
    } else { 
      X <- trans.X(X, eigenK$vectors[,seq_len(p)], mean(Y))
    }
  } else {
    X <- trans.X(X, mean.y = mean(Y))
  }
  
  # check dimensions before anything
  n <- nrow(x)
  if(!missing(K)) {
    if(n != nrow(K) | n != ncol(K)) stop("K and x dimensions don't match")
  }
  if(!missing(eigenK)) {
    if(n != nrow(eigenK$vectors) | n != ncol(eigenK$vectors) | n != length(eigenK$values)) 
      stop("eigenK and x dimensions don't match")
  }

  # random effect
  if(match.arg(method) == "lmm") { 
     # if(response == "binary" & test == "score") {
      # warning('Binary phenotype and method = "lmm" force test = "score"')
      # test <- "score"
     # }

    if(test == "score" | response == "binary") {
      if(missing(K)) stop("For a score test and for binary traits, argument K is mandatory")
      # avec le score test on peut gérer les données manquantes dans Y
      if( any(is.na(Y)) ) {
        w <- !is.na(Y)
        X <- as.matrix(X[w,])
        Y <- Y[w]
        if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
        warning(sum(!w), 'individuals with missing phenotype are ignored.\n')
      } 
    } else {
      if(missing(eigenK)) 
        stop("For quantitative Wald and LRT tests, argument eigenK is mandatory")
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
      if( any(is.na(E)) )
        stop("Can't handle missing data in E, please recompute eigenK for the individuals with non-missing phenotype")
    }

    if(response == "quantitative") { # score (argument K), wald ou lrt (eigen K) possibles
      if(test == "score") {
        if (df==1) {
          stop("For 1 df interaction test, score is not judicious (computationnaly heavy). Try LRT or Wald tesst.")
        } else if (df==2) {
          X <- cbind(X, E)
          model <- lmm.aireml(Y, X = X, K, get.P = TRUE, ... )
          t <- .Call("gg_GxE_lmm_score_2df", PACKAGE = "GEnX", x@bed, model$Py, model$P, x@mu, E, beg-1, end-1)
          t$p <- pchisq( t$score, df = 2, lower.tail=FALSE)
        } else if (df==3) {
          model <- lmm.aireml(Y, X = X, K, get.P = TRUE, ... )
          t <- .Call("gg_GxE_lmm_score_3df", PACKAGE = "GEnX", x@bed, model$Py, model$P, x@mu, E, beg-1, end-1)
          t$p <- pchisq( t$score, df = 3, lower.tail=FALSE)
        } else stop("df must be equal to 1, 2, or 3.")		
      } else if(test == "wald") {
        X <- cbind(X, E, 0, 0) # space for the SNP, E and SNPxE
        if (df==1) {
          t <- .Call("gg_GxE_lmm_wald_1df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
          t$Wald <- (t$beta_ExSNP/t$sd_ExSNP)**2
          t$p <- pchisq( t$Wald, df = 1, lower.tail=FALSE)
        } else if (df==2) {
          t <- .Call("gg_GxE_lmm_wald_2df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
          t$p <- pchisq( t$Wald, df = 2, lower.tail=FALSE)
        } else if (df==3) {
          t <- .Call("gg_GxE_lmm_wald_3df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
          t$p <- pchisq( t$Wald, df = 3, lower.tail=FALSE)
        }
      } else { # test == "lrt"
        X <- cbind(X, E, 0, 0) # space for the SNP and interaction
        if (df==1) {
          t <- .Call("gg_GxE_lmm_lrt_1df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
          t$p <- pchisq( t$LRT, df = 1, lower.tail=FALSE)
        } else if (df==2) {
          t <- .Call("gg_GxE_lmm_lrt_2df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
          t$p <- pchisq( t$LRT, df = 2, lower.tail=FALSE)
        } else if (df==3) {
          t <- .Call("gg_GxE_lmm_lrt_3df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, p, eigenK$values, eigenK$vectors, beg-1, end-1, tol)
          t$p <- pchisq( t$LRT, df = 3, lower.tail=FALSE)
        } else stop("df must be equal to 1, 2, or 3.")
      }
    } else { # response == "binary", seulement le score test, avec argument K
      if(test == "score") {
        if (df==1) {
          stop("For 1 df interaction test, score is not judicious (computationnaly heavy). Try 2 ou 3 df tests.")
        } else if (df==2) {
          X <- cbind(X, E)
          model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
          omega <- model$BLUP_omega
          if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
          pi <- 1/(1+exp(-omega))
          t <- .Call("gg_GxE_lmm_score_2df", PACKAGE = "GEnX", x@bed, Y-pi, model$P, x@mu, E, beg-1, end-1)
          t$p <- pchisq( t$score, df = 2, lower.tail=FALSE)
        } else if (df==3) {
          model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
          omega <- model$BLUP_omega
          if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
          pi <- 1/(1+exp(-omega))
          t <- .Call("gg_GxE_lmm_score_3df", PACKAGE = "GEnX", x@bed, Y-pi, model$P, x@mu, E, beg-1, end-1)
          t$p <- pchisq( t$score, df = 3, lower.tail=FALSE)
        } else stop("df must be equal to 1, 2, or 3.")		
      } else if(test == "wald") {
        #stop("Wald tests for binary trait not available")
        X <- cbind(X, E, 0, 0) # E and space for the SNP and SNPxE
        if (df==1) {
          t <- .Call("gg_GxE_logitmm_wald_1df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, K, beg-1, end-1, tol)
          t$Wald <- (t$beta_ExSNP/t$sd_ExSNP)**2
          t$p <- pchisq( t$Wald, df = 1, lower.tail=FALSE)
        } else if (df==2) {
          t <- .Call("gg_GxE_logitmm_wald_2df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, K, beg-1, end-1, tol)
          t$p <- pchisq( t$Wald, df = 2, lower.tail=FALSE)
        } else if (df==3) {
          t <- .Call("gg_GxE_logitmm_wald_3df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, K, beg-1, end-1, tol)
          t$p <- pchisq( t$Wald, df = 3, lower.tail=FALSE)
        } else stop("df must be equal to 1, 2, or 3.")
      } else stop("LRT tests for binary trait not available")
    }
  }

  # only fixed effects
  if(match.arg(method) == "lm") {
    if(test != "wald") warning('Method = "lm" force test = "wald"')
    if(response == "quantitative") {
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
      X <- cbind(X, E, 0, 0)
      if (df==1) {
        t <- .Call("gg_GxE_lm_quanti_1df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1);
        t$Wald <- (t$beta_ExSNP/t$sd_ExSNP)**2
        t$p <- pchisq( t$Wald, df = 1, lower.tail=FALSE)
      } else if (df==2) {
        t <- .Call("gg_GxE_lm_quanti_2df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1);
        t$p <- pchisq( t$Wald, df = 2, lower.tail=FALSE)
      } else if (df==3) {
        t <- .Call("gg_GxE_lm_quanti_3df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1);
        t$p <- pchisq( t$Wald, df = 3, lower.tail=FALSE)
      }
    }
    if(response == "binary") {
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
      X <- cbind(X, E, 0, 0)
      if (df==1) {
        t <- .Call("gg_GxE_logit_wald_1df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1, tol);
        t$Wald <- (t$beta_ExSNP/t$sd_ExSNP)**2
        t$p <- pchisq( t$Wald, df = 1, lower.tail=FALSE)
      } else if (df==2) {
        t <- .Call("gg_GxE_logit_wald_2df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1, tol);
        t$p <- pchisq( t$Wald, df = 2, lower.tail=FALSE)
      } else if (df==3) {
        t <- .Call("gg_GxE_logit_wald_3df", PACKAGE = "GEnX", x@bed, x@mu, Y, X, beg-1, end-1, tol);
        t$p <- pchisq( t$Wald, df = 3, lower.tail=FALSE)
      }
    }
  }
  L <- list(chr = x@snps$chr, pos = x@snps$pos, id  = x@snps$id)
  if(beg > 1 | end < ncol(x))  # avoid copy
  L <- lapply(L, function(l) l[beg:end])

  data.frame( c( L, t) )
}



# Check AND replaces by QR decomposition...
trans.X <- function(X, PCs = matrix(0, nrow=nrow(X), ncol=0), mean.y = 1) {
  if(any(is.na(X)))
    stop("Covariates can't be NA")

  PCs <- as.matrix(PCs) # cas où p = 1
  n.X  <- ncol(X)
  n.pc <- ncol(PCs)
  n <- n.X + n.pc

  qr.X <- qr( cbind(PCs, X) );
  if(qr.X$rank < n) {
    warning("Covariate matrix X is not full rank, removing col(s) ", paste(qr.X$pivot[ seq(qr.X$rank+1,n) ] - n.pc , collapse = ", "))
    X <- X[ , qr.X$pivot[seq(n.pc+1, qr.X$rank)] - n.pc]
    qr.X <- qr(X)
  }
  if(mean.y > 1e-4) {
    X1 <- cbind(1,X);
    qr.X1 <- qr(X1);
    if(qr.X1$rank == ncol(X1)) {
      warning("An intercept column was added to the covariate matrix X")
      X <- X1;
      qr.X <- qr.X1
    }
  }
  if( qr.X$rank == ncol(X) )
    qr.Q(qr.X)
  else
    qr.Q(qr(X))
}
