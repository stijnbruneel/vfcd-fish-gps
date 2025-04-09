'manovRDa_personal' <- function(Y, s, tt, S.mat=NULL, T.mat=NULL, STI.mat=NULL, Sfixed=TRUE, Tfixed=TRUE, S.test=TRUE, T.test=TRUE, STI.test=TRUE, model = "5", nperm=999)
  
{
  
  n <- nrow(Y)
  p <- ncol(Y)
  a <- ncol(S.mat)
  b <- ncol(T.mat)
  
  if(!is.null(STI.mat)) {cc <- ncol(STI.mat)}
  else {cc = 0}
  
  A <- S.mat
  B <- T.mat
  if(!is.null(STI.mat)) AxB <- STI.mat
  
  # Compute projector of A and Yfit.A
  invA <- ginv(t(A) %*% A)
  projA <- A %*% invA %*% t(A)
  Yfit.A <- projA %*% Y
  
  # Compute projector of B and Yfit.B
  invB <- ginv(t(B) %*% B)
  projB <- B %*% invB %*% t(B)
  Yfit.B <- projB %*% Y
  
  # Compute projector of AxB and Yfit.AxB
  if(!is.null(STI.mat)) {
    invAxB <- ginv(t(AxB) %*% AxB)
    projAxB <- AxB %*% invAxB %*% t(AxB)
    Yfit.AxB <- projAxB %*% Y
  } else {
    projAxB=NULL  
  }
  
  # Create a "compound matrix" to obtain R-square and adjusted R-square
  if(!is.null(STI.mat)) {
    ABAxB <- cbind(A,B,AxB)
  } else {
    ABAxB <- cbind(A,B)
  }
  
  # Compute projector of ABAxB and Yfit.ABAxB
  invABAxB <- ginv(t(ABAxB) %*% ABAxB)
  projABAxB <- ABAxB %*% invABAxB %*% t(ABAxB)
  Yfit.ABAxB <- projABAxB %*% Y
  
  # Compute Sums of squares (SS) and Mean squares (MS)
  SS.Y <- sum(Y*Y)
  SS.Yfit.ABAxB <- sum(Yfit.ABAxB*Yfit.ABAxB)
  SS.Yfit.A <- sum(Yfit.A*Yfit.A)
  SS.Yfit.B <- sum(Yfit.B*Yfit.B)
  if(!is.null(STI.mat)) SS.Yfit.AxB <- sum(Yfit.AxB*Yfit.AxB)
  MS.A <- SS.Yfit.A/a
  MS.B <- SS.Yfit.B/b
  if(!is.null(STI.mat)) MS.AxB <- SS.Yfit.AxB/cc
  MS.Res <- (SS.Y-SS.Yfit.ABAxB)/(n-(a+b+cc)-1)
  
  
  if(STI.test==TRUE) { # Test interaction (unrestricted permutations)
    nPGE.AxB = 1
    Fref.AxB <- MS.AxB/MS.Res
    
    nPGE.AxB=.Call("sti_loop",nperm,Y,s,tt,a,b,cc,SS.Y,Fref.AxB,projAxB,projABAxB)  ### NM
    P.AxB <- nPGE.AxB/(nperm+1)
    
    R2 <- SS.Yfit.AxB/SS.Y
    R2a <- 1-((n-1)/(n-dim(STI.mat)[2]-1))*(1-R2)
    
    testSTI <- list(MS.num=MS.AxB, MS.den=MS.Res, R2=R2, R2.adj=R2a, F=Fref.AxB, Prob=P.AxB)
  } else {
    testSTI <- NULL
  }
  
  
  if(S.test==TRUE) { # Test factor A (space) using restricted permutations within time blocks
    
    #########################
    nPGE.A=1
    if(Tfixed==FALSE && !is.null(STI.mat)) { # Time random factor in crossed design with interaction
      Fref.A <- MS.A/MS.AxB
      MS.den <- MS.AxB
      T_fixed <- 1
    } else if(Tfixed==FALSE && model=="6b") { # Time random factor in nested design
      Fref.A <- MS.A/MS.B
      MS.den <- MS.B
      T_fixed <- 2
    } else {
      Fref.A <- MS.A/MS.Res
      MS.den <-  MS.Res
      T_fixed <- 3
    }
    
    nPGE.A=.Call("s_loop",nperm,Y,s,tt,a,b,cc,SS.Y,Fref.A,projA,projB,projAxB,projABAxB,T_fixed)  ### NM
    
    P.A <- nPGE.A/(nperm+1)
    
    R2 <- SS.Yfit.A/SS.Y
    R2a <- 1-((n-1)/(n-a-1))*(1-R2)
    
    testS <- list(MS.num=MS.A, MS.den=MS.den, R2=R2, R2.adj=R2a, F=Fref.A, Prob=P.A)
  } else {
    testS <- NULL
  }
  
  if(T.test==TRUE) { # Test factor B (time) using restricted permutations within time blocks
    nPGE.B = 1
    if(Sfixed==FALSE && !is.null(STI.mat)) { # Space random factor in crossed design with interaction
      Fref.B <- MS.B/MS.AxB
      MS.den <- MS.AxB
      T_fixed=1   
    } else if(Sfixed==FALSE && model=="6a") { # Space random factor in nested design
      Fref.B <- MS.B/MS.A
      MS.den <- MS.A
      T_fixed=2    
    } else {
      Fref.B <- MS.B/MS.Res
      MS.den <- MS.Res
      T_fixed=3     
    }
    
    
    nPGE.B <- .Call("t_loop",nperm,Y,s,tt,a,b,cc,SS.Y,Fref.B,projA,projB,projAxB,projABAxB,T_fixed)  # NM
    
    P.B <- nPGE.B/(nperm+1)
    
    R2 <- SS.Yfit.B/SS.Y
    R2a <- 1-((n-1)/(n-b-1))*(1-R2)
    
    testT <- list(MS.num=MS.B, MS.den=MS.den, R2=R2, R2.adj=R2a, F=Fref.B, Prob=P.B)
  } else {
    testT <- NULL
  }
  
  return(list(testSTI = testSTI, testS = testS, testT=testT))
}