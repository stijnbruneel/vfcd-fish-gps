'stimodels_personal' <- function(Y, S, Ti,  model="5", nperm=999, nS=-1, nT=-1, Sfixed=TRUE, Tfixed=TRUE, COD.S=NULL, COD.T=NULL, print.res=TRUE, Y.type="community")
{
  if(!is.logical(print.res)) {
    stop("Wrong operator; 'print.res' should be either 'FALSE' or 'TRUE'", call.=FALSE)
  }
  if(model!="2" && model!="3a" && model!="3b" && model!="4" && model!="4" && model!="5" && model!="6a" && model!="6b" && model!="7") {
    stop(paste("Unrecognized model ",model,"; 'model' should be '2', '3a', '3b', '4', '5', '6a','6b' or '7'.", sep=""), call.=FALSE)
  }
  
  aa <- system.time( {
    # Sets the number and location of spatial and temporal points
    S <- as.matrix(S)
    if(dim(S)[1]==1 && dim(S)[2]==1) {
      s <- S[1,1]
      sitesX <- c(1:s)
    } else {
      s <- dim(S)[1]
      sitesX <- S
    }	
    weightsTi=aem.weight.time(Ti) 
    Ti <- as.matrix(Ti)
    if(dim(Ti)[1]==1 && dim(Ti)[2]==1) {
      tt <- Ti[1,1]
      timesX <- c(1:tt) 
    } else {
      tt <- dim(Ti)[1]
      timesX <- Ti
    }
    
    if (Y.type!="community"){
      Y.dist=vegdist(Y,method="bray",correction="cailliez")
      Y=pcoa(Y.dist)$vectors
    }
    
    # Total number of rows in matrix
    n <- s*tt			
    p <- dim(Y)[2]
    
    # Check response data file containing species data
    Y <- as.matrix(Y)
    p <- dim(Y)[2]
    
    if(dim(Y)[1] != n) stop("The number of rows in species file is not (S x Ti)", call.=FALSE) 
    
    # Center response data
    Y <- scale(Y, center=TRUE, scale=FALSE)
    
    
    if(print.res) {
      cat("=======================================================\n")
      cat("        Space-time ANOVA without replicates\n")
      cat("                                                  \n")
      cat("  Pierre Legendre, Miquel De Caceres, Daniel Borcard\n")
      cat("=======================================================\n\n")
      cat(" Number of space points (s) =", s,'\n')
      cat(" Number of time points (tt) =", tt,'\n')
      cat(" Number of observations (n = s*tt) =", n,'\n')
      cat(" Number of response variables (p) =", p,'\n','\n')
    }
    
    # Generates space dbMEM variables (if necessary)
    if(model=="3a"|| model=="6a" || model=="7"|| model=="4"|| model=="5") {
      if(is.null(COD.S)) {	# Generate spatial dbMEMs if not given by user
        if(print.res) cat(" Computing dbMEMs to code for space\n")
        if(s==2) {
          dbMEM.S <- as.matrix(rep(c(-1,1),tt))
          cat("\n There are only two space points. In this case, conputation of dbMEMs")
          cat("\n is useless. They are replaced by a vector of Helmert contrasts\n\n")
        } 
        else {
          #dbMEM.S.tmp <- dbmem(sitesX, MEM.autocor="positive")
    
          thresh<-give.thresh(dist(sitesX))
          nbnear <- dnearneigh(sitesX, 0, thresh)
          distgab <- nbdists(nbnear, sitesX)
          fdist <- lapply(distgab, function(x) 1 - x/max(dist(sitesX)))
          listwgab <- nb2listw(nbnear, glist = fdist, style = "B")
          dbMEM.S.tmp<-mem(listwgab,MEM.autocor="positive")
          
          SS <- as.matrix(dbMEM.S.tmp)
          dbMEM.S.thresh <- give.thresh(dist(sitesX))
          if(print.res) cat(" Truncation level for space dbMEMs =", dbMEM.S.thresh, "\n")
          dbMEM.S <- SS
          for(j in 2:tt) dbMEM.S <- rbind(dbMEM.S,SS)
          if(nS==-1) nS <- ncol(SS)
          else {
            if(nS > ncol(SS)) {
              cat("\n Number of requested spatial variables nS =", nS, "larger than available\n")
              cat(" The", ncol(SS), "available spatial dbMEM will be used\n\n")
              nS <- ncol(SS)
            }
            else dbMEM.S <- dbMEM.S[, 1:nS, drop = FALSE]
          }
        }
      } else {
        dbMEM.S <- apply(as.matrix(COD.S),2,scale,center=TRUE,scale=TRUE)
        nS <- dim(dbMEM.S)[2]
        if(nS>=s) stop("The number of spatial coding functions must be lower than s", call.=FALSE) 
        if(nrow(as.matrix(COD.S))!=dim(Y)[1]) stop("The number of rows in COD.S must be equal to the number of observations", call.=FALSE) 
      }
    }
    
    # Generate dbMEM variables for time (if necessary)
    if(model=="3b" || model=="6b"|| model=="7"||model=="4" || model=="5") {
      if(is.null(COD.T)) {      # Generate temporal dbMEM variables if not given by user
        if(tt==2) {  
          dbMEM.T=as.matrix(c(rep(-1, s), rep(1,s)))  #DB
          cat("\n There are only two time points. In this case, conputation of dbMEMs")
          cat("\n is useless. They are replaced by a vector of Helmert contrasts\n\n")
        } else {
          if(print.res) cat(" Computing dbMEMs to code for time\n")
          dbMEM.T.tmp <- aem.time(nrow(Ti),weightsTi,moran=TRUE)
          dbMEM.T.tmp <- dbMEM.T.tmp$aem[,which(dbMEM.T.tmp$Moran$obs>0)]
          TT <- as.matrix(dbMEM.T.tmp)
          dbMEM.T.thresh <- give.thresh(dist(timesX))
          if(print.res) cat(" Truncation level for time dbMEMs =", dbMEM.T.thresh, "\n\n")
          T.temp <- TT[1,]
          for(i in 2:s) T.temp <- rbind(T.temp,TT[1,])
          dbMEM.T <- as.matrix(T.temp)
          for(j in 2:tt) {
            T.temp <- TT[j,]
            for(i in 2:s) T.temp <- rbind(T.temp,TT[j,])
            dbMEM.T <- as.matrix(rbind(dbMEM.T,T.temp))
          }
          if(nT==-1)  nT <- ncol(TT)
          else {
            if(nT > ncol(TT)) {
              cat("\n Number of requested temporal variables nT =", nT, "larger than available\n")
              cat(" The", ncol(TT), "available temporal dbMEM will be used\n\n")
              nT <- ncol(TT)
            }
            else dbMEM.T <- dbMEM.T[, 1:nT, drop = FALSE]
          }
          if(nT) {
            if(nT > ncol(TT)) {
              cat("Number of requested temporal variables nT=", nT, "larger than available\n")
              cat("The", ncol(TT), "available temporal dbMEM will be used\n\n")
              nT <- ncol(TT)
            }
            else dbMEM.T <- dbMEM.T[, 1:nT, drop = FALSE]
          }
          else  nT <- ncol(TT)
        }
        
        
      } else {
        dbMEM.T=apply(as.matrix(COD.T),2,scale,center=TRUE,scale=TRUE)
        nT <- dim(dbMEM.T)[2]
        if(nT>=tt) stop("The number of temporal coding functions must be lower than t", call.=FALSE)
        if(nrow(as.matrix(COD.T))!=dim(Y)[1]) stop("The number of rows in COD.T must be equal to the number of observations", call.=FALSE) 
      }
    }
    
    if(print.res) {
      if(model!="2" && model!="3b" && model!="6b") { 
        if(s==2)   #DB To avoid number = -1
          cat("Number of space coding functions = 1 \n")
        else
          cat(" Number of space coding functions =", nS,'\n\n')
      }
      if(model!="2" && model!="3a" && model!="6a") {
        if(tt==2)  #DB To avoid number = -1
          cat("Number of time coding functions = 1 \n")
        else
          cat(" Number of time coding functions =", nT, "\n\n")}            		
    }
    
    
    # Generates space and time helmert contrasts
    A <- as.factor(rep(1:s,tt))
    B <- rep(1,s)
    for(i in 2:tt) B <- c(B,rep(i,s))
    B <- as.factor(B)
    HM <- model.matrix(~ A + B, contrasts = list(A="contr.helmert", B="contr.helmert"))
    HM.S <- as.matrix(HM[,2:s])
    HM.T <- as.matrix(HM[,(s+1):(s+tt-1)])
    
    test.STI = TRUE
    test.S = TRUE
    test.T = TRUE
    
    # Defines X (variables for the factor of interest) and W (covariables) for each test to be done
    if(model=="5") {
      XSTI <- dbMEM.S*dbMEM.T[,1]
      if(dim(dbMEM.T)[2]>1) for(j in 2:dim(dbMEM.T)[2]) XSTI <- cbind(XSTI,dbMEM.S*dbMEM.T[,j])
      XS <- HM.S
      XT <- HM.T
      if(print.res) {
        cat(" MODEL V: HELMERT CONTRAST FOR TESTING MAIN FACTORS. \n")
        cat("          SPACE AND TIME dbMEMs FOR TESTING INTERACTION.",'\n')
      }					
    } else if(model=="4") {
      XSTI <- dbMEM.S*dbMEM.T[,1]
      if(dim(dbMEM.T)[2]>1) for(j in 2:dim(dbMEM.T)[2]) XSTI = cbind(XSTI,dbMEM.S*dbMEM.T[,j])
      XS <- dbMEM.S
      XT <- dbMEM.T
      if(print.res) {
        cat(" MODEL IV: dbMEMs FOR BOTH SPACE AND TIME.",'\n')
      }		
    } else if(model=="3a") {
      XSTI <- dbMEM.S*HM.T[,1]
      if(tt>2) for(j in 2:(tt-1)) XSTI <- cbind(XSTI,dbMEM.S*HM.T[,j])
      XS <- dbMEM.S
      XT <- HM.T
      if(print.res) {
        cat(" MODEL IIIa: dbMEMs FOR SPACE AND HELMERT CONTRASTS FOR TIME.",'\n')
      }	
    } else if(model=="3b") {
      XSTI <- dbMEM.T*HM.S[,1]
      if(dim(HM.S)[2]>1)
        for(j in 2:(s-1)) XSTI <- cbind(XSTI,dbMEM.T*HM.S[,j])
        XS <- HM.S
        XT <- dbMEM.T
        if(print.res) {
          cat(" MODEL IIIb: HELMERT CONTRASTS FOR SPACE AND dbMEMs FOR TIME.",'\n')
        }		
    } else if(model=="7") {
      XSTI <- HM.S*HM.T[,1]
      if(dim(HM.T)[2]>1) for(j in 2:dim(HM.T)[2]) XSTI <- cbind(XSTI,HM.S*HM.T[,j])
      XS <- dbMEM.S
      XT <- dbMEM.T
      if(print.res) {
        cat(" MODEL VII: dbMEMs FOR BOTH SPACE AND TIME BUT HELMERT CONTRAST FOR INTERACTION.",'\n')
      }		
    } else if(model=="6a") {
      XS <- dbMEM.S
      for(j in 1:(tt-1)) XS <- cbind(XS,dbMEM.S*HM.T[,j])
      XT <- HM.T
      XSTI = NULL
      if(print.res) {
        cat(" MODEL VIa: NESTED MODEL.\n")
        cat("            TESTING FOR THE EXISTENCE OF SPATIAL STRUCTURE (COMMON OR SEPARATE)",'\n')
      }	
      test.STI = FALSE
    } else if(model=="6b") {
      XT <- dbMEM.T
      for(j in 1:(s-1)) XT <- cbind(XT,dbMEM.T*HM.S[,j])
      XS <- HM.S						
      XSTI = NULL
      if(print.res) {
        cat(" MODEL VIb: NESTED MODEL.\n")
        cat("            TESTING FOR THE EXISTENCE OF TEMPORAL STRUCTURE (COMMON OR SEPARATE).",'\n')
      }	
      test.STI = FALSE
    } else if(model=="2") {
      XS <- HM.S
      XT <- HM.T
      XSTI = NULL
      if(print.res) {
        cat(" MODEL II: HELMERT CONTRAST FOR SPACE AND TIME. NO INTERACTION TERM.",'\n')
      }	
      test.STI = FALSE
    } 		
    if(print.res) {
      cat("   Number of space variables =", dim(XS)[2],'\n')
      cat("   Number of time variables =", dim(XT)[2],'\n')
      if(test.STI) {
        cat("   Number of interaction variables =", dim(XSTI)[2],'\n')
        cat("   Number of residual degrees of freedom =", (s*tt-dim(XS)[2]-dim(XT)[2]-dim(XSTI)[2]-1),"\n\n")
        if((s*tt-dim(XS)[2]-dim(XT)[2]-dim(XSTI)[2]-1) <= 0) stop("Not enough residual degrees of freedom for testing. \nTry with a lower number of space or time variables.", call.=FALSE)
      } else {
        cat("   Number of residual degrees of freedom =", (s*tt-dim(XS)[2]-dim(XT)[2]-1),"\n\n")			
        if((s*tt-dim(XS)[2]-dim(XT)[2]-1) <= 0) stop("Not enough residual degrees of freedom for testing. \nTry with a lower number of space or time variables.", call.=FALSE)
        
      }
    }		
    
    
    res <- manovRDa_personal(Y=Y,s=s,tt=tt,S.mat=XS,T.mat=XT,STI.mat=XSTI, Sfixed= Sfixed, Tfixed=Tfixed, S.test=test.S, T.test=test.T, STI.test=test.STI, model = model, nperm=nperm)
    
    if(test.STI==TRUE) {
      cat(' Interaction test:   R2 =', round(res$testSTI$R2, 4),'  F =',round(res$testSTI$F, 4),'  P(',nperm,'perm) =',res$testSTI$Prob,'\n')
    }
    if(test.S==TRUE) {
      cat(' Space test:         R2 =', round(res$testS$R2, 4),'  F =',round(res$testS$F, 4),'  P(',nperm,'perm) =',res$testS$Prob,'\n')
    }
    if(test.T==TRUE) {
      cat(' Time test:          R2 =', round(res$testT$R2, 4),'  F =',round(res$testT$F, 4),'  P(',nperm,'perm) =',res$testT$Prob,'\n')
    }
    
  })
  aa[3] <- sprintf("%2f",aa[3])
  if(print.res) {
    cat("-------------------------------------------------------\n")
    cat("      Time for this computation =",aa[3]," sec",'\n')
    cat("=======================================================\n\n")
  }
  class(res) <- "sti"
  invisible(res)
}