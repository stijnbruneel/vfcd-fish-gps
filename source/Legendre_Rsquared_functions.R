# Three R functions (in alphabetic order), written for the analysis of the 
# case study in the Practicals on Temporal eigenfunction methods.
# The whole file of functions can be loaded to the R console,  
# using the "Source R Code..." menu for Windows clients 
# or the "Source File..." menu for MacOS X clients,
# or the command source(file.choose()) written in the R console.

# ========

drop.levels <- function(data) 
{
  for (i in 1:ncol(data)) {
    if (is.factor(data[,i])) { data[,i] <- data[,i][, drop=TRUE] } 
  }
  return(data)
}

# ========

R2.by.variable <- function(Y, X, scale.Y=FALSE, nperm=99, random.var=FALSE)
  # Partial R-squares for scalogram -- 
  # Compute R-square and test of significance for each explanatory variable or MEM 
  # eigenfunction in the presence of all others.
  #
  # Arguments --
  # Y : Matrix or data frame of response data.
  # X : Matrix or data frame of explanatory variables (random variables or MEM).
  # scale.Y : =TRUE : standardize the Y variables by columns;
  #           =FALSE : do not standardize the Y variables.
  # nperm : Number of permutations for the tests of the semipartial R2.
  # random.var : =TRUE when X contains random explanatory variables;
#              =FALSE when X contains fixed explanatory variables 
#               like MEM eigenfunctions or factors.
#
# Value (elements in the output list) --
#
# License: GPL-2 
# Author:: Pierre Legendre, July 2013
#
{
  # Preliminaries
  X <- as.data.frame(X)
  n <- nrow(X)
  m <- ncol(X)
  out <- matrix(NA,m,4)
  colnames(out) <- c("Rsquare","Adj.Rsquare","F","P.value")
  #
  # Compute a partial RDA for each MEM in presence of all the others (Condition).
  # The P-values are conservative. They correspond to the semipartial R-squares.
  # In Legendre and Legendre (2012, Fig. 14.5), the p-values obtained by forward 
  # selection were used instead of the partial R2 p-values.
  for(j in 1:m) {
    XX <- cbind(X[,j], X[,-j])
    colnames(XX) <- paste("V",1:m, sep="")
    rda.out <- rda(Y ~ V1 + Condition(-V1), data=XX, scale=scale.Y)
    out[j,1] <- RsquareAdj(rda.out)$r.squared
    out[j,2] <- RsquareAdj(rda.out)$adj.r.squared
    temp <- anova(rda.out, step=(nperm+1), perm.max=(nperm+1))
    out[j,3] <- temp[1,3]   # F-statistic
    out[j,4] <- temp[1,4]   # P-value
  }
  
  if(is.null(colnames(X))) {
    rownames(out) <- paste("Var",1:m,sep=".")
  } else {
    rownames(out) <- colnames(X)
  }
  if(!random.var) out <- out[,-2]
  out
}

R2.by.variable_PERSONAL <- function(Y, X, scale.Y=FALSE, nperm=99, random.var=FALSE)
  # Partial R-squares for scalogram -- 
  # Compute R-square and test of significance for each explanatory variable or MEM 
  # eigenfunction in the presence of all others.
  #
  # Arguments --
  # Y : Matrix or data frame of response data.
  # X : Matrix or data frame of explanatory variables (random variables or MEM).
  # scale.Y : =TRUE : standardize the Y variables by columns;
  #           =FALSE : do not standardize the Y variables.
  # nperm : Number of permutations for the tests of the semipartial R2.
  # random.var : =TRUE when X contains random explanatory variables;
#              =FALSE when X contains fixed explanatory variables 
#               like MEM eigenfunctions or factors.
#
# Value (elements in the output list) --
#
# License: GPL-2 
# Author:: Pierre Legendre, July 2013
#
{
  # Preliminaries
  X <- as.data.frame(X)
  n <- nrow(X)
  m <- ncol(X)
  out <- matrix(NA,m,4)
  colnames(out) <- c("Rsquare","Adj.Rsquare","F","P.value")
  #
  # Compute a partial RDA for each MEM in presence of all the others (Condition).
  # The P-values are conservative. They correspond to the semipartial R-squares.
  # In Legendre and Legendre (2012, Fig. 14.5), the p-values obtained by forward 
  # selection were used instead of the partial R2 p-values.
  for(j in 1:m) {
    XX <- cbind(X[,j], X[,-j])
    colnames(XX) <- paste("V",1:m, sep="")
    #rda.out <- rda(Y ~ V1 + Condition(-V1), data=XX, scale=scale.Y)
    rda.out <- dbrda(Y ~ V1 + Condition(-V1), data=XX, scale=scale.Y,dist="bray")
    out[j,1] <- RsquareAdj(rda.out)$r.squared
    out[j,2] <- RsquareAdj(rda.out)$adj.r.squared
    temp <- anova(rda.out, step=(nperm+1), perm.max=(nperm+1))
    out[j,3] <- temp[1,3]   # F-statistic
    out[j,4] <- temp[1,4]   # P-value
  }
  
  if(is.null(colnames(X))) {
    rownames(out) <- paste("Var",1:m,sep=".")
  } else {
    rownames(out) <- colnames(X)
  }
  if(!random.var) out <- out[,-2]
  out
}

# ========

temperature.trend <- function(waterqual, site.ID.tags, sampling) 
{
  n = length(site.ID.tags)
  res = matrix(NA,n,4)
  rownames(res) = paste("Site",site.ID.tags, sep=".")
  colnames(res)= c("slope","Rsquare","adjRsquare","p-value")
  for(i in 1:n) {
    curr.station <- site.ID.tags[i]
    SI <- sampling[sampling$STATION == curr.station, ]
    WQ <- waterqual[sampling$STATION == curr.station, ]
    re <- lm(WQ$WTEMP ~ SI$SAMPLE_DATE)
    res[i,1] = summary(re)$coefficients[2,1]
    res[i,2] = summary(re)$r.squared
    res[i,3] = summary(re)$adj.r.squared
    f.vector = summary(re)$fstatistic
    res[i,4] = pf(f.vector[1],f.vector[2],f.vector[3], lower.tail=FALSE)
  }
  res
}

# ========
