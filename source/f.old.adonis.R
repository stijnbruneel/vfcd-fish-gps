`adonis2` <-
  function(formula, data, permutations = 999, method = "bray",
           sqrt.dist = FALSE, add = FALSE, by = "terms",
           parallel = getOption("mc.cores"), na.action = na.fail,
           strata = NULL, ...)
  {
    ## handle missing data
    if (missing(data))
      data <- model.frame(delete.response(terms(formula)),
                          na.action = na.action)
    ## we accept only by = "terms", "margin" or NULL
    if (!is.null(by))
      by <- match.arg(by, c("terms", "margin", "onedf"))
    ## evaluate lhs
    YVAR <- formula[[2]]
    lhs <- eval(YVAR, environment(formula), globalenv())
    environment(formula) <- environment()
    ## Take care that input lhs are dissimilarities
    if ((is.matrix(lhs) || is.data.frame(lhs)) &&
        isSymmetric(unname(as.matrix(lhs))))
      lhs <- as.dist(lhs)
    if (!inherits(lhs, "dist"))
      lhs <- vegdist(as.matrix(lhs), method=method, ...)
    ## adjust distances if requested
    if (sqrt.dist)
      lhs <- sqrt(lhs)
    if (is.logical(add) && isTRUE(add))
      add <- "lingoes"
    if (is.character(add)) {
      add <- match.arg(add, c("lingoes", "cailliez"))
      if (add == "lingoes") {
        ac <- addLingoes(as.matrix(lhs))
        lhs <- sqrt(lhs^2 + 2 * ac)
      }
      else if (add == "cailliez") {
        ac <- addCailliez(as.matrix(lhs))
        lhs <- lhs + ac
      }
    }
    ## adonis0 & anova.cca should see only dissimilarities (lhs)
    if (!missing(data)) # expand and check terms
      formula <- terms(formula, data=data)
    if (is.null(attr(data, "terms"))) # not yet a model.frame?
      data <- model.frame(delete.response(terms(formula)), data,
                          na.action = na.action)
    formula <- update(formula, lhs ~ .)
    sol <- adonis0(formula, data = data, method = method)
    ## handle permutations
    perm <- getPermuteMatrix(permutations, NROW(data), strata = strata)
    out <- anova(sol, permutations = perm, by = by,
                 parallel = parallel)
    ## attributes will be lost when adding a new column
    att <- attributes(out)
    ## add traditional adonis output on R2
    out <- rbind(out, "Total" = c(nobs(sol)-1, sol$tot.chi, NA, NA))
    out <- cbind(out[,1:2], "R2" = out[,2]/sol$tot.chi, out[,3:4])
    ## Fix output header to show the adonis2() call instead of adonis0()
    att$heading[2] <- deparse(match.call(), width.cutoff = 500L)
    att$names <- names(out)
    att$row.names <- rownames(out)
    attributes(out) <- att
    out
  }

`adonis0` <-
  function(formula, data=NULL, method="bray")
  {
    ## First we collect info for the uppermost level of the analysed
    ## object
    Trms <- terms(data)
    sol <- list(call = match.call(),
                method = "adonis",
                terms = Trms,
                terminfo = list(terms = Trms))
    sol$call$formula <- formula(Trms)
    TOL <- 1e-7
    lhs <- formula[[2]]
    lhs <- eval(lhs, environment(formula)) # to force evaluation
    formula[[2]] <- NULL                # to remove the lhs
    rhs <- model.matrix(formula, data) # and finally the model.matrix
    assign <- attr(rhs, "assign") ## assign attribute
    sol$terminfo$assign <- assign[assign > 0]
    rhs <- rhs[,-1, drop=FALSE] # remove the (Intercept) to get rank right
    rhs <- scale(rhs, scale = FALSE, center = TRUE) # center
    qrhs <- qr(rhs)
    ## input lhs should always be dissimilarities
    if (!inherits(lhs, "dist"))
      stop("internal error: contact developers")
    if (any(lhs < -TOL))
      stop("dissimilarities must be non-negative")
    ## if there was an na.action for rhs, we must remove the same rows
    ## and columns from the lhs (initDBRDA later will work similarly
    ## for distances and matrices of distances).
    if (!is.null(nas <- na.action(data))) {
      lhs <- as.matrix(lhs)[-nas,-nas, drop=FALSE]
      n <- nrow(lhs)
    } else
      n <- attr(lhs, "Size")
    ## G is -dmat/2 centred
    G <- initDBRDA(lhs)
    ## preliminaries are over: start working
    Gfit <- qr.fitted(qrhs, G)
    Gres <- qr.resid(qrhs, G)
    ## collect data for the fit
    if(!is.null(qrhs$rank) && qrhs$rank > 0)
      CCA <- list(rank = qrhs$rank,
                  qrank = qrhs$rank,
                  tot.chi = sum(diag(Gfit)),
                  QR = qrhs)
    else
      CCA <- NULL # empty model
    ## collect data for the residuals
    CA <- list(rank = n - max(qrhs$rank, 0) - 1,
               u = matrix(0, nrow=n),
               tot.chi = sum(diag(Gres)))
    ## all together
    sol$tot.chi <- sum(diag(G))
    sol$adjust <- 1
    sol$Ybar <- G
    sol$CCA <- CCA
    sol$CA <- CA
    class(sol) <- c("adonis2", "dbrda", "rda", "cca")
    sol
  }

`initDBRDA` <-
  function(Y)
  {
    ## check
    Y <- as.matrix(Y)
    dims <- dim(Y)
    if (dims[1] != dims[2] || !isSymmetric(unname(Y)))
      stop("input Y must be distances or a symmetric square matrix")
    ## transform
    Y <- -0.5 * GowerDblcen(Y^2)
    attr(Y, "METHOD") <- "DISTBASED"
    Y
  }

GowerDblcen <- function(x, na.rm = TRUE)
{
  cnt <- colMeans(x, na.rm = na.rm)
  x <- sweep(x, 2L, cnt, check.margin = FALSE)
  cnt <- rowMeans(x, na.rm = na.rm)
  sweep(x, 1L, cnt, check.margin = FALSE)
}

### Internal functions to find additive constants to non-diagonal
### dissimilarities so that there are no negative eigenvalues. The
### Cailliez constant is added to dissimilarities and the Lingoes
### constant is added to squared dissimilarities. Legendre & Anderson
### (Ecol Monogr 69, 1-24; 1999) recommend Lingoes, but
### stats::cmdscale() only provides Cailliez. Input parameters: d are
### a matrix of dissimilarities.

addCailliez <- function(d)
{
  n <- nrow(d)
  q1 <- seq_len(n)
  q2 <- n + q1
  ## Cailliez makes a 2x2 block matrix with blocks of n x n elements.
  ## Blocks anti-clockwise, upper left [0]
  z <- matrix(0, 2*n, 2*n)
  diag(z[q2,q1]) <- -1
  z[q1,q2] <- -GowerDblcen(d^2)
  z[q2,q2] <- GowerDblcen(2 * d)
  ## Largest real eigenvalue
  e <- eigen(z, symmetric = FALSE, only.values = TRUE)$values
  out <- max(Re(e))
  max(out, 0)
}

addLingoes <- function(d)
{
  ## smallest negative eigenvalue (or zero)
  d <- -GowerDblcen(d^2)/2
  e <- eigen(d, symmetric = TRUE, only.values = TRUE)$values
  out <- min(e)
  max(-out, 0)
}

`getPermuteMatrix` <-
  function(perm, N,  strata = NULL)
  {
    ## 'perm' is either a single number, a how() structure or a
    ## permutation matrix
    if (length(perm) == 1) {
      perm <- how(nperm = perm)
    }
    ## apply 'strata', but only if possible: ignore silently other cases
    if (!missing(strata) && !is.null(strata)) {
      if (inherits(perm, "how") && is.null(getBlocks(perm)))
        setBlocks(perm) <- strata
    }
    ## now 'perm' is either a how() or a matrix
    if (inherits(perm, "how"))
      perm <- shuffleSet(N, control = perm)
    else { # matrix: check that it *strictly* integer
      if(!is.integer(perm) && !all(perm == round(perm)))
        stop("permutation matrix must be strictly integers: use round()")
    }
    ## now 'perm' is a matrix (or always was). If it is a plain
    ## matrix, set minimal attributes for printing. This is a dirty
    ## kluge: should be handled more cleanly.
    if (is.null(attr(perm, "control")))
      attr(perm, "control") <-
        structure(list(within=list(type="supplied matrix"),
                       nperm = nrow(perm)), class = "how")
    perm
  }