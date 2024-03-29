predict.NaiveBayes <- function (object, newdata, threshold = 0.001, ...)
{
    if (missing(newdata))
        newdata <- object$x
    if (sum(is.element(colnames(newdata), object$varnames)) < length(object$varnames))
        stop("Not all variable names used in object found in newdata") 
    ## (both colnames & varnames are given) & (varnames is a subset of colnames):
    newdata <- data.frame(newdata[, object$varnames])
    nattribs <- ncol(newdata)
    islogical <- sapply(newdata, is.logical)
    isnumeric <- sapply(newdata, is.numeric)
#   L <- sapply(
#      1:nrow(newdata),
#      function(i)
#      {
#         ndata <- as.numeric(newdata[i, ])
#         L <-  sapply(
#               1:nattribs,
#               function(v)
#               {
#                  nd <- ndata[v]
#                  if (is.na(nd))
#                  {
#                     rep(1, length(object$apriori))
#                  } else {
#                     prob <- if (isnumeric[v])
#                     {
#                        msd <- object$tables[[v]]
#                        if (object$usekernel) sapply(
#                           msd,
#                           FUN = function(y)
#                           {
#                              dkernel(x = nd, kernel = y, ...)
#                           })
#                        else dnorm(nd, msd[, 1], msd[, 2])
#                     }
#                     else object$tables[[v]][, nd]
#
#                     prob
#                  }
#               }
#               )
#
#         L <- ifelse(L < threshold, threshold, L)
#         # normalize by p(x) = p(x_1|y) + ... + p(x_p|y)
#         Lnorm <- apply(L, 2, function(x, y) x/sum(x * y), y = as.vector(object$apriori))
#         # get product
#         Lprod <- apply(Lnorm, 1, prod)
#         # normalize by posterior
#         Lpost <- object$apriori * Lprod
#         Lpost <- Lpost/sum(Lpost)
#         Lpost
#      }
#   )
    isfactor <- sapply(object$x, is.factor)
    if(any(isfactor)){
        factors <- object$varnames[isfactor]
        newdata[factors] <- lapply(factors, function(f) {
            undf <- unique(newdata[[f]])
            contained <- undf %in% levels(object$x[[f]])
            if(any(!contained)){
                newlevels <- undf[!contained]
                warning(paste0("NAs produced for new level(s) (", paste(sQuote(newlevels, "'"), collapse=", "), ") of newdata$", f))
            }
            factor(newdata[[f]], levels=levels(object$x[[f]]))
        })
    }

    
    newdata <- data.matrix(newdata)
    Lfoo <- function(i) {
        tempfoo <- function(v) {
            nd <- ndata[v]
            if (is.na(nd))
                return(rep(1, length(object$apriori)))
            prob <-
                if (isnumeric[v]) {
                    msd <- object$tables[[v]]
                    if (object$usekernel)
                        sapply(msd, FUN = function(y)
                            dkernel(x = nd, kernel = y, ...))
                    else dnorm(nd, msd[, 1], msd[, 2])
                } else if (islogical[v]) {
                    object$tables[[v]][, nd + 1]
                } else {
                    object$tables[[v]][, nd]
                }
            prob[prob == 0] <- threshold
            return(prob)
        }

        ndata <- newdata[i, ]
        tempres <- log(sapply(1:nattribs, tempfoo))
        L <- log(object$apriori) + rowSums(tempres)

#        L <- exp(L)
#        L/sum(L)

        if(isTRUE(all.equal(sum(exp(L)), 0)))
            warning("Numerical 0 probability for all classes with observation ", i)
        L
    }
    L <- sapply(1:nrow(newdata), Lfoo)

    classdach <- factor(object$levels[apply(L, 2, which.max)],
                        levels = object$levels)
    posterior <- t(apply(exp(L), 2, function(x) x/sum(x)))

#                  print(str(posterior))
    colnames(posterior) <- object$levels
    rownames(posterior) <- names(classdach) <- rownames(newdata)
    return(list(class = classdach, posterior = posterior))
}
