pvs <- function(x, ...) 
{
  UseMethod("pvs")
}

pvs.default <- function(x, grouping, prior=NULL, method="lda", vs.method=c("ks-test","stepclass","greedy.wilks"), 
        niveau=0.06, fold=10, impr=0.1, direct="backward", out=FALSE, ...){ 
    cl <- match.call()
    cl[[1]] <- as.name("pvs")
    vs.method <- match.arg(vs.method)
    classes <- levels(grouping)
    if(length(classes) < 2) stop("at least 2 classes required")
    classcombins <- t(combinations(length(classes), 2, classes))
    if(is.null(prior)) prior <- table(grouping)/sum(table(grouping))
    names(prior) <- classes
    
    ## for a given pair of classes (classvec) test whether variables differ,  
    ## estimate distribution of resulting subspace for both classes
    calc.classcompare <- function(classvec, pval=niveau, ...){
        class1 <- which(grouping==classvec[1])
        class2 <- which(grouping==classvec[2])     
        apriori <- as.vector(prior[classvec]/sum(prior[classvec]))   
        
        switch(vs.method,
            "ks-test" = {
                pvalues <- apply(x, 2, function(vec) 
                    ks.test(vec[class1], vec[class2], alternative="two.sided", exact=TRUE)$p.value)
                relevant <- pvalues <= pval
                ## if all variables are discarded, build model of those with minimal p.values...
                if(!sum(relevant)) relevant <- pvalues <= min(pvalues)  
            
                subgrouping <- factor(grouping[c(class1, class2)], labels=classvec)  
                model <- try(do.call(method, list(as.matrix(x[c(class1, class2), relevant]), 
                    subgrouping, prior=apriori, ...)))
                if(inherits(model, "try-error"))
                    stop("method ", sQuote(method), " resulted in the error reported above")
                return(list(classpair=classvec, subspace=list(variables=which(relevant), p.values=pvalues), 
                    model=model))
                },
            "stepclass" = {
                which.relevant <- stepclass(x=x[c(class1, class2), ], 
                    grouping=grouping[c(class1, class2)], method=method, prior=apriori, fold=fold, improve=impr, 
                    direction=direct, output=out, ...)   
                details <- which.relevant$process
                which.relevant <- which.relevant$model$nr
                relevant <- !logical(ncol(x)) # in clude all variables (ovoid no variable will be selected)
                if(length(which.relevant)) relevant[-which.relevant] <- FALSE # if any variable selected
            
                subgrouping <- factor(grouping[c(class1, class2)],  labels=classvec)  
                model <- try(do.call(method, list(as.matrix(x[c(class1, class2), relevant]), 
                    subgrouping, prior=apriori, ...)))
                if(inherits(model, "try-error"))
                    stop("method ", sQuote(method), " resulted in the error reported above")
                return(list(classpair=classvec, subspace=list(variables=which(relevant), details=details),  
                    model=model))
                },
            "greedy.wilks" = {
                which.relevant <- greedy.wilks(X=x[c(class1, class2), ], grouping=grouping[c(class1, class2)], 
                    method=method, prior=apriori, niveau=niveau)
                results <- which.relevant$results
                which.relevant <- as.numeric(which.relevant[[1]]$vars) # indexes of chosen variables
                relevant <- !logical(ncol(x)) # in clude all variables (ovoid no variable will be selected)
                if(length(which.relevant)) relevant[-which.relevant] <- FALSE # if any variable selected
                model <- try(do.call(method, list(as.matrix(x[c(class1, class2), relevant]), 
                    subgrouping, prior=apriori, ...)))
                if(inherits(model, "try-error"))
                    stop("method ", sQuote(method), " resulted in the error reported above")
                subgrouping <- factor(grouping[c(class1, class2)],  labels=classvec)  
                return(list(classpair=classvec,  subspace=list(variables=which(relevant), results=results), 
                    model=model))     
                } 
        )
    }
        
    models <- apply(classcombins, 2, calc.classcompare, ...)
    result <- list(classes=classes, prior=prior, method=method, vs.method=vs.method, submodels=models, call=cl)    
    class(result) <- "pvs"
    return(result)
}


pvs.formula <- function(formula, data = NULL, ...)
{
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
        x <- x[, -xint, drop = FALSE]
    res <- pvs(x, grouping, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("pvs")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}



predict.pvs <- function(object, newdata, quick = FALSE, detail = FALSE, ...){
  # calculates for 1 n
    pred.obj <- function(newob, ob=object){
        comp.result <- diag(0, length(ob$classes))
        colnames(comp.result) <- rownames(comp.result) <- ob$classes 
        for(i in ob$submodels){  
            pairclas <- try(predict(i$model, newob[i$subspace$variables]), silent=TRUE)
            if(is.list(pairclas)) pairclas<-pairclas$posterior                        # weighted by prior for classpair
            for(j in 1:2) comp.result[colnames(pairclas)[j], colnames(pairclas)[3-j]] <- pairclas[j]     
        } # compares every pair of classes (of obb$modell) and replaces 0 by one following: rows indicate 'winner',  coloumns loser          
        return(comp.matrix=comp.result)
    }  
    
      
    classifier <- function(comp.mat, quick = FALSE){
        comp.mat<-matrix(comp.mat, ncol=length(object$classes))       
        colnames(comp.mat) <- rownames(comp.mat) <- object$classes     
        ## quick-and-dirty-estimate if chosen in function call quick==TRUE
        result <- if(quick) rowMeans(comp.mat) else pfromp(comp.mat)
        return(result)
    }         

    if (!inherits(object, "pvs")) 
        stop("object not of class", " 'pvs'")
    if (!is.null(Terms <- object$terms)) {
        if (missing(newdata)) 
            newdata <- model.frame(object)
        else {
            newdata <- model.frame(as.formula(delete.response(Terms)), 
                newdata, na.action = function(x) x, xlev = object$xlevels)
        }
        x <- model.matrix(delete.response(Terms), newdata, contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(x), nomatch = 0)
        if (xint > 0) 
            x <- x[, -xint, drop = FALSE]
    }
    else {
        if (missing(newdata)) {
            if (!is.null(sub <- object$call$subset)) 
                newdataa <- eval.parent(parse(text = paste(deparse(object$call$x, 
                  backtick = TRUE), "[", deparse(sub, backtick = TRUE), 
                  ",]")))
            else newdata <- eval.parent(object$call$x)
            if (!is.null(nas <- object$call$na.action)) 
                newdata <- eval(call(nas, newdata))
        }
        if (is.null(dim(newdata))) 
            dim(newdata) <- c(1, length(newdata))
        x <- as.matrix(newdata)
    }


    pairwise.comparisons <- apply(x, 1, pred.obj, ob=object)

    posterior <- apply(pairwise.comparisons, 2, classifier, quick = quick)
    rownames(posterior) <- object$classes
    classes <- apply(posterior, 2, function(x) return(names(x)[which.max(x)]))
    details <- vector(mode="list", length=ncol(pairwise.comparisons))
    for(i in 1:ncol(pairwise.comparisons)) {
            details[[i]] <- matrix(pairwise.comparisons[, i], ncol=length(object$classes), nrow=length(object$classes))
            colnames(details[[i]])  <- rownames(details[[i]]) <- as.factor(object$classes)
    }        
    return(c(list(class=factor(classes), posterior=t(posterior)), if(detail) list(details=details)))
}




#Function for estimating class probabilities from probabilities of pairwise KLOGREG
#see: Hastie,  Tibshirani: Classification by pairwise coupling. In Jordan,  Kearns,  Solla,  editors, 
#Advances in Neural Information Processing Systems,  volume 10. The MIT Press,  1998.
#Author: Marcos Marin-Galiano,  Dept. Of Statistics,  University Of Dortmund,  Germany
#        modified by Gero Szepannek
#email: marcos.marin-galiano@uni-dortmund.de.
#input: matr is the matrix of r_ij = p_i / (p_i + p_j). matr must have the same number of
#rows and columns,  matnum is the matrix of the n_ij
pfromp <- function(matr, tolerance=1.E-4, matnum=matrix(rep(1, dim(matr)[1]^2), nrow=dim(matr)[1]))
{
    classnum<-dim(matr)[1]
    #initial set for the matrix of mu´s.
    matmue<-matrix(0, nrow=classnum, ncol=classnum)
    #initial estimate for class probabilities = laplace probabilities
    pestnew<-rep(1/classnum, classnum)
    wmatr<-diag(matnum%*%t(matr))


    #algorithm loop
    repeat
        {pest<-pestnew

        #computation of new matrix mu and the new estimate for p
        for (i in 2:classnum){
                temp <- pest[i]/(pest[i]+pest[1:(i-1)])
                matmue[i, 1:(i-1)] <- temp
                matmue[1:(i-1), i] <- 1 - temp
        }
        ##normalization of p
        divisor <- diag(matnum%*%t(matmue))
        pestnew<-pest*wmatr/divisor
        if(any(divisor==0)) pestnew[divisor==0] <- 0
        pestnew <- pestnew/sum(pestnew)
        #breaking rule
        if (sum(abs(pestnew-pest))<tolerance) break
    }
    return(pestnew)
}


print.pvs <- function(x, ...){
    cat("Used classifier: ", x$method, "\n\n")
    dummy <- x$vs.method
    if(is.null(dummy)) dummy <- "Kolmogorov Smirnov - test"
    cat("Used variable selection: ", dummy, "\n\n")
    cat("Pairwise subspaces: \n")
    lapply(x$submodels, function(x) cat("classes: ", x$classpair, "\t variable subset: ", x$subspace$variables, "\n"))
    invisible(x)
}
