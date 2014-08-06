cond.index <- function(formula, data, ...)
{
reg<-lm(formula, data, x=TRUE, ...)
x.scale <- reg$x%*%diag(1/sqrt(apply(reg$x^2,2,sum)))
x.scale.svd<- svd(t(x.scale)%*%x.scale)
ci <- sqrt(max(x.scale.svd$d)/x.scale.svd$d)
return(ci)
}