\name{pvs}
\alias{pvs}
\alias{pvs.default}
\alias{pvs.formula}
\alias{print.pvs}
\title{Classification by pairwise variable selection}
\description{Computes a classification model for any classification method by 
    performing class pair wise variable selection and thus a maximal 
    dimension reduction to avoid sparsity of the data.}
\usage{
pvs(x, ...)
\method{pvs}{default}(x, grouping, prior=NULL, method="lda", 
        vs.method=c("ks-test","stepclass","greedy.wilks"), 
        niveau=0.06, fold=10, impr=0.1, direct="backward", out=FALSE, ...)
\method{pvs}{formula}(formula, data=NULL, ...)
}
\arguments{
  \item{x}{A matrix or data frame containing the explanatory variables.}
  \item{grouping}{A factor specifying the class for each observation.}
  \item{prior}{The prior probabilities of class membership. If unspecified, the class proportions for the training set are used. If present, the probabilities should be specified in the order of the factor levels.}
  \item{formula}{A formula of the form groups ~ x1 + x2 + ... That is, the response is the grouping factor and the right hand side specifies the (non-factor) discriminators.}
  \item{data}{Data frame from which variables specified in formula are preferentially to be taken.}
  
  \item{method}{The desired classification method to be used.}
  \item{vs.method}{The desired variable selection procedure that should be implemented, either "ks-test", "stepclass" or "greedy.wilks".}
  \item{niveau}{Specifying p value for using ks-test or F test (using Wilk's lambda) as criterion for variable selection.}
  \item{fold}{Number of data partitions for cross-validation when using stepclass variable selection.}
  \item{impr}{Specifying the necessary percentage of improvement of an additional variable to be selected into the model.}
  \item{direct}{Specifying the direction of variable selection. One of "forward", "backward" or "both".}
  \item{output}{Specifying if an output (if vs.method=="stepclass" is chosen) should be printed while computing.}
  \item{\dots}{Further arguments passed to the used classification method.}
}
\details{
}
\value{
  Returns an object of class \emph{pvs}. 
  \item{classes}{The class levels.}
  \item{prior}{The prior probabilities used.}
  \item{method}{The implemented classification method.}
  \item{vs.method}{The implemented variable selection method.}
  \item{submodels}{A list of $classpair $subspace $model containing detailed information about the class pair wise variable subset and their associated classification models.}
}
\references{Szepannek, G. and Weihs, C., \emph{Variable selection for discrimination of more than two classes where the data is sparse}. In: Proceedings of GfKl 2005 Annual Conference.}
To calculate the posterior probabilities the Pairwise Coupling algorithm (Hastie and Tibshirani, 1998) is used.

Hastie, T. and Tibshirani, R. (1998): Classification by Pairwise Coupling. {\emph Annals of Statistics, 26(1), 451--471}

Breiman, L. Friedman, J., Olshen, R., and Stone, C. (1984):
Classification and regression trees. {\emph Chapman \& Hall, NY}
}

\author{Gero Szepannek}
\seealso{\code{\link{predict.pvs}}, \code{\link{stepclass}}}

\examples{
library(mlbench)
library(MASS)
trainset <- mlbench.waveform(300) # See Breiman et al. for details.
pvsmodel <- pvs(trainset$x, trainset$classes, niveau=0.05) # default: using method="lda"
 
trainset2 <- data.frame(trainset$x, classes=trainset$classes)
pvsmodel2 <- pvs(classes ~ ., data = trainset2, niveau=0.05) # default: using method="lda"
 
testset <-  mlbench.waveform(500) 
prediction <- predict(pvsmodel, testset$x)

testset2 <- data.frame(testset$x, classes=testset$classes)
prediction2 <- predict(pvsmodel2, newdata = testset2)
 
# test error rate
1-sum(testset$classes==prediction$classes)/length(testset$classes)
# Bayes error is 0.149

# comparison to just using lda
ldamodel <- lda(trainset$x, trainset$classes)
LDAprediction <- predict(ldamodel, testset$x)
}
\keyword{classif}
\keyword{multivariate}
