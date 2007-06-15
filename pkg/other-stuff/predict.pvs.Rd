\name{predict.pvs}
\alias{predict.vs}
\title{Classification by pairwise variable selection}
\description{Performs classification of new objects by a classification rule given}
\usage{
predict.pvs(obj, x, quick = FALSE, detail = FALSE, ...)
}
\arguments{
  \item{obj}{An object of class {\em pvs} used for classification of the new data.}
  \item{x}{New data frame or matrix to be classified.}
  \item{quick}{Specifies whether a quick but less accurate estimation should be calculated, 
  if the classification is more important than an estimator of posterior probabilities. Default = FALSE.}
  \item{detail}{}
}
\details{To calculate the posterior probabilities the Pairwise Coupling algorithm (Hastie and Tibshirani, 1998) is used. Here, finding the posterior 
probabilities from pair wise probabilities is formulated as an optimization problem of Kullback-Leibler distance. A simple non-iterative
estimate of the posterior probabilities (Friedman, 1996) can be obtained by setting quick=TRUE, By averaging (and normalizing) all pair wise estimates 
over the classes. This formally leads to the same classification rule but the posterior probability estimates tend to be too close to an equi-probability 
of all classes.  
}
\value{
  \item{classes}{The class levels.}
  \item{prior}{The prior probabilities used.}
  \item{details}{List containing the matrix pair wise posterior probabilities for every object.}
}
\references{Szepannek, G. and Weihs, C., \emph{Variable selection for discrimination of more than two classes where the data is sparse}. In: Proceedings of GfKl 2005 Annual Conference.

Hastie, T. and Tibshirani, R. (1998): Classification by Pairwise Coupling. {\emph Annals of Statistics, 26(1), 451--471}

Friedmann, J. and Tibshirani, R. (1998): Another approach to polychotomous classification. {\emph Technical report, Stanford University}}


\author{Gero Szepannek}
\seealso{\code{\link{pvs}}}

\examples{
library(mlbench)
trainset <- mlbench.waveform(300) # See Breiman et al. for details.
pvsmodel <- pvs(trainset$x, trainset$classes, niveau=0.05) # default: using method="lda"
 
testset <-  mlbench.waveform(500) 
prediction <- predict(pvsmodel, testset$x)

# test error rate
1-sum(testset$classes==prediction$classes)/length(testset$classes)
# Bayes error is 0.149

# comparison to just using lda
ldamodel <- lda(trainset$x, trainset$classes)
LDAprediction <- predict(ldamodel, testset$x)
}
\keyword{classif}
\keyword{multivariate}
