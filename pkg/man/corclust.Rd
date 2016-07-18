\name{corclust}
\alias{corclust}
\title{Function to identify groups of highly correlated variables 
for removing correlated features from the data for further analysis.}
\description{
A hierarchical clustering of variables using \code{hclust} is performed using 
1 - the absolute correlation as a distance measure between tow variables.  

}
\usage{
corclust(x, cl = NULL, method = "complete")
}
\arguments{
  \item{x}{Either a data frame or a matrix consisting of numerical attributes.}
  \item{cl}{Optional vector of ty factor indicating class levels, if class specific correlations should to be considered.}
  \item{method}{Linkage to be used for clustering. Default is \code{complete} linkage.}
}
\details{
Each cluster consists of a set of correlated variables according to the chosen clustering criterion. 
The default criterion is \sQuote{\code{complete}}. This choice is meaningful as it represents the 
\emph{minimum absolute correlation} between all variables of a cluster.\\
The data set is split into numerics and factors two separate clustering models are built, depending on the variable type. 
For factors distances are computed based on 1-Cramer's V statistic using \code{\link{chisq.test}}. 
For a large number of factor variables this might take some time. 
The resulting trees can be plotted using \code{\link{plot.corclust}}.\\
Further proceeding would consist in chosing one variable of each cluster to obtain a 
subset of rather uncorrelated variables for further analysis. An automatic variable selection can be done using \code{\link{cvtree}} 
and \code{\link{xtractvars}}.\\
If an additional class vector \code{cl} is given to the function for any two variables their minimum 
correlation over all classes is used.  
}
\value{Object of class \code{corclust}.
  \item{cor}{Correlation matrix of numeric variables.}
  \item{crv}{Matrix of Cramer's V for factor variables.}
  \item{cluster.numerics}{Resulting hierarchical \code{hclust} model for numeric variables.}
  \item{cluster.factors}{Resulting hierarchical \code{hclust} model for factor variables.}
  \item{id.numerics}{Variable IDs of numeric variables in \code{x}.}
  \item{id.factors}{Variable IDs of factor variables \code{x}.}

}

\references{Roever, C.and Szepannek, G. (2005): Application of a genetic algorithm
to variable selection in fuzzy clustering. In C. Weihs and W. Gaul (eds), Classification 
- The Ubiquitous Challenge, 674--681, Springer.}


\author{Gero Szepannek}

\seealso{See also \code{\link{plot.corclust}}, \code{\link{cvtree}}, and \code{\link{xtractvars}} as well as \code{\link{hclust}}, for details on the clustering algorithm.}
\examples{
    data(iris)
    classes <- iris$Species
    variables <- iris[,1:4]
    ccres <- corclust(variables, classes)
    plot(ccres, mincor = 0.6)
    
}

\keyword{cluster}
\keyword{multivariate}
\keyword{classif}
\keyword{manip}
\keyword{attribute}
