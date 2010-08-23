\name{getQuantileScores}
\alias{getQuantileScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Assigns quantile score to each sample.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  For each gene, generate quantile scores for each sample based on the
  distribution of expression values.
}
\usage{
getQuantileScores(exprs, grp1, grp2, n.quantiles = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprs}{data frame or matrix of expression values, with
    probes/probesets as rows and samples as columns.
%%     ~~Describe \code{exprs} here~~
}
  \item{grp1}{vector of indices or sample names for study group 1.
%%     ~~Describe \code{grp1} here~~
}
  \item{grp2}{vector of indices or sample names for study group 2.
%%     ~~Describe \code{grp2} here~~
}
  \item{n.quantiles}{number of quantiles used to distribute sample
    expression values for each gene.
%%     ~~Describe \code{n.quantiles} here~~
}
}
\details{
  n.quantiles quantiles are generated using all sample expression values
  for a given gene. Using these quantiles as thresholds, each sample is then
  assigned one of n.quantiles+1 quantiles scores, depending on which
  quantile thresholds the sample's expression value is >= to and <=
  to. Quantile scores are centered on 0.

  The quantile algorithm used is the same as the default used by the
  quantile() function from the stats package (type 7). The number of
  probabilities is equal to n.quantiles. Each probability is equal
  to 1/(n.quantiles+1) x i, where i is the index position in the
  probabilities vector. For example, if n.quantiles=4, the
  probabilities are 0.2, 0.4, 0.6, and 0.8.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
  Data frame with probes/probesets as rows and samples as columns, where
  each expression value has been replaced with a quantile score.
}
\references{which paper?
%% ~put references to the literature/web site here ~
}
\author{Sadik A. Khuder, Peter S. Bazeley
%%  ~~who you are~~
}
\note{Additional notes?
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{do we want anything here?
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##START
library(QuantCombine)

#Example 1

exprs = data.frame(matrix(rnorm(50),ncol=10))
exprs

#arbitrary groups in this case
grp1 = 1:5
grp2 = 6:10

sc = getQuantileScores(exprs,grp1=grp1,grp2=grp2)
sc

#Example 2 - use 6 quantiles instead of the default 4

sc = getQuantileScores(exprs,grp1,grp2,n.quantiles=6)
sc


#Example 3 - use data from NCBI GEO

library(GEOquery)

#GSE6872 represents 21 samples hybridized to Affymetrix HGU133 Plus 2.0
# arrays. 
gse = getGEO("GSE6872",GSEMatrix=T)

#to view the organization of the gse object:
#str(gse)

#the sample groups
grps = gse[[1]]@phenoData@data$title
grps

grp1 = grep("Normo",grps)
grp1
grp2 = grep("Terato",grps)
grp2

#the expression data
exprs = gse[[1]]@assayData$exprs
exprs[1:5,]

sc = getQuantileScores(exprs,grp1,grp2)
sc[1:5,]

##END
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
