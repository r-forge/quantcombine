\name{collateGSMs}
\alias{collateGSMs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Collate expression value tables for a set of GSMs.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Expression value tables are often available in NCBI's GEO
  database. This function takes a vector of GSM names, downloads the
  expression table for each, and collates them into a single data frame.
}
\usage{
collateGSMs(gsms, sample.names = gsms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsms}{vector of GSM names.
    %%     ~~Describe \code{gsms} here~~
  }
  \item{sample.names}{(Optional) vector of sample names corresponding to
    the GSMs in the gsms argument.
    %%     ~~Describe \code{sample.names} here~~
  }
}
\details{By default, the column names in the output data frame are the
  same as the GSMs. If, for example, group labels are preferred, the
  sample.names argument can be specified.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
  Data frame containing 1 row per probe/probeset and 1 column for each
  GSM/sample.
}
\references{paper?
  %% ~put references to the literature/web site here ~
}
\author{Sadik A. Khuder, Peter S. Bazeley
  %%  ~~who you are~~
}
\note{additional notes?
  %%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##START

# this example obtains data found in data(spermatogenesis)
data(spermatogenesis)

#this table lists the samples/groups; the samples (FileName column)
# will be used
gse6968.gsms # 10 samples using Illumina Sentrix HumanRef8 Beadchip

#get data - this will take several minutes, since each GSM is
# downloaded separately

gse6968 = collateGSMs(gse6968.gsms$FileName)

#show first 10 transcripts
gse6968[1:10,]

#alternatively, to use the group labels for the sample names
gse6968.g = collateGSMs(gse6968.gsms$FileName,
                        sample.names=gse6968.gsms$Group)

gse6968.g[1:10,]

#Note that this second example could be reproduced from the first
# without re-downloading all the samples
gse6968.g = gse6968
names(gse6968.g) = gse6968.gsms$Group

gse6968.g[1:10,]

##END
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
