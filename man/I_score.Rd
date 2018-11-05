\name{I_score}
\alias{I_score}

\title{
%%  ~~function to do ... ~~
The I-score :an alternative statistic to the X2 designed to recognize strongly predictive variables
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function makes it possible to define a score for each SNP according to a genotype matrix and a phenotype vector. This score reflects the predictivity of this SNP for this phenotype.
}
\usage{
I_score(genotypes, phenotypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{
%%     ~~Describe \code{x} here~~
should be a $S \times N$ matrix containing 0, 1, 2 or NA;
}
  \item{phenotypes}{
should be a vector of $S$ values for the phenotype; use numeric values for a quantitative trait and booleans for a case/control association study.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Return a matrix containing the I-score values associated with each SNP
}
\references{
%% ~put references to the literature/web site here ~

"Why significant variables aren{\textquoteright}t automatically good predictors"
Adeline Lo, Herman Chernoff, Tian Zheng, Shaw-Hwa Lo
Proceedings of the National Academy of Sciences Nov 2015, 112 (45) 13892-13897; DOI: 10.1073/pnas.1518285112

}
\author{
NOIREL Josselin, SAFA-TAHAR-HENNI Safia
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function function (genotypes, phenotypes) {
  b  = ! is.na(genotypes)
  Ym = colSums(b * phenotypes)/colSums(b)

  n0 = colSums((M0 = genotypes == 0), na.rm=TRUE)
  n1 = colSums((M1 = genotypes == 1), na.rm=TRUE)
  n2 = colSums((M2 = genotypes == 2), na.rm=TRUE)
  n  = n0 + n1 + n2

  S  = (colSums(phenotypes * M0, na.rm=TRUE) - n0 * Ym)^2 +
    (colSums(phenotypes * M1, na.rm=TRUE) - n1 * Ym)^2 +
    (colSums(phenotypes * M2, na.rm=TRUE) - n2 * Ym)^2

  S/sum((phenotypes - weighted.mean(Ym, n))^2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
