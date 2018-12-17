\name{I_score}
\alias{I_score}

\title{
  function to calculate an I-score for each SNP. The I-score :an alternative statistic to the X2 designed to recognize strongly predictive variables.
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
    should be a $S \times N$ matrix containing 0, 1, 2 or NA;
  }
  \item{phenotypes}{
    should be a vector of $S$ values for the phenotype; use numeric values for a quantitative trait and booleans for a case/control association study.
  }
}

\value{
  Return a matrix containing the I-score values associated with each SNP
}
\references{
  "Why significant variables aren{\textquoteright}t automatically good predictors"
  Adeline Lo, Herman Chernoff, Tian Zheng, Shaw-Hwa Lo
  Proceedings of the National Academy of Sciences Nov 2015, 112 (45) 13892-13897; DOI: 10.1073/pnas.1518285112

}
\author{
  NOIREL Josselin, SAFA-TAHAR-HENNI Safia
}
