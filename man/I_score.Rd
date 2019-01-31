\name{I_score}
\alias{I_score}
\title{Function to calculate the vector of I scores from a genotype matrix and a phenotype vector}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  This function using vectorised computations to speed up the calculation of the
  vector of I scores from a genotype matrix and a phenotype vector.
}
\usage{
  I_score(genotypes, phenotypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotypes}{a $N \times k$ matrix containing genotypes encoded as 0, 1, 2 or NA
  }
  \item{phenotypes}{a vector of $k$ values for the phenotype (quantitative or case/control)}
}
\value{
  Return a vector containing the I-score values associated with each SNP
}
\examples{
case_probs = c(.274, .5, .226)
ctrl_probs = rev(case_probs)

n = 500   # Number of cases and controls; cohort is twice that
p = 10000 # Number of SNPs; replicates

genotypes =
  matrix(c(sample(0:2, size=n * p, prob=ctrl_probs, replace=TRUE),
           sample(0:2, size=n * p, prob=case_probs, replace=TRUE)),
         byrow=TRUE, ncol=p)

phenotypes = rep(0:1, each=n)

I_samples = I_score(genotypes, phenotypes)
}
\references{
  Lo et al. Why significant variables aren't automatically good predictors, PNAS 2015 112 (45) 13892-13897
}
\author{
  Josselin Noirel, Safia Safa-tahar-henni
}
