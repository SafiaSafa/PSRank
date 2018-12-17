\name{PSRank}
\alias{PSRank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Evaluate the predictivity GWAS of according to the ranks of SNPs by different methods and a "polygenic risk score".
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The population is separate into two sets of data:
    A learning set where will be estimated associations through linear univariate model. This set of SNPs is going to be sorted. The packages allows to sorte the SNPs according differentes scores. By P-values, β or I-score.
    And a test set where will be calculated a cumulative score, called "polygenic risk score"  and on which we will evaluate the predictivity
}
\usage{
PSRank(plink,
       frac,
       ranking,
       method,
       sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plink}{
  Filename of the PLINK binary file (BED/BIM/FAM)
}
  \item{frac}{
  Fraction of the dataset that's used for training
}
  \item{ranking}{
  What method should be used to rank the SNPs: "Z" for ranking by Z-score, "beta" for ranking by Beta, "iscore" for ranking by I-score and "combined" for ranking by Z-score+I-score
}
  \item{method}{
  Cross-validation method (subsampling of kfold)
}
  \item{sample}{
  Number of replication
}
}

\value{

This function return a matrix of AUC values for each SNPs. The number of rows in the matrix depends on the number of replication.
}

\author{
NOIREL Josselin, SAFA-TAHAR-HENNI Safia
}
