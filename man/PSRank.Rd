\name{PSRank}
\alias{PSRank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
This functon do GWA analysis.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
We split the data into two subsets: training and testing. The associations are calculated on the trainig set. SNPs are ranked according to the different scores. We then calculate the polygenic score. AUC values are then computed. To obtain AUC values that are representative of the observed phenotype, the analyzes are performed 100 times.
}
\usage{
PSRank(plink="filename",
       frac=.5,
       ranking="Z")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frac}{
%%     ~~Describe \code{x} here~~
Fraction of the dataset that's used for training
}
  \item{ranking}{
What method should be used to rank
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
We split the data into two subsets: training and testing.  It is done in such a way as to preserve the number of cases and controls in all simulations.The associations are calculated on the trainig set. SNPs are ranked according to the different scores :$p$-values (or equivalently the $Z$ score), $beta, $I$-score or the combination of two scores: $Z$ score + $I$-score. If there are missing data, we replace the missing values with the expected values (since we assume independence); this could be achieved by centring the genotypes matrix and replacing missing values by zeroes. We then calculate the polygenic score. AUC values are then computed. To obtain AUC values that are representative of the observed phenotype, the analyzes are performed 100 times.
The use of the polygenic score and a method of evaluating predictive models, the AUC. To evaluate the associations obtained according to the different classification methods. This is to determine which methods can improve predictivity.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
This function return a matrix of AUC values for each. It's a matrix of 100 rows for each replication.
}
\references{
%% ~put references to the literature/web site here ~
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

library("pROC")
library("matrixStats")
library("caTools")
library("snpStats")
library("ggplot2")
fam = "output/unif_Laplace_omni.fam"
bim = "output/unif_Laplace_omni.bim"
bed = "output/unif_Laplace_omni.bed"
dat = read.plink(bed, bim, fam)

frac = .5

res = PSRank(dat,
             frac=.5,
             ranking="beta")
# Write CSV in R
write.table(res, file = "res_AUC.csv",row.names=FALSE, col.names=FALSE, sep=",")
# Read CSV into R
MyData <- read.csv(file="res_AUC.csv", header=FALSE, sep=",")



# plot means
summary_res = as.data.frame(cbind(SNPs = 1:dim(res)[2],Means = colMeans(res),sd = colSds(res)))

ggplot(summary_res, aes(x=SNPs, y=Means)) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=Means-sd, ymax=Means+sd), width=.2,
                position=position_dodge(0.05))


## The function is currently defined as
function (filename,frac,ranking){
  N = nrow(filename$map)

  size = nrow(filename$fam)

  genotypes = as.numeric(filename$genotypes)
  genotypes[genotypes == 0] = NA
  genotypes = matrix(genotypes, nrow=size)
  genotypes = genotypes - 1

  phenotypes = filename$fam$affected == 2

  controls = !is.na(phenotypes) & !phenotypes
  cases    = !is.na(phenotypes) &  phenotypes

  controls = (1:size)[controls]
  cases    = (1:size)[cases]

  t(replicate(100,{
    controls_train = ceiling(frac * length(controls))
    controls_test  = length(controls) - controls_train
    cases_train    = ceiling(frac * length(cases))
    cases_test     = length(cases)    - cases_train

    training = c(sample(controls, size=controls_train),
                 sample(cases,    size=cases_train))
    training = sort(training)

    testing = setdiff(c(controls, cases), training)

    associations = snp.rhs.estimates(affected ~ 1,
                                     family="binomial",
                                     data=dat$fam[training, ],
                                     snp.data=dat$genotypes[training, ])

    beta    = sapply(associations, FUN=function (x) x$beta)


    if (ranking == "Z"){
      zscores = sapply(associations, FUN=function (x) x$beta/sqrt(x$Var.beta))
      ranking = order(zscores, decreasing=TRUE)
    }else if (ranking == "beta"){
      ranking = order(abs(beta), decreasing=TRUE)
    }else if (ranking == "iscore"){
      iscores = I_score(genotypes[training, ], phenotypes[training])
      ranking = order(iscores, decreasing=TRUE)
    }else if (ranking == "combined"){
      zscores = sapply(associations, FUN=function (x) x$beta/sqrt(x$Var.beta))
      iscores = I_score(genotypes[training, ], phenotypes[training])
      ranking = order(zscores + iscores, decreasing=TRUE)
    }

    rgenotypes = genotypes[testing, ranking]
    rbeta = beta[ranking]

    rbeta_sign = sign(rbeta)

    rgenotypes = t(t(rgenotypes) - colMeans(rgenotypes, na.rm=TRUE))
    rgenotypes[is.na(rgenotypes)] = 0

    polygenic = rowCumsums(t(t(rgenotypes) * rbeta))

    c(0.5, colAUC(polygenic, phenotypes[testing])[1, ])
  }))


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
