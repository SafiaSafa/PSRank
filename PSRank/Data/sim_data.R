library("snpStats")

params = list(seed=1,           # Random seed
              size=2000,        # Size of cohorts
              n=1000,           # Number of causal SNPs
              recycle=TRUE,     # Only one cohort
              AUC=FALSE,        # Correlation with liability vs AUC with phenotype
              prefix='sim9d',   # What prefix to use
              replicates=1000,  # Number of replicates
              distr='rlaplace') # Effect size distribution

set.seed(params$seed)

rlaplace = function (n, ...) {
  rexp(n, ...) * sample(c(-1, 1), size=n, replace=TRUE)
}

m = 10000
n = params$n
f = runif(m, min=.05, max=.95)
# w = rlaplace(n)
w = eval(parse(text=params$distr))(n)
prevalence = .25

h2 = .5

rs = c(paste("c", 1:n, sep=""), paste("d", 1:(m - n), sep=""))

size = 10000

genotypes = rbinom(size * n, size=2, prob=f[1:n])
genotypes = matrix(genotypes, nrow=size, byrow=TRUE)
liabilities = colSums(t(genotypes) * w)

sigma = sd(liabilities) * sqrt((1 - h2)/h2)
noise = rnorm(size, sd=sigma)
liabilities = liabilities + noise
thrs = quantile(liabilities, 1 - prevalence, names=FALSE)

size = params$size
ID = paste("ID", 1:size, sep="")

genotypes = rbinom(size * m, size=2, prob=f)
genotypes = matrix(genotypes, nrow=size, byrow=TRUE)
liabilities = colSums(t(genotypes[, 1:n]) * w) + rnorm(size, sd=sigma)
phenotypes = liabilities > thrs

rownames(genotypes) = ID
colnames(genotypes) = rs

data = data.frame(affected=phenotypes)
rownames(data) = ID

subject_data = data.frame(FID=ID, IID=ID,
                          sex=c(1, 2), phenotype=phenotypes)

rownames(subject_data) = ID

snp_data = data.frame(chromosome=1, position=1:m, A="A", B="G")

rownames(snp_data) = rs

write.plink("sim9d",
            snps=as(genotypes, "SnpMatrix"),
            subject.data=subject_data,
            phenotype=phenotype,
            sex=sex,
            snp.data=snp_data,
            chromosome=chromosome,
            allele.1=A,
            allele.2=B,
            position=position)
