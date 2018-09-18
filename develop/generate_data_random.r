# Initialisation ----

library("snpStats")

set.seed(1)

# Parameters ----

size = 1000  # Size of cohorts
N    = 50000 # Number of SNPs (all independent)

frequency = runif(N, min=0.05, max=.95) # MAF is > 5%
weight    = rexp(N) * sample(c(-1, 1), size=N, replace=TRUE)

# Genotypes ----

mat = rbinom(size * N, size=2, prob=frequency)
mat = matrix(mat, nrow=size, byrow=TRUE)

# Phenotypes ----

liability = colSums(t(mat) * weight)
liability = liability + rnorm(size, mean=0, sd=sd(liability))

threshold = quantile(liability, 2/3) # Prevalence is 33%

phenotype = liability > threshold

# Add missing values (1%) ----

mat[sample(c(T, F), size=N * size, replace=TRUE, prob=c(1, 99))] = NA
phenotype[sample(c(T, F), size=size, replace=TRUE, prob=c(1, 99))] = NA

# Conversion to PLINK format ----

snp_names = paste("rs", 1:N, sep="")
ind_names = paste("ID", 1:size, sep="")
colnames(mat) = snp_names
rownames(mat) = ind_names

genotypes = as(mat, "SnpMatrix")

phenotype = ifelse(is.na(phenotype), 0, phenotype + 1)

subject_data = data.frame(FID=ind_names, IID=ind_names,
                          sex=c(1, 2), phenotype=phenotype)

rownames(subject_data) = ind_names

snp_data = data.frame(chromosome=0, position=1:N, A="A", B="G")

rownames(snp_data) = snp_names

write.plink("output/unif_Laplace_omni",
            snps=genotypes,
            subject.data=subject_data,
            phenotype=phenotype,
            sex=sex,
            snp.data=snp_data,
            chromosome=chromosome,
            allele.1=A,
            allele.2=B,
            position=position)
