I_score = function (genotypes, phenotypes) {
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

PSRank = function (filename,frac,ranking){
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

# Progress bar 
#install.packages("svMisc")
#require(svMisc)
#for (i in 0:101){
#  progress(i, progress.bar=T)
#  Sys.sleep(0.01)
#  if(i==101) cat("Done!\n")
#}
