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
  
  as.list(c(0.5, colAUC(polygenic, phenotypes[testing])[1, ]))
  
}




fam = "output/unif_Laplace_omni.fam"
bim = "output/unif_Laplace_omni.bim"
bed = "output/unif_Laplace_omni.bed"
dat = read.plink(bed, bim, fam)

frac = .5


mat = NULL
for (i in 1:100){
  if (length(mat) == 0){
    print(i)
    mat = PSRank(dat,
                 frac=.5,
                 ranking="beta")
  }else{
    print(i)
    mat = rbind(mat, PSRank(dat,
                   frac=.5,
                   ranking="beta"))
  }
}

ggplot(mat, aes(y = AUC, x = SNPs)) + geom_smooth() +
  geom_point() + theme_bw()