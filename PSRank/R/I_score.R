# Iscore

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
