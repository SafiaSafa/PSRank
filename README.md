# PSRank

Genome-wide association studies (GWASs) are designed to identify genetic
susceptibility factors for complex diseases. They are effective approaches for
identifying genetic polymorphisms associated with disease risk. Early GWASs of
major psychiatric disorders (schizophrenia, bipolar disorder, major depressive
disorder, autism spectrum disorders, and attention-deficit hyperactivity
disorder) revealed few hits that passed the significance thresholds. Very little
of the 80% of heritable variance can thus be captured (Huang et al., 2010).
This, among other things, has contributed to the so-called "missing
heritability" problem which states that unique genetic variations can not
explain the vast majority of the heritability of complex diseases. Several
explanations for missing heritability have been proposed (Wray et al., 2014)
ranging from a larger number of weaker effects polymorphisms, to rarer
polymorphisms (with possible greater effects) that would be poorly detected by
available genotyping methods (which focus on polymorphisms present in 5% or more
of the population).

Overall, much of the problem stems from a low statistical power. This depends on
four factors: (1) size of the cohort, (2) effect size distribution, (3)
significance level and (4) the statistic used for testing for an association.
The PSRank package concentrates on the fourth point, the metric, in order
statistical power.

For this, the population is split into two sets of data:

1. A learning set from which associations are sought using an additive model. An
effect size, a p-value and other summary statistics are calculated for each SNP.
SNPs can then be sorted by decreasing importance based on those summary
statistics, — traditionally, by increasing p-value or by decreasing (absolute)
effect size. Lo et al. (PNAS 2015) introduce an alternative statistic, the
I score: designed to recognize strongly predictive variables rather than
merely significant ones. The I score can be used as well by PSRank.
2. A test set from which the cumulative predictive power of SNPs can be assessed
using a version of the polygenic risk score developed by Purcell et al. using
AUC measures of accuracy.

To install the package `PSRank` and load it into R (version 3.4.2):

```
devtools::install_github("SafiaSafa/PSRank")
library(PSRank)
```
