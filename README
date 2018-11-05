# PSRank


Genome-wide association studies (GWASs) are designed to identify genetic susceptibility factors for complex diseases.
They are effective approaches for identifying genetic polymorphisms associated with disease risk. 
Early GWASs of major psychiatric disorders (schizophrenia, bipolar disorder, major depressive
disorder, autism spectrum disorders, and attention-deficit hyperactivity disorder) revealed
few "hits" that passed the significance thresholds. We are still far from the hoped-for 80% of
explained variance due to the genetic heritability expected (Huang et al., 2010). This has
contributed to the so-called "Missing heritability" problem which states that unique genetic
variations can not explain the vast majority of the heritability of complex diseases. Several
explanations for missing heritability have been proposed (Wray et al., 2014) ranging from
a larger number of weaker effects polymorphisms, to rarer polymorphisms (with possible
greater effects) that would be poorly detected by available genotyping methods (which focus
on polymorphisms present in 5% or more of the population) .
Overall, it boils down to a relatively low statistical power as heist. Which depends on three
factors: 
1. The size of the cohort. 
2. The size of the effect. 
3. The metric used.

THe `PSRank`'s package concentrate on the metric used. To recapture missing signal it creates a polygenic's model. 
For this, the population is separate into two sets of data: 
1. A learning set where will be estimated associations through linear univariate model. This set of SNPs is going to be sorted. 
The packages allows to sorte the SNPs according differentes scores. Classically, by *P-values* or *β*.
The Lo et al., 2015 ’s article offer an alternative statistic to the χ2, the *I-score*: 
designed to recognize strongly predictive variables rather than highly significant. This score can be used to sort a set of SNPs.
2. And a test set where will be calculated a cumulative score, called "polygenic risk score" developed by Purcell et al. and
on which we will evaluate the predictivity

The data used is produce by simulation according to a Laplace distribution.

To install the package `PSRank` and load it into R :

```
devtools::install_github("SafiaSafa/PSRank")
library(PSRank)
```
