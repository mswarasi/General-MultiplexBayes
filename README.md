# General-MultiplexBayes

This repository contains R programs of the article, "Estimating the prevalence of two or more diseases using outcomes from multiplex group testing." Two main R functions (mult.gt.bayes for L=1 assay and mult.gt.bayes_L2 for L=2 assays) are provided, each of which can implement the posterior sampling algorithm and the EM algorithm proposed in Warasi et al. (2022+). This article uses multiplex group testing data for estimating the coinfection probabilities and the assay accuracy probabilities (sensitivity and specificity).

R code of the simulation examples in the article is split into three files:

1. Simulation1.R -- to reproduce the estimation results shown in Table D.1 (in Web Appendix D).

2. Simulation2.R -- to reproduce the estimation results shown in Tables 2-3 (in the article).

3. Simulation3.R -- to reproduce the estimation results shown in Tables D.3-D.8 (in Web Appendix D). 



Reference

Warasi, M., Tebbs, J., McMahan, C., and Bilder, C. (2022+). Estimating the prevalence of two or more diseases using outcomes from multiplex group testing. Under review.


