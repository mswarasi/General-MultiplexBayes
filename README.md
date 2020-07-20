# Multiplex-Group-Testing

This repository contains R programs of the article, "Estimating the prevalence of two or more diseases using outcomes from multiplex group testing." An R function "multDiseaseBayes" is provided to fit the proposed Bayesian estimation method and maximum a posteriori (MAP) estimation method, which can accommodate ANY group testing data involving multiplex assays and provide cost-effective estimates for the prevalence of multiple diseases as well as the assay accuracies (sensitivity and specificity). 

MainProgram.R - file that consists of the main R function combining the R subfunctions and FORTRAN subroutines. Documentation with illustrative (simulation) examples is also provided. 

SupportPrograms.txt - file that consists of standalone, executable R functions.

Note: For efficient execution, compute-intensive parts of the program are written in FORTRAN 95 and called from R through three DLL files, gbbstwodisgen.dll, mapacrtwodgen.dll, and ytiltwodbayes.dll, which work in a 64-bit R package. 


Reference

Warasi, M., McMahan, C., Tebbs, J., and Bilder, C. (2020+). Estimating the prevalence of two or more diseases using outcomes from multiplex group testing. Under review.
