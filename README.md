# General-Framework-for-Prevalence-Estimation

This repository provides R programs (and associated FORTRAN DLL's) of the article titled, "Estimating the prevalence of two or more diseases using outcomes from multiplex group testing." A user-friendly R function "multDiseaseBayes" is provided to fit the proposed estimation methods, Bayesian estimation and maximum a posteriori (MAP) estimation, for any group testing algorithms involving multiplex assays. For efficient execution, compute-intensive parts of the program are written in FORTRAN and called from R through three DLL files, gbbstwodisgen.dll, mapacrtwodgen.dll, and ytiltwodbayes.dll, which work in a 64-bit R package. The R function with documentation and illustrative examples is  provided in the file: MainProgram.R

MainProgram.R - file that consists of the main R function combining the R subfunctions and FORTRAN subroutines.
SupportPrograms.txt - file that consists of standalone, executable R functions.

Note: The DLL files, gbbstwodisgen.dll, mapacrtwodgen.dll, and ytiltwodbayes.dll, are combiled subroutines written in a FORTRAN 95 compiler. 


