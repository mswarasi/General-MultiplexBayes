
## The tables shown below have been reported in the manuscript:

Table_4 <- readRDS("./data_analysis/intermediate_results/Table_4.rds")
Table_5 <- readRDS("./data_analysis/intermediate_results/Table_5.rds")
Table_6 <- readRDS("./data_analysis/intermediate_results/Table_6.rds")

## Table 4
Table_4


## Table 5
protocols <- paste("                          H2                 H3                  AT")
cline <-     paste("                   ============    ===================   ===================")

noquote(protocols); noquote(cline); Table_5


## Table 6
protocols <- paste("                               H2            H3            AT")
cline <-     paste("                         ============   ===========   ===========")

noquote(protocols); noquote(cline); Table_6

