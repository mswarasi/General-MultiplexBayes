
## The tables shown below have been reported in the manuscript:

Table_1 <- readRDS("./simulation/intermediate_results/Table_1.rds")
Table_2 <- readRDS("./simulation/intermediate_results/Table_2.rds")
Table_3 <- readRDS("./simulation/intermediate_results/Table_3.rds")


## Table 1
config <- paste("   Configuration I            Configuration II   ")
cline <- paste("=========================  ==========================")
noquote(config); noquote(cline); Table_1


## Table 2
protocols <- paste("                                  H2                       H3                           H4                         AT")
cline <-     paste("                         ===================    =========================   ==========================   ===========================")

noquote(protocols); noquote(cline); Table_2



## Table 3
protocols <- paste("                                   H2                  H3                  H4                  AT")
cline <-     paste("                           =================   =================   ==================   =================")

noquote(protocols); noquote(cline); Table_3

