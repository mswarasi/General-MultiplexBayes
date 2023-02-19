
## The tables shown below have been reported 
## in the Supporting Information (Appendix D):

Table_D1 <- readRDS("./simulation/intermediate_results/Table_D1.rds")
Table_D2 <- readRDS("./simulation/intermediate_results/Table_D2.rds")
Table_D3 <- readRDS("./simulation/intermediate_results/Table_D3.rds")
Table_D4 <- readRDS("./simulation/intermediate_results/Table_D4.rds")
Table_D5 <- readRDS("./simulation/intermediate_results/Table_D5.rds")
Table_D6 <- readRDS("./simulation/intermediate_results/Table_D6.rds")
Table_D7 <- readRDS("./simulation/intermediate_results/Table_D7.rds")
Table_D8 <- readRDS("./simulation/intermediate_results/Table_D8.rds")
Table_D9 <- readRDS("./simulation/intermediate_results/Table_D9.rds")
Table_D10 <- readRDS("./simulation/intermediate_results/Table_D10.rds")
Table_D11 <- readRDS("./simulation/intermediate_results/Table_D11.rds")


## Table D.1
protocols <- paste("                                 MPT                     H2                    H3                    H4                    AT")
cline <-     paste("                          ===================    ===================   ===================   ===================   ===================")

noquote(protocols); noquote(cline); Table_D1


## Table D.2
config <- paste("   Configuration I            Configuration II   ")
cline <- paste("=========================  ==========================")

noquote(config); noquote(cline); Table_D2


## Table D.3
protocols <- paste("                                   H2                    H3                    H4                   AT")
cline <-     paste("                          ===================   ===================   ===================   ====================")

noquote(protocols); noquote(cline); Table_D3


## Table D.4
protocols <- paste("                         H2                  H3                 H4                  AT")
cline <-     paste("                 =================   =================   =================   =================")

noquote(protocols); noquote(cline); Table_D4


## Table D.5
protocols <- paste("                         H2                  H3                 H4                  AT")
cline <-     paste("                 =================   =================   =================   =================")

noquote(protocols); noquote(cline); Table_D5


## Table D.6
protocols <- paste("                                   H2                    H3                    H4                   AT")
cline <-     paste("                          ===================   ===================   ===================   ====================")

noquote(protocols); noquote(cline); Table_D6


## Table D.7
protocols <- paste("                         H2                  H3                 H4                  AT")
cline <-     paste("                 =================   =================   =================   =================")

noquote(protocols); noquote(cline); Table_D7


## Table D.8
protocols <- paste("                         H2                  H3                 H4                  AT")
cline <-     paste("                 =================   =================   =================   =================")

noquote(protocols); noquote(cline); Table_D8


## Table D.9
protocols <- paste("                                    MPT                     H2                    H3                    H4                    AT")
cline <-     paste("                            ===================    ===================   ===================   ===================   ====================")

noquote(protocols); noquote(cline); Table_D9


## Table D.10
protocols <- paste("                                     H2                    H3                    H4                   AT")
cline <-     paste("                            ===================   ===================   ===================   ====================")

noquote(protocols); noquote(cline); Table_D10


## Table D.11
protocols <- paste("                                       H2                  H3                 H4                  AT")
cline <-     paste("                               =================   =================   =================   =================")

noquote(protocols); noquote(cline); Table_D11

