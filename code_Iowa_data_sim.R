
## The code below simulates individual test true 
## statuses (1 = positive, 0 = negative) at the  
## coinfection prevalence vectors shown below. 
## These vectors are obtained from individual case
## identification outcomes in the Iowa CT/NG data.  

set.seed( 123 )

p.urine <- c( 0.9082, 0.0811, 0.0057, 0.0050 )
p.swab <-  c( 0.9086, 0.0812, 0.0054, 0.0048 )

## Urine stratum:
N.urine <- 4402   # number of urine specimens tested 
Ymul <- rmultinom(N.urine, 1, p.urine)
Ytil1 <- Ymul[2, ] + Ymul[4, ]
Ytil2 <- Ymul[3, ] + Ymul[4, ]
Yt.urine <- cbind( Ytil1, Ytil2 )
saveRDS(Yt.urine, "./data_analysis/Iowa_data/data_urine_stratum.rds ")

## Swab stratum:
N.swab <- 10048   # number of swab specimens tested 
Ymul <- rmultinom(N.swab, 1, p.swab)
Ytil1 <- Ymul[2, ] + Ymul[4, ]
Ytil2 <- Ymul[3, ] + Ymul[4, ]
Yt.swab <- cbind( Ytil1, Ytil2 )
saveRDS(Yt.swab, "./data_analysis/Iowa_data/data_swab_stratum.rds ")

