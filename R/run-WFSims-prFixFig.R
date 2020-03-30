################################################################
#  RUN W-F SIMULATIONS FOR FIG.1 -- Pr(fix | x)
#
#  R code for W-F forward simulations. Generates output data
#  as .csv files saved to ./output/data/simResults/SA/
#
#
#  Author: Colin Olito
#
#  NOTES:  
#		




#####################
##  Dependencies
rm(list=ls())
source('R/functions-Neutral-WFSims.R')



############
# Y -linked
############


############
# Neutral

N_y = 100
par.list  <-  list(Ny   = 10^2,
				   sd   = 0.05,
				   Ud   = 0.2)
neutralFwdSimY(par.list=par.list)

# N_y = 500
par.list  <-  list(Ny   = 5*10^2,
				   sd   = 0.05,
				   Ud   = 0.2)
neutralFwdSimY(par.list=par.list)

# N_y = 750
par.list  <-  list(Ny   = 7.5*10^2,
				   sd   = 0.05,
				   Ud   = 0.2)
neutralFwdSimY(par.list=par.list)

# N_y = 1000
par.list  <-  list(Ny   = 10^3,
				   sd   = 0.05,
				   Ud   = 0.2)
neutralFwdSimY(par.list=par.list)

# N_y = 10000
par.list  <-  list(Ny   = 10^4,
				   sd   = 0.05,
				   Ud   = 0.2)
neutralFwdSimY(par.list=par.list)



############
# Beneficial
rm(list=ls())
source('R/functions-Beneficial-WFSims.R')

# N_y = 1000
par.list  <-  list(Ny   = 10^3,
				   sI   = 0.02,
				   sd   = 0.02,
				   Ud   = 0.2)
beneficialFwdSimY(par.list=par.list)

# N_y = 1000
par.list  <-  list(Ny   = 10^3,
				   sI   = 0.02,
				   sd   = 0.03,
				   Ud   = 0.2)
beneficialFwdSimY(par.list=par.list)

# N_y = 10000
par.list  <-  list(Ny   = 10^4,
				   sI   = 0.02,
				   sd   = 0.02,
				   Ud   = 0.2)
beneficialFwdSimY(par.list=par.list)



############
# SA Selection
rm(list=ls())
source('R/functions-SA-Y-2Locus-determSims.R')
source('R/functions-SA-Y-2Locus-WFSims.R')

Ud          <-  0.1
sd          <-  0.01
hf          <-  1/2
sf          <-  0.05
hm          <-  1/2
sm          <-  0.05
r_linked    <-  0.0101
r_unlinked  <-  1/2
A           <-  1
P           <-  0.05

## SA locus initially unlinked with SDL 
# N_y = 10^3
 par.list  <-  list(
					gen     =  10^4,
					N       =  10^4,
					sm      =  sm,
					sf      =  sf,
					hm      =  hm,
					hf      =  hf,
					r       =  r_unlinked,
					Ud      =  Ud,
					sd      =  sd,
					x       =  NA,
					selType =  "SA",
					selPars =  c(hf, sf, hm, sm),
					A       =  A,
					P       =  P
					)
SAFwdSimYpFixFig(par.list=par.list)




## SA locus initially LINKED with SDL 
# N_y = 10^3
 par.list  <-  list(
					gen     =  10^4,
					N       =  10^4,
					sm      =  sm,
					sf      =  sf,
					hm      =  hm,
					hf      =  hf,
					r       =  r_linked,
					Ud      =  Ud,
					sd      =  sd,
					x       =  NA,
					selType =  "SA",
					selPars =  c(hf, sf, hm, sm),
					A       =  A,
					P       =  P
					)
SAFwdSimYpFixFig(par.list=par.list)
