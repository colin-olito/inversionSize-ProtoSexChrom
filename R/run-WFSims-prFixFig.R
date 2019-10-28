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

# N_y = 10
par.list  <-  list(Ny   = 1*10^1,
				   sd   = 0.05,
				   Ud   = 0.1)
neutralFwdSimY(par.list=par.list)

# N_y = 100
par.list  <-  list(Ny   = 1*10^2,
				   sd   = 0.05,
				   Ud   = 0.1)
neutralFwdSimY(par.list=par.list)

# N_y = 1000
par.list  <-  list(Ny   = 1*10^3,
				   sd   = 0.05,
				   Ud   = 0.1)
neutralFwdSimY(par.list=par.list)



############
# Beneficial
rm(list=ls())
source('R/functions-Beneficial-WFSims.R')

# N_y = 10
par.list  <-  list(Ny   = 1*10^1,
				   sI   = 0.02,
				   sd   = 0.05,
				   Ud   = 0.1)
beneficialFwdSimY(par.list=par.list)

# N_y = 50
par.list  <-  list(Ny   = 5*10^1,
				   sI   = 0.02,
				   sd   = 0.05,
				   Ud   = 0.1)
beneficialFwdSimY(par.list=par.list)

# N_y = 100
par.list  <-  list(Ny   = 1*10^2,
				   sI   = 0.02,
				   sd   = 0.05,
				   Ud   = 0.1)
beneficialFwdSimY(par.list=par.list)

# N_y = 1000
par.list  <-  list(Ny   = 1*10^3,
				   sI   = 0.02,
				   sd   = 0.05,
				   Ud   = 0.1)
beneficialFwdSimY(par.list=par.list)
