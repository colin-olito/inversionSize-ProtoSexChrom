################################################################
#  Wright-Fisher forward simulations to validate approximation 
#  of the probability of invasion for a Y-LINKED inversions
#  capturing the sex-determining locus, and a second sexually
#  antagonistic locus
#
#  R code for W-F forward simulations. Generates output data
#  as .csv files saved to ./output/data/simResults/SA/
#
#
#  Author: Colin Olito
#
#  NOTES:  
#		



rm(list=ls())
#####################
##  Dependencies
source('R/functions-SA-Y-WFSims.R')
source('R/functions-SA-Y-determSims.R')


##############################################################################################
##############################################################################################
# Run sims until pcrit corresponding to invasion probability of 0.9997

######################################
##	Small population Size (N = 30k)	##
##	Gradient of sf values			##
##	no deleterious mutations		##
##	h = 1/2, sm = 0.01, r = 1/2		##
######################################
	sm     <-  0.01
	lower  <-  sm / (1 + sm)
	upper  <-  sm / (1 - sm)
	sfs    <-  seq((lower - 1*(upper - lower)),(upper + 1*(upper - lower)), length=15)
	par.list  <-  list(
					   gen    =  10000,
					   N      =  30000,
					   sm     =  sm,
					   sf     =  sfs,
					   hm     =  1/2,
					   hf     =  1/2,
					   r      =  1/2
					  )
runWFSimsYinversionValidatePinvApprox(par.list=par.list)




######################################
##	Small population Size (N = 30k)	##
##	Gradient of r values			##
##	no deleterious mutations		##
##	h = 1/2, s = 0.01				##
######################################
	rs  <-  seq(0,1/20,length=15)
	par.list  <-  list(
					   gen    =  10000,
					   N      =  30000,
					   sm     =  0.01,
					   sf     =  0.01,
					   hm     =  1/2,
					   hf     =  1/2,
					   r      =  rs
					  )
runWFSimsYinversionValidatePinvApprox(par.list=par.list)






