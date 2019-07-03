################################################################
#  Wright-Fisher forward simulations to calculate Pr(fix | x) 
#  for Y-LINKED inversions capturing the sex-determining locus, 
#  and an arbitrary number of sexually antagonistic loci 
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
#source('R/functions-SA-Y-WFSims.R')
#source('R/functions-SA-Y-determSims.R')


##############################################################################################
##############################################################################################


#' Wright-Fisher simulations to estimate Pr(fix | x) for multi-locus SA model of a Y-linked 
#' inversion
#'
#' @title Wright-Fisher simulations to estimate Pr(fix | x) for multi-locus SA model of a Y-linked 
#' 		  inversion
#' 
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'					  N    =  1000,		# Population size
#'					  sm   =  0.1,		# selection coefficient in males
#'					  sf   =  0.1,		# selection coefficient in females
#'					  hm   =  1/2,		# dominance coefficient in males
#'					  hf   =  1/2,		# dominance coefficient in males
#'					  A    =  0.1,		# frequency of SA loci on the chromosome (Poisson rate parameter=)
#'					  P    =  0.1,		# Size of sl-PARl (fraction of chromosome)
#' 					  Ud   =  0.1,		# Deleterious mutation rate
#' 					  sd   =  0.05		# Selection coefficient for deleterious mutations.
#'				   	  )
#' @return Returns a list with timeseries for each genotype, 
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' detSimYInversion(par.list, Fx.init, Fy.init, eqStart = TRUE, threshold = 1e-6) 


simPrInvMultiLocusSA  <-  function(par.list)