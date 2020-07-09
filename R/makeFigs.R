#  Functions to generate figures for: 
#    
#  Title: 	Local adaptation and the evolution of  
#         	X-linked vs. autosomal inversions 
#
#			A Contributed paper for Phil. Trans.
#			Roy. Soc. Theme Issue put together by 
#		   	the ESEB Special Topics Network: linking 
#			local adaptation with the evolution of 
#			sex-differences.
#
#
#  Author: Colin Olito
#
#
#  NOTES: Run this file, either from terminal using Rscript,
#		  or interactively in R. This should create all the 
#		  figures needed to correctly compile the mansucript
#		  LaTeX file.  
#          

rm(list=ls())
###############
# Dependencies
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)
library(data.table)

source('R/functions-figures.R')


########################
# Figures for the paper
########################


 toPdf(fixationProbabilityFigure(), figPath(name='Fig1.pdf'), width=10, height=4)
 embed_fonts(figPath(name='Fig1.pdf'))

# Pr(fix | x, SDR)
 toPdf(PrCatchFixFigure(SDRloc = c(1/2, 1/10)), figPath(name='PrCatchFixFig.pdf'), width=10, height=8)
 embed_fonts(figPath(name='PrCatchFixFig.pdf'))

 toPdf(expectedDistributionFig(SDRloc = c(1/2, 1/10)), figPath(name='expDistFig.pdf'), width=7, height=7)
 embed_fonts(figPath(name='expDistFig.pdf'))

source('R/functions-figures.R')
 toPdf(determInvFreqPlot(wesMovie = 'Zissou1'), figPath(name='detEqInvFreqFig.pdf'), width=10, height=10)
 embed_fonts(figPath(name='detEqInvFreqFig.pdf'))

source('R/functions-figures.R')
 toPdf(recombEffectSATwoLocus(), figPath(name='recombEffect.pdf'), width=10, height=6)
 embed_fonts(figPath(name='recombEffect.pdf'))


# toPdf(Fig3Alt(), 
#             figPath(name='Fig3Alt.pdf'), width=7, height=7)
# embed_fonts(figPath(name='Fig3Alt.pdf'))

# toPdf(Fig4(), 
#             figPath(name='Fig4.pdf'), width=7, height=7)
# embed_fonts(figPath(name='Fig4.pdf'))

# toPdf(propEstSuppFigs(), 
#             figPath(name='ProportionEstablishSuppFig.pdf'), width=18, height=10)
# embed_fonts(figPath(name='ProportionEstablishSuppFig.pdf'))

# toPdf(finalFreqSuppFig(), 
#             figPath(name='FinalInvFreqSuppFig.pdf'), width=7, height=7)
# embed_fonts(figPath(name='FinalInvFreqSuppFig.pdf'))

########################
# Supplementary Figures
########################






######################
# Exploratory Figures
######################

# Validation of analytic approximations for the SA model 
# invasion probability of Y-linked inversions
toPdf(SA_Y_validatePinv(s.df = "SA_Y_validatePinv_sfGrad_N30000_hf0.5_hm0.5_sm0.01_r0.5.csv",
				  r.df = "SA_Y_validatePinv_rGrad_N30000_hf0.5_sf0.01_hm0.5_sm0.01.csv"), 
      figPath(name='SA_Y_validatePinv.pdf'), width=10, height=5)
embed_fonts(figPath(name='SA_Y_validatePinv.pdf'))
