#  Estimate Genome Size from k-mer distribution data: 
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


source('R/functions-estGenomeSize.R')


##########################
# Estimate Genome Size 
# from k-mer distribution
##########################

