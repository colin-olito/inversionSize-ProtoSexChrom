################################################################
#  RUN SIMULATIONS and CREATE OUTPUT DATA FOR PLOTTING
#  
#  R code for forward simulations. Generates output data
#  as .csv files saved to ./output/data
#
#
#  Author: Colin Olito
#
#  NOTES:  
#		

######################################
#' Create output directories if they
#' do not already exist
dataDirectoryExists  <-  dir.exists("./output/data/simResults")

if(!dataDirectoryExists) {
	dir.create("./output/data/simResults")
}

figuresDirectoryExists  <-  dir.exists("./output/figures")

if(!figuresDirectoryExists) {
	dir.create("./output/figures")
}

########################################
# Fixation Probability ~ Inv. Size
########################################
#' Note: these simulations create data
#' to produce Fig. X showing the fixation 
#' probability of different sized neutral 
#' inversions expanding the SLR on Y 
#' chromosomes under partially recessive
#' deleterious mutation pressure. 
#' Uses multilocus recursions.

###################################
## Sheltering (neutral inversions)
rm(list=ls())
source('R/functions-Neutral-WFSims-CORRECTION.R')

makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
					 nTot = 10^4, N = 10^4, Nfname="_N10k")

makeDataPrFixInvSizeHaploid(h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
					 nTot = 10^4, N = 10^5, Nfname="_N100k")



#########################
## Beneficial Inversions
rm(list=ls())
#source('R/functions-beneficial-recDel.R')
source('R/functions-Beneficial-WFSims-CORRECTION.R')



# makeDataPrFixInvSize(sI = 0.02, h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
# 					 nTot = 10^4, N.vals = c(10^4))
makeDataBeneficialPrFixInvSizeHaploid(sI = 0.02, h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
 					 nTot = 10^4, N = 10^4, Nfname="_N10k")

#makeDataPrFixInvSize(sI = 0.02, h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
#					 nTot = 10^4, N.vals = c(10^5))
makeDataBeneficialPrFixInvSizeHaploid(sI = 0.02, h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
 					 nTot = 10^4, N = 10^5, Nfname="_N100k")


################
## SA Selection

rm(list=ls())
source('R/functions-SA-recDel.R')

makeDataPrFixInvSize(sf = 0.05, hf = 1/2, sm = 0.05, hm = 1/2, rSA=1/2,
					 h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = 10^4, Nfname="_N10k_deterministic_q_r0.5")

makeDataPrFixInvSize(sf = 0.05, hf = 1/2, sm = 0.05, hm = 1/2, rSA=1/2,
					 h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = 10^5, Nfname="_N100k_deterministic_q_r0.5")






#############################
#' Simulate invasion of SA using 
#' invasion fitness approximation
#' s_I
rm(list=ls())
source('R/functions-Beneficial-WFSims-CORRECTION.R')

#' Evenly spaced log-scale points on recombination axis
#' 
log(0.002)
log(0.498)
logSeq   <-  seq(from=-6.214608, to=-0.6971552, len=5)
rSeq     <-  round(exp(logSeq), digits=3)
rSeq[4]  <-  rSeq[4] - 0.001



##' Additive SA fitness (h_f = h_m = 1/2)
## Equal selection

# Import sI data from Mathematica
sI.dat  <-  read.csv('./output/data/simResults/SA-sI-EqualSel_sf0.05_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
# subset sI data to these values
sI.subdat  <-  sI.dat[sI.dat$r %in% rSeq,]
sI.subdat  <-  sI.subdat[sI.subdat$sf == 0.05,]

for(i in 1:nrow(sI.subdat)) {
	fileid  <-  paste("_SA_sIApprox_N100k_r", rSeq[i], "_sf", sI.subdat$sf[1], "_sm" , sI.subdat$sm[1], "_hf" , sI.subdat$hf[1], "_hm" , sI.subdat$hm[1], sep="")
	makeDataBeneficialPrFixInvSizeHaploid(sI = sI.subdat$sI[i], h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
										  nTot = 10^4, N = 10^5, Nfname=fileid)
}


#' Dominance Reversal SA fitness (h_f = h_m = 1/4)
# Import sI data from Mathematica
sI.dat  <-  read.csv('./output/data/simResults/SA-sI-EqualSel_sf0.05_sm0.05_hf0.25_hm0.25.csv', header=TRUE)
# subset sI data to these values
sI.subdat  <-  sI.dat[sI.dat$r %in% rSeq,]

for(i in 1:nrow(sI.subdat)) {
	fileid  <-  paste("_SA_sIApprox_N100k_r", rSeq[i], "_sf", sI.subdat$sf[1], "_sm" , sI.subdat$sm[1], "_hf" , sI.subdat$hf[1], "_hm" , sI.subdat$hm[1], sep="")
	makeDataBeneficialPrFixInvSizeHaploid(sI = sI.subdat$sI[i], h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
					 nTot = 10^4, N = 10^5, Nfname=fileid)
}


# Sex biased selection
sI.dat  <-  read.csv('./output/data/simResults/SA-sI-BiasedSel_sf0.0526315789473684_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
# subset sI data to these values
sI.subdat  <-  sI.dat[sI.dat$r %in% rSeq,]

for(i in 1:nrow(sI.subdat)) {
	fileid  <-  paste("_SA_sIApprox_N100k_r", rSeq[i], "_sf", sI.subdat$sf[1], "_sm" , sI.subdat$sm[1], "_hf" , sI.subdat$hf[1], "_hm" , sI.subdat$hm[1], sep="")
	makeDataBeneficialPrFixInvSizeHaploid(sI = sI.subdat$sI[i], h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
					 nTot = 10^4, N = 10^5, Nfname=fileid)
}







################################################################
################################################################
##  WARNING: All code after this is functional but deprecated!
##  		 It will produce results from a much earlier 
##			 version of the manuscript which focused on 
##			 additive deleterious mutations.
################################################################
################################################################



#############################
#' Simulate invasion of SA using 
#' explicit recursions with equilibrium frequencies
#' s_I

# rm(list=ls())
# source('R/functions-SA-recDel.R')

#' Evenly spaced log-scale points on recombination axis
#' 
# log(0.002)
# log(0.498)
# logSeq   <-  seq(from=-6.214608, to=-0.6971552, len=5)
# rSeq     <-  round(exp(logSeq), digits=3)
# rSeq[4]  <-  rSeq[4] - 0.001



##' Additive SA fitness (h_f = h_m = 1/2)
## Equal selection

# Import equilibrium frequency data from Mathematica
# eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-EqualSel_sf0.05_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
# subset sI data to these values
# eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]

# for(i in 1:nrow(eqFreq.subdat)) {
# 	makeDataPrFixInvSize(sf = eqFreq.subdat$sf[i], hf = eqFreq.subdat$hf[i], sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], rSA=eqFreq.subdat$r[i], eq.df=eqFreq.subdat,
# 						 h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
# 					 	 nTot = 10^4, N.vals = c(10^5), Nfname="_NumEqFreq_N100k")
# }

# sI.vec  <-  c()
# for(i in 1:nrow(eqFreq.subdat)) {
# 	sI.vec[i]  <-  sI_SA(sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], Yhat = eqFreq.subdat$Y[i], Xfhat = eqFreq.subdat$Xf[i])
# }

##' Dominance Reversal SA fitness (h_f = h_m = 1/4)

# Import equilibrium frequency data from Mathematica
# eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-EqualSel_sf0.05_sm0.05_hf0.25_hm0.25.csv', header=TRUE)
# subset sI data to these values
# eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]

# for(i in 1:nrow(eqFreq.subdat)) {

# 	makeDataPrFixInvSize(sf = eqFreq.subdat$sf[i], hf = eqFreq.subdat$hf[i], sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], rSA=eqFreq.subdat$r[i], eq.df=eqFreq.subdat,
# 						 h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
# 					 	 nTot = 10^4, N.vals = c(10^4), Nfname="_NumEqFreq_N10k")
# }


##' Sex-biased selection
# Import equilibrium frequency data from Mathematica
# eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-BiasedSel_sf0.0526316_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
# subset sI data to these values
# eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]

# for(i in 1:nrow(eqFreq.subdat)) {
# #	fileid  <-  paste("eqFreq_N10k_r", rSeq[i], "_sf", eqFreq.subdat$sf[i], "_sm" , eqFreq.subdat$sm[i], "_hf" , eqFreq.subdat$hf[i], "_hm" , eqFreq.subdat$hm[i], sep="")
# 	makeDataPrFixInvSize(sf = eqFreq.subdat$sf[i], hf = eqFreq.subdat$hf[i], sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], rSA=eqFreq.subdat$r[i], eq.df=eqFreq.subdat,
# 						 h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
# 					 	 nTot = 10^4, N.vals = c(10^4), Nfname="_NumEqFreq_N10k")
# }
# sI.vec  <-  c()
# for(i in 1:nrow(eqFreq.subdat)) {
# 	sI.vec[i]  <-  sI_SA(sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], Yhat = eqFreq.subdat$Y[i], Xfhat = eqFreq.subdat$Xf[i])
# }

