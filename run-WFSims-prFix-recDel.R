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
#dataDirectoryExists  <-  dir.exists(".output/data")

#if(!dataDirectoryExists) {
#	dir.create("./data")
#}

#figuresDirectoryExists  <-  dir.exists("./output/figures")

#if(!figuresDirectoryExists) {
#	dir.create("./output/figures")
#}

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
source('R/functions-sheltering-recDel.R')

makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_deterministic_q")

makeDataPrFixInvSize(h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_deterministic_q")



#########################
## Beneficial Inversions
rm(list=ls())
source('R/functions-beneficial-recDel.R')

makeDataPrFixInvSize(sI = 0.02, h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_deterministic_q")

makeDataPrFixInvSize(sI = 0.02, h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^5), Nfname="_N100k_deterministic_q")


################
## SA Selection

rm(list=ls())
source('R/functions-SA-recDel.R')

makeDataPrFixInvSize(sf = 0.05, hf = 1/2, sm = 0.05, hm = 1/2, SAunlinked=TRUE,
					 h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname="_N10k_deterministic_q_r0.5")






#############################
#' Simulate invasion of SA using 
#' invasion fitness approximation
#' s_I
rm(list=ls())
source('R/functions-beneficial-recDel.R')

#' Evenly spaced log-scale points on recombination axis
#' 
log(0.002)
log(0.498)
logSeq  <-  seq(from=-6.214608, to=-0.6971552, len=5)
rSeq  <-  round(exp(logSeq), digits=3)
rSeq[4]  <-  rSeq[4] - 0.001



##' Additive SA fitness (h_f = h_m = 1/2)
## Equal selection

# Import sI data from Mathematica
sI.dat  <-  read.csv('./output/data/simResults/SA-sI-EqualSel_sf0.05_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
# subset sI data to these values
sI.subdat  <-  sI.dat[sI.dat$r %in% rSeq,]

for(i in 1:nrow(sI.subdat)) {
	fileid  <-  paste("SA_sIApprox_N10k_r", rSeq[i], "_sf", sI.subdat$sf[1], "_sm" , sI.subdat$sm[1], "_hf" , sI.subdat$hf[1], "_hm" , sI.subdat$hm[1], sep="")
	makeDataPrFixInvSize(sI = sI.subdat$sI[i], h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname=fileid)
}


#' Dominance Reversal SA fitness (h_f = h_m = 1/4)
# Import sI data from Mathematica
sI.dat  <-  read.csv('./output/data/simResults/SA-sI-EqualSel_sf0.05_sm0.05_hf0.25_hm0.25.csv', header=TRUE)
# subset sI data to these values
sI.subdat  <-  sI.dat[sI.dat$r %in% rSeq,]

for(i in 1:nrow(sI.subdat)) {
	fileid  <-  paste("SA_sIApprox_N10k_r", rSeq[i], "_sf", sI.subdat$sf[1], "_sm" , sI.subdat$sm[1], "_hf" , sI.subdat$hf[1], "_hm" , sI.subdat$hm[1], sep="")
	makeDataPrFixInvSize(sI = sI.subdat$sI[i], h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname=fileid)
}


# Sex biased selection
sI.dat  <-  read.csv('./output/data/simResults/SA-sI-BiasedSel_sf0.0526316_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
# subset sI data to these values
sI.subdat  <-  sI.dat[sI.dat$r %in% rSeq,]

for(i in 1:nrow(sI.subdat)) {
	fileid  <-  paste("SA_sIApprox_N10k_r", rSeq[i], "_sf", sI.subdat$sf[1], "_sm" , sI.subdat$sm[1], "_hf" , sI.subdat$hf[1], "_hm" , sI.subdat$hm[1], sep="")
	makeDataPrFixInvSize(sI = sI.subdat$sI[i], h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 nTot = 10^4, N.vals = c(10^4), Nfname=fileid)
}






#############################
#' Simulate invasion of SA using 
#' explicit recursions with equilibrium frequencies
#' s_I
rm(list=ls())
source('R/functions-SA-recDel.R')

#' Evenly spaced log-scale points on recombination axis
#' 
log(0.002)
log(0.498)
logSeq  <-  seq(from=-6.214608, to=-0.6971552, len=5)
rSeq  <-  round(exp(logSeq), digits=3)
rSeq[4]  <-  rSeq[4] - 0.001



##' Additive SA fitness (h_f = h_m = 1/2)
## Equal selection

# Import equilibrium frequency data from Mathematica
eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-EqualSel_sf0.05_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
# subset sI data to these values
eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]

for(i in 1:nrow(eqFreq.subdat)) {
	makeDataPrFixInvSize(sf = eqFreq.subdat$sf[i], hf = eqFreq.subdat$hf[i], sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], rSA=eqFreq.subdat$r[i], eq.df=eqFreq.subdat,
						 h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 	 nTot = 10^4, N.vals = c(10^4), Nfname="_NumEqFreq_N10k")
}

sI.vec  <-  c()
for(i in 1:nrow(eqFreq.subdat)) {
	sI.vec[i]  <-  sI_SA(sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], Yhat = eqFreq.subdat$Y[i], Xfhat = eqFreq.subdat$Xf[i])
}

##' Dominance Reversal SA fitness (h_f = h_m = 1/4)

# Import equilibrium frequency data from Mathematica
eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-EqualSel_sf0.05_sm0.05_hf0.25_hm0.25.csv', header=TRUE)
# subset sI data to these values
eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]

for(i in 1:nrow(eqFreq.subdat)) {

	makeDataPrFixInvSize(sf = eqFreq.subdat$sf[i], hf = eqFreq.subdat$hf[i], sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], rSA=eqFreq.subdat$r[i], eq.df=eqFreq.subdat,
						 h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 	 nTot = 10^4, N.vals = c(10^4), Nfname="_NumEqFreq_N10k")
}


##' Sex-biased selection
# Import equilibrium frequency data from Mathematica
eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-BiasedSel_sf0.0526316_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
# subset sI data to these values
eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]

for(i in 1:nrow(eqFreq.subdat)) {
#	fileid  <-  paste("eqFreq_N10k_r", rSeq[i], "_sf", eqFreq.subdat$sf[i], "_sm" , eqFreq.subdat$sm[i], "_hf" , eqFreq.subdat$hf[i], "_hm" , eqFreq.subdat$hm[i], sep="")
	makeDataPrFixInvSize(sf = eqFreq.subdat$sf[i], hf = eqFreq.subdat$hf[i], sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], rSA=eqFreq.subdat$r[i], eq.df=eqFreq.subdat,
						 h = 0.25, s = 0.01, Us.factor.vals = c(2, 5, 10),
					 	 nTot = 10^4, N.vals = c(10^4), Nfname="_NumEqFreq_N10k")
}
sI.vec  <-  c()
for(i in 1:nrow(eqFreq.subdat)) {
	sI.vec[i]  <-  sI_SA(sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], Yhat = eqFreq.subdat$Y[i], Xfhat = eqFreq.subdat$Xf[i])
}



# >> CONTINUE RUNNING SIMULATIONS FROM HERE


################################################################
################################################################
##  End of new sim code for partially recessive del. mut.
################################################################
################################################################



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
