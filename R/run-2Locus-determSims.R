################################################################
#  RUN Deterministic Simulations for FIG.2
#
#  R code for deterministic forward simulations. Generates 
#  output data as .csv files saved to ./output/data/simResults/SA/
#
#
#  Author: Colin Olito
#
#  NOTES:  
#		




#####################
##  Dependencies
rm(list=ls())


############
# Y -linked
############
source('R/functions-SA-Y-2Locus-determSims.R')

# r = 1/2
detYSAInversionEQFreqs(selRange = c(0,1), by=0.01, generations=5000,
					   dom = c(1/2, 1/2), r = 1/2)
# r = 0.02
detYSAInversionEQFreqs(selRange = c(0,1), by=0.01, generations=10000, 
					   dom = c(1/2, 1/2), r = 0.02)
# r = 0.002
detYSAInversionEQFreqs(selRange = c(0,1), by=0.01, generations=10000, 
					   dom = c(1/2, 1/2), r = 0.002)



############
# X-linked
############
rm(list=ls())
source('R/functions-SA-X-2Locus-determSims.R')

# r = 1/2
detXSAInversionEQFreqs(selRange = c(0,1), by=0.01, generations = 30000,
					   dom = c(1/2, 1/2), r = 1/2, eq.threshold = 1e-7)
# r = 0.02
detXSAInversionEQFreqs(selRange = c(0,1), by=0.01, generations = 20000,
					   dom = c(1/2, 1/2), r = 0.02, eq.threshold = 1e-7)
# r = 0.002
detXSAInversionEQFreqs(selRange = c(0,1), by=0.01, generations = 10000,
					   dom = c(1/2, 1/2), r = 0.002, eq.threshold = 1e-7)



selRange = c(0.02,1)
by=0.02
generations = 10000
dom = c(1/2, 1/2)
r = 0.5
eq.threshold = 1e-7
sfs      <-  seq(from=selRange[1], to=selRange[2], by=by)
hf=dom[1]
hm=dom[2]

sf  <-  sfs[6]
sm  <-  sfs[11]

sf  <-  sfs[30]
sm  <-  sfs[33]
 par.list  <-  list(
				   selType  =  "SA",
				   selPars  =  c(hf, sf, hm, sm),
				   gen      =  10000,
				   r        =  r
				   )

test  <-  detSimXInversion(par.list=par.list)
test$eqInvFreq
test$gen
str(test)
initInvFreq  <-  (2*sum(test$Fx.gen[1,      c(2,3,5)]*c(1/2, 1 ,1/2)) + sum(test$Fy.gen[1,      c(2,4)]))/3
invFreqDyn  <-  (2*((test$Fx.gen[,2]*1/2) + (test$Fx.gen[,3]) + (test$Fx.gen[,5]*1/2)) + (test$Fy.gen[,2] + test$Fy.gen[,4]))/3
initInvFreq - invFreqDyn[length(invFreqDyn)]

par(mfrow=c(1,2))
plot(invFreqDyn, type='l', lwd=2, ylim=c(0,1), col=2)
lines(test$Fx.gen[,1], lwd=2, col=1)
lines(test$Fx.gen[,2], lwd=2, col=3) #AI
lines(test$Fx.gen[,3], lwd=2, col=4)# II
lines(test$Fx.gen[,4], lwd=2, col=5)
lines(test$Fx.gen[,5], lwd=2, col=6) #aI
lines(test$Fx.gen[,6], lwd=2, col=9)
lines(invFreqDyn, lwd=2, col=2)


plot(invFreqDyn, type='l', lwd=2, ylim=c(0,1), col=2)
lines(test$Fy.gen[,1], lwd=2, col=1)
lines(test$Fy.gen[,2], lwd=2, col=2) #AI
lines(test$Fy.gen[,3], lwd=2, col=3)
lines(test$Fy.gen[,4], lwd=2, lty=2, col=2) #aI
lines(test$Fy.gen[,5], lwd=2, col=4) 
lines(test$Fy.gen[,6], lwd=2, col=5)
