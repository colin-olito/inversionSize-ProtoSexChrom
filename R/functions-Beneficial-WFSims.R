################################################################
#  BENEFICIAL INVERSIONS
#  Wright-Fisher forward simulations to calculate Pr(fix | x) 
#  for an unconditionally beneficial inversion capturing the 
#  the sex-determining locus on the Y chromosome.
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

#####################
##  Functions

pr.fix.Ybeneficial = function(x, sI, sd, Ud){
	2*sI*(1 + ((Ud*x) / (1 - (1 - sI)*exp(-sd))) )*exp(-Ud*x/sd)
}

beneficialFwdSimY  <-  function(par.list) {

	#unpack par.list
	Ny   <-  par.list$Ny
	sI   <-  par.list$sI
	sd   <-  par.list$sd
	Ud   <-  par.list$Ud

	#number of simulation runs per parameter set and x
#	reps = max(10^5*Ny*sI)
	reps = 10^5

	#initial inversion frequency
	q.0 = 1/Ny

	# Storage vectors
	invSizes = 1:11/11
	simulated = rep(0, length(invSizes))

	# Calculate threshold frequency for inversion establishment
		pcrit        <-  4/(Ny*sI)

	for(j in 1:length(invSizes)){
		#inversion size
		x = invSizes[j]

		#fixation events (initially 0)
		successes = 0
  
		for(i in 1:reps){
			t = 0
			q = q.0
    
			while(q*(1 - q) > 0){
				#expected frequency after selection
#				w.avg = q*(1 + sI)*exp(-Ud*x*(1 - exp(-sd*t))) + (1 - q)*exp(-Ud*x)
#				F.sel = q + q*(1 - q)*((1 + sI)*exp(-Ud*x*(1 - exp(-sd*t))) - exp(-Ud*x))/w.avg
				w.avg = q*(1 + sI)*exp(-Ud*x*(2 - exp(-sd*t))) + (1 - q)*exp(-2*Ud*x)
				F.sel = q + q*(1 - q)*((1 + sI)*exp(-Ud*x*(2 - exp(-sd*t))) - exp(-2*Ud*x))/w.avg
      
				#binomial sampling
				q = rbinom(1, Ny, F.sel)/(Ny)
				t = t + 1
			}
    
			if(q > 0){successes = successes + 1} else{successes = successes + 0}
		}

		simulated[j] = successes/reps
		cat('\r', paste(j, "/ ", length(invSizes)))
	}
	PrMutFree = exp(-Ud*invSizes/sd)
	ben.sim = PrMutFree*simulated

	res  <-  data.frame(
				  "Ny"           =  Ny,
				  "sI"           =  sI,
				  "Ud"           =  Ud,
				  "sd"           =  sd,
				  "invSizes"     =  invSizes,
				  "PrFixSim"     =  simulated,
				  "reps"         =  reps,
				  "PrMutFree"    =  PrMutFree,
				  "PrFixSimDel"  =  ben.sim
				  )
	# create file name
	filename  <-  paste("./output/data/simResults/beneficial-Y-pFix", "_Ny", Ny, "_sI", sI, "_Ud", Ud, "_sd", sd, ".csv", sep="")

	# export dataframe
	write.csv(res, file=filename)
}