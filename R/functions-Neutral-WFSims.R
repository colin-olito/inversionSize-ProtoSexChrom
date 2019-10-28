################################################################
#  NEUTRAL INVERSIONS
#  Wright-Fisher forward simulations to calculate Pr(fix | x) 
#  for a neutral Y-LINKED inversion capturing the sex-determining
#  locus
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

pr.fix.neutral = function(x, Ny, sd, Ud){
	1 / (Ny*exp(Ud*x/sd))
}

neutralFwdSimY  <-  function(par.list) {

	#unpack par.list
	Ny   <-  par.list$Ny
	sd   <-  par.list$sd
	Ud   <-  par.list$Ud

	#number of simulation runs per parameter set and x
	reps = max(200*Ny)

	#initial inversion frequency
	q.0 = 1/Ny

	# Wright-Fisher simulations
	# Storage vectors
	invSizes = 1:11/11
	simulated = rep(0, length(invSizes))

	for(j in 1:length(invSizes)){
		#inversion size
		x = invSizes[j]
  	
		#fixation events (initially 0)
		successes = 0

		if(j == 1)
			cat('\r', paste(0, "/ ", length(invSizes)))
		for(i in 1:reps){
			t = 0
			q = q.0

			while(q*(1 - q) > 0){
				#expected frequency after selection
				w.avg = q*exp(-Ud*x*(1 - exp(-sd*t))) + (1 - q)*exp(-Ud*x)
				F.sel = q + q*(1 - q)*(exp(-Ud*x*(1 - exp(-sd*t))) - exp(-Ud*x))/w.avg

				#binomial sampling
				q = rbinom(1, Ny, F.sel)/(Ny)
				t = t + 1
			}

			if(q > 0){successes = successes + 1} else{successes = successes + 0}
		}

		simulated[j] = successes/reps
		cat('\r', paste(j, "/ ", length(invSizes)))
	}

	PrMutFree = exp(-2*Ud*invSizes/sd)
	neutral.sim = PrMutFree*simulated

	res  <-  data.frame(
				  "Ny"           =  Ny,
				  "Ud"           =  Ud,
				  "sd"           =  sd,
				  "invSizes"     =  invSizes,
				  "PrFixSim"     =  simulated,
				  "reps"         =  reps,
				  "PrMutFree"    =  PrMutFree,
				  "PrFixSimDel"  =  neutral.sim
				  )
	# create file name
	filename  <-  paste("./output/data/simResults/neutral-Y-pFix", "_Ny", Ny, "_Ud", Ud, "_sd", sd, ".csv", sep="")

	# export dataframe
	write.csv(res, file=filename)
}