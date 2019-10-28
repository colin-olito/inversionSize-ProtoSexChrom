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
source('R/functions-SA-Y-2Locus-determSims.R')


##############################################################################################
##############################################################################################


#' Wright-Fisher simulations to estimate Pr(fix | x) for multi-locus SA model of a Y-linked 
#' inversion
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
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
simPrInvMultiLocusSA  <-  function(par.list) {

	# unpack par.list
	reps  =  par.list$reps
	N     =  par.list$N
	gen   =  par.list$gen
	sm    =  par.list$sm
	sf    =  par.list$sf
	hm    =  par.list$hm
	hf    =  par.list$hf
	r     =  par.list$r
	A     =  par.list$A
	P     =  par.list$P
	Ud    =  par.list$Ud
	sd    =  par.list$sd

	# storage vectors
	invSizes    <-  1:11/11
	simulated   <-  rep(0, length(invSizes))
	PrFixApprox <-  rep(0, length(invSizes))

	for(i in 1:length(invSizes)){

		# inversion size
		x  <-  invSizes[i]

		# fixation events (initially 0)
		successes  <-  0

			# equilibrium frequencies at SA loci
			par.list.A    <-  par.list
			par.list.A$r  <-  1/2
			SL.freqs  <-  getEqFreqsSAY(par.list, eq.threshold=1e-6)
			A.freqs   <-  getEqFreqsSAY(par.list.A, eq.threshold=1e-6)

			pFixApp  <-  c()

		for(j in 1:reps) {

			# time counter
			t  <-  0

			# number of locally SA loci captured by inversion
			n        <-  rpois(1, A*x)
			SAPos    <-  runif(n)
			n.slPAR  <-  sum(SAPos <= P)
			while(n.slPAR > 1) {
				SAPos    <-  runif(n)
				n.slPAR  <-  sum(SAPos <= P)
			}

			# overall selection coefficient on inversion
			qHat.aPAR   <-  A.freqs$aFreq.eq[3]
			qHat.slPAR  <-  SL.freqs$aFreq.eq[3]
			capturedAlleles.aPAR   <-  sum(rbinom(n = (n - n.slPAR), prob=qHat.aPAR, size=1))
			capturedAlleles.slPAR  <-  sum(rbinom(n = n.slPAR, prob=qHat.slPAR, size=1))

			s_nM  <-  ((n-n.slPAR)*sm*(1-qHat.aPAR)*(1 - hm*(1-2*qHat.aPAR)-qHat.aPAR) + 
						n.slPAR*sm*(1-qHat.slPAR)*(1 - hm*(1-2*qHat.slPAR)-qHat.slPAR)) -
					  ((n-n.slPAR-capturedAlleles.aPAR)*sm*(1 - hm*(1-2*qHat.aPAR)-qHat.aPAR) + 
						(n.slPAR-capturedAlleles.slPAR)*sm*(1 - hm*(1-2*qHat.slPAR)-qHat.slPAR))
			pFixApp[j]  <-  2*s_nM*(1 + ((Ud*x)/(1 - (1 - s_nM)*exp(-sd))))

			# Introduce inversion at low frequency
			q.0  <-  2/N
			q    <-  q.0
			
			# forward simulation loop
			while(q*(1 - q) > 0){
				# expected frequency after selection
				w.avg  <- q*(1 + s_nM)*exp(-Ud*x*(1 - exp(-sd*t))) + (1 - q)*exp(-Ud*x)
				F.sel  <- q + q*(1 - q)*((1 + s_nM)*exp(-Ud*x*(1 - exp(-sd*t))) - exp(-Ud*x))/w.avg

				#binomial sampling
				q  <-  rbinom(1, N, F.sel)/(N)
				t  <-  t + 1
			}
			if(q > 0){successes  <-  successes + 1} else{successes  <-  successes + 0}

		}
		print(c(x, successes/reps))
		simulated[i]    <-  successes/reps
		PrFixApprox[i]  <-  mean(pFixApp)
	}

	# calculate simulated probability of fixation
	Pr.mut.free  <-  exp(-2*Ud*invSizes/sd)
	PrFixSim     <-  simulated*Pr.mut.free
	PrFixApproxDel  <-  PrFixApprox*Pr.mut.free

	# Analytic Approximation
#	PrFixApprox  <-  2*s_nM*(1 + ((Ud*invSizes)/(1 - (1 - s_nM)*exp(-sd))))*qHat.aPAR^(n-1)*qHat.slPAR*(invSizes^(n+1))*A*exp(-x*A)*Pr.mut.free
#	PrFixApprox  <-  2*s_nM*(1 + ((Ud*invSizes)/(1 - (1 - ss_nM)*exp(-sd))))*Pr.mut.free

	# results
	res  <-  data.frame("invSize"      =  invSizes,
						"simPrFix"     =  simulated,
						"simPrFixDel"  =  PrFixSim,
						"approxPrFix"  =  PrFixApprox,
						"approxPrFixDel"  =  PrFixApproxDel,
						"Pr.mut.free"  =  Pr.mut.free
		)
	return(res)
}








#' compare multilocus s_I with 2-locus models
#' requires new arg in par.list: locSA, a vector of length two and sum = 2
#' providing the number of SA loci in the slPAR and aPar respectively 
#' (e.g., locSA can be c(0,2), c(1,1), c(2,0))
simPrInvMultiLocusSATwoLocus  <-  function(par.list) {

	# unpack par.list
	N      =  par.list$N
	gen    =  par.list$gen
	reps   =  par.list$reps
	locSA  = par.list$locSA
	sm     =  par.list$sm
	sf     =  par.list$sf
	hm     =  par.list$hm
	hf     =  par.list$hf
	r      =  par.list$r
	A      =  par.list$A
	P      =  par.list$P
	Ud     =  par.list$Ud
	sd     =  par.list$sd

	# storage vectors
	invSizes    <-  1:11/11
	simulated   <-  rep(0, length(invSizes))

	for(i in 1:length(invSizes)){

		# inversion size
		x  <-  invSizes[i]

		# fixation events (initially 0)
		successes  <-  0

			# equilibrium frequencies at SA loci
			par.list.A    <-  par.list
			par.list.A$r  <-  1/2
			SL.freqs  <-  getEqFreqsSAY(par.list, eq.threshold=1e-6)
			A.freqs   <-  getEqFreqsSAY(par.list.A, eq.threshold=1e-6)

			pFixApp  <-  c()

		for(j in 1:reps) {

			# time counter
			t  <-  0

			# overall selection coefficient on inversion
			n  <-  2
			qHat.aPAR   <-  A.freqs$aFreq.eq[3]
			qHat.slPAR  <-  SL.freqs$aFreq.eq[3]
			capturedAlleles.slPAR  <-  sum(rbinom(n = locSA[1], prob=qHat.slPAR, size=1))
			capturedAlleles.aPAR   <-  sum(rbinom(n = locSA[2], prob=qHat.aPAR, size=1))

			s_nM  <-  (locSA[2]*sm*(1-qHat.aPAR)*(1 - hm*(1-2*qHat.aPAR)-qHat.aPAR) + 
						locSA[1]*sm*(1-qHat.slPAR)*(1 - hm*(1-2*qHat.slPAR)-qHat.slPAR)) -
					  ((n-locSA[1]-capturedAlleles.aPAR)*sm*(1 - hm*(1-2*qHat.aPAR)-qHat.aPAR) + 
						(locSA[1]-capturedAlleles.slPAR)*sm*(1 - hm*(1-2*qHat.slPAR)-qHat.slPAR))

			# Introduce inversion at low frequency
			q.0  <-  2/N
			q    <-  q.0
			
			# forward simulation loop
			while(q*(1 - q) > 0){
				# expected frequency after selection
				w.avg  <- q*(1 + s_nM)*exp(-Ud*x*(1 - exp(-sd*t))) + (1 - q)*exp(-Ud*x)
				F.sel  <- q + q*(1 - q)*((1 + s_nM)*exp(-Ud*x*(1 - exp(-sd*t))) - exp(-Ud*x))/w.avg

				#binomial sampling
				q  <-  rbinom(1, N, F.sel)/(N)
				t  <-  t + 1
			}
			if(q > 0){successes  <-  successes + 1} else{successes  <-  successes + 0}

		}
		print(c(x, successes/reps))
		simulated[i]    <-  successes/reps
		PrFixApprox[i]  <-  mean(pFixApp)
	}

	# calculate simulated probability of fixation
	Pr.mut.free  <-  exp(-2*Ud*invSizes/sd)
	PrFixSim     <-  simulated*Pr.mut.free

	# Analytic Approximation
#	PrFixApprox  <-  2*s_nM*(1 + ((Ud*invSizes)/(1 - (1 - s_nM)*exp(-sd))))*qHat.aPAR^(n-1)*qHat.slPAR*(invSizes^(n+1))*A*exp(-x*A)*Pr.mut.free
#	PrFixApprox  <-  2*s_nM*(1 + ((Ud*invSizes)/(1 - (1 - ss_nM)*exp(-sd))))*Pr.mut.free

	# results
	res  <-  data.frame("invSize"      =  invSizes,
						"simPrFix"     =  simulated,
						"simPrFixDel"  =  PrFixSim,
						"Pr.mut.free"  =  Pr.mut.free
		)
	return(res)
}
