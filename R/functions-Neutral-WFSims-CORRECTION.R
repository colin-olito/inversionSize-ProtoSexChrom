################################################################
#' 	"Sheltering: CORRECTED"
#' 
#'  This file contains CORRECTED functions to perform deterministic and W-F simulations 
#'  for Y-LINKED inversions capturing the sex-determining locus, 
#'  and n selected loci with segregating partially recessive
#'  deleterious mutations.
#'
#'
#'  Author: Colin Olito
#'
#'  NOTES:  Functions based on supplementary code from 
#' 
#' 			Olito, Ponnikas, Hansson, Abbott (2024). Consequences
#' 			of partially recessive deleterious genetic variation 
#' 			for the evolution of inversions suppressing 
#' 			recombination between sex chromosomes. 
#' 			Evolution 78: 1499--1510. doi: 10.1093/evolut/qpae060.
#'
#'			NOTE that Olito et al. (2024) corrected and replaced 
#'			Olito et al. (2022), which contained an error in the
#'			functions for neutral Y-linked inversions.




######################################
#' Equilibrium Approximation WF model

#' Deterministic fitness expressions for
#' equilibrium approximation WF model 
w.YI.X  <-  function(n, r, t, u, s, h, qHat) {
	((qHat*(1 - s) + (1 - qHat)*(1 - s*h))^r) * (1 - u*(2 - exp(-s*h*t)))^(n - r)
}
w.Y.X  <-  function(n, u) {
	(1 - 2*u)^(n)
}


####################################################################
#' Function to estimate Pr(fix | x) for different inversion sizes (x)
#' Using equilibrium approximation WF simulations
makeDataPrFixInvSizeHaploid  <-  function(h = 0.25, s = 0.01, U.vals = c(0.02, 0.05, 0.1),
										  nTot = 10^4, N = 10^4, Nfname = "") {

	# Containers
	PrFix      <-  c()
	rFixedInv  <-  c()
	YI.t       <-  c()

	rFixedInvTab  <-  data.frame()
	rInvSize      <-  c()
	rNs           <-  c()
	rUs           <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Initial frequency: single-copy
	YI.0   <-  2/N

	#number of simulations 
	sims  = 200*N/2


	# Loop over Us factor values
	for(j in 1:length(U.vals)) {
		# mutation rate & initial qHat
		U     <-  U.vals[j]
		u     <-  U/nTot
		qHat  <-  u/(h*s)

		# Loop over inversion size
		for(k in 1:length(invSize)) {


			# Draw random values for # loci captured by inversion
			ns   <-  rpois(sims, lambda=nTot*invSize[k])

			# Draw random value for # del. mutations captured by inversion (r) given x
			rs  <-  c()
			for(m in 1:length(ns)) {
				rs[m]  <-  sum(rbinom(n=ns[m], size=1, prob=qHat))
			}

			# counter for fixations
			fix   = 0
			rFixedInv  <-  c()
			
			# Loop over replicate simulations
			for(l in 1:sims){

				# take randomly drawn n and r values
				n  <-  ns[l]
				r  <-  rs[l]

				# Assign frequencies
				# Note implicit assumption of equal initial
				# frequencies in XOv, XSp, Y chromosomes
				YI.t  <-  YI.0
				t     <-  1

				# Run forward simulation
				while(YI.t*(1 - YI.t) > 0) {
					# Fitness expressions
							w_YI.X  <-  w.YI.X(n=n, r=r, t=t, u=u, s=s, h=h, qHat=qHat)
							w_Y.X   <-  w.Y.X(n=n, u=u)
							w_avg   <-  YI.t*w_YI.X + (1 - YI.t)*w_Y.X
          		YI.sel  <-  YI.t*w_YI.X/w_avg

					#binomial sampling
					YI.t     <-  rbinom(1, (N/2), YI.sel) / (N/2)
					# time counter
					t  <-  t + 1
				}
				if(YI.t == 1) {
					fix        <-  fix + 1
					rFixedInv  <-  c(rFixedInv, rs[l])
				}
			}

			PrFix         <-  c(PrFix, (fix/sims))
			rTab          <-  as.data.frame(table(rFixedInv))
			rTab$Freq     <-  rTab$Freq/sims
			rFixedInvTab  <-  rbind(rFixedInvTab, rTab)
			rNs           <-  c(rNs, rep(N, times = nrow(rTab)))
			rUs           <-  c(rUs, rep(U.vals[j], times = nrow(rTab)))
			rInvSize      <-  c(rInvSize, rep(invSize[k], times = nrow(rTab)))
								cat('\r', paste("U: ", j, "/", length(U.vals), 
								", x: ", k, "/", length(invSize), " complete", sep=""))
		}
	}

	# Index variables
	Ns        <-  rep(N, times=(length(U.vals)*length(invSize)))
	Us        <-  rep(U.vals, each=(length(invSize)))
	invSizes  <-  rep(invSize, times=length(U.vals))
	
	# Export Results Dataframe
	filename  <-  paste("./output/data/simResults/PrFixFig_Haploid_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

	filename  <-  paste("./output/data/simResults/PrFixFig_Haploid_rFixedInv_h", h, "_s", s, Nfname, ".csv", sep="")
	r.d       <-  as.data.frame(cbind(rNs, rUs, rInvSize, rFixedInvTab))
	write.csv(r.d, file=filename, row.names=FALSE)

}





