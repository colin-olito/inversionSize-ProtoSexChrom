##########################################################
#  Y-linked inversion spanning the SDR and a second locus
#
#  Necessary functions for stochastic WF simulation
#  of adult genotypic frequencies for the SEXUALLY
#  ANTAGONISTIC SELECTION MODEL of the evolution of 
#  new inversions spanning the male sex-determining 
#  locus and a second (SA) locus within the PAR region.
#
#  Author: Colin Olito
#
#  NOTES:  These simulations examine how well the 
#			invasion fitness of the new inversion, 
#			which is approximated from the leading
#			eigenvalue of the local stability analysis,
#			approximates the probability of invasion
#			under weak selection using  
#          

##########################
##  Key to adult genotypes
#
# Females
#  x1  =  XA | XA  =  Fx[1]
#  x2  =  XA | Xa  =  Fx[2]
#  x3  =  Xa | Xa  =  Fx[3]
#
# Males
#  y1    =  XA | YA   =  Fy[1]
#  y2c   =  XA | Ya   =  Fy[2]
#  y2cI  =  XA | YaI  =  Fy[3]
#  y2t   =  Xa | YA   =  Fy[4]
#  y3    =  Xa | Ya   =  Fy[5]
#  y3I   =  Xa | YaI  =  Fy[6]
#
# 'c' denotes 'cis-', and indicates that the 'A' allele is on the X in the heterozygote
# 't' denotes 'trans-', and indicates that the 'a' allele is on the X in the heterozygote
# 'I' indicates an inverted genotype



#######################################################
## Necessary Functions
#######################################################
source('./R/functions-SA-Y-determSims.R')


#################################
##  Simulation functions
#################################


#' Forward deterministic simulation of genotypic recursions for
#' inversion
#'
#' @title Forward deterministic simulation of adult genotypic recursions for 2-locus model
#' 		  of PAR dynamics with an inversion genotype spanning the male SDR and the second locus
#' 
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   gen    =  5000,
#'				   N      =  1000,				
#'				   sm     =  0.1,
#'				   sf     =  0.1,
#'				   hm     =  1/2,
#'				   hf     =  1/2,
#'				   r      =  1/2
#'				   )
#' @param Fx.init A vector of length 3 containing initial genotypic frequencies for females
#' 				  Fx = c(x1, x2, x3)
#' @param Fy.init A vector of length 6 containing initial genotypic frequencies for males
#' 				  Fy = c(y1, y2c, y2cI, y2t, y3, y3I)
#' @return Returns a list with final frequencies of all genotypes and invasion info. 
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' invFwdSimYInversion(par.list, Fx.init, Fy.init) 

invFwdSimYInversion  <-  function(par.list, Fx.init, Fy.init, ...) {
	##  Warnings
#	if(any(par.list[3:6] < 0) | any(par.list[3:6] > 1) | par.list$r > 0.5)
#		stop('The chosen parameter values fall outside of biologically plausible bounds')

	if(any(c(round(sum(Fx.init), digits=4), round(sum(Fy.init), digits=4)) != 1))
		stop('Initial frequencies must sum to 1')

	# unpack par.list
	N    =  par.list$N
	sm   =  par.list$sm
	sf   =  par.list$sf
	hm   =  par.list$hm
	hf   =  par.list$hf
	r    =  par.list$r

	# Sexually antagonistic fitness expressions
	 Wf  <-  c(1,      1 - hf*sf, 1 - sf)
	 Wm  <-  c(1 - sm, 1 - hm*sm, 1)

	##  Initilize data storage structures
	Fx.gen  <-  matrix(0, ncol=3, nrow=par.list$gen)
	colnames(Fx.gen)  <-  c('x1', 'x2', 'x3')
	Fy.gen  <-  matrix(0, ncol=6, nrow=par.list$gen)
	colnames(Fy.gen)  <-  c('y1', 'y2c', 'y2cI', 'y2t', 'y3', 'y3I')

	# Find equilibrium frequencies in absence of inversion
	Fx.eq  <-  Fx.init
	Fy.eq  <-  Fy.init
	aFreq.eq  <-  c(Xf(Fx=Fx.eq), Xm(Fy=Fy.eq, r=par.list$r), Y(Fy=Fy.eq, r=r))

	# Evaluate Eigenvalue associated with invasion 
	# of the inversion for initial equilibrium
	lambdaInv  <-  lambdaYWkSel(Fx=Fx.eq, Fy=Fy.eq, r=r, hm=hm, sm=sm)

	# Define threshold frequency for establishment of inversion
	pcrit  <-  4/(N*(lambdaInv - 1))
	if(r == 0) {
		pcrit  <-  0.15
	}

	# warning for high pcrit
	if(pcrit > 0.5)
		stop('pcrit is quite high... reconsider pop. size and other parameters')

	# Introduce inversion at low frequency
	invMut  <-  c(2,5)[as.vector(rmultinom(1,prob=Fy.eq[c(2,5)], size=1)) == 1]
	Fy.eq[invMut]  <-  Fy.eq[invMut] - (2/N)
	if(invMut == 2) {
		Fy.eq[3]  <-  Fy.eq[3] + (2/N)		
	}
	if(invMut == 5) {
		Fy.eq[6]  <-  Fy.eq[6] + (2/N)		
	}

	# Initialize .gen storage
	Fx.gen  <-  Fx.eq
	Fy.gen  <-  Fy.eq

	# Initial inversion frequency 
	InvFreq    <-  sum(Fy.eq[c(3,6)])
	E.InvFreq  <-  InvFreq
		
	## Start forward simulation with newly introduced inversion
	gen  <-  1

	# use threshold frequency for establishment of inversion
	while(InvFreq > 0 & InvFreq < pcrit) {		
		
		## Step through recursions:
		# 1) Meiosis & random mating  --> Offspring genotype frequencies
		aFreq.eq  <-  c(Xf(Fx=Fx.gen), Xm(Fy=Fy.gen, r=r), Y(Fy=Fy.gen, r=r))
		Fxm  <-  round(c(x1m(  Fx=Fx.gen, Fy=Fy.gen, r=r),
					  x2m(  Fx=Fx.gen, Fy=Fy.gen, r=r),
					  x3m(  Fx=Fx.gen, Fy=Fy.gen, r=r)), digits=8)
		Fym  <-  round(c(y1m(  Fx=Fx.gen, Fy=Fy.gen, r=r),
					  y2cm( Fx=Fx.gen, Fy=Fy.gen, r=r),
					  y2cIm(Fx=Fx.gen, Fy=Fy.gen, r=r),
					  y2tm( Fx=Fx.gen, Fy=Fy.gen, r=r),
					  y3m(  Fx=Fx.gen, Fy=Fy.gen, r=r),
					  y3Im( Fx=Fx.gen, Fy=Fy.gen, r=r)), digits=8)
		# 4) Expected frequencies
		FxPr  <-  round(c(x1Pr( Fxm=Fxm, Wf=Wf),
					  x2Pr(  Fxm=Fxm, Wf=Wf),
					  x3Pr(  Fxm=Fxm, Wf=Wf)), digits=8)
		FyPr  <-  round(c(y1Pr( Fym=Fym, Wm=Wm),
					  y2cPr( Fym=Fym, Wm=Wm),
					  y2cIPr(Fym=Fym, Wm=Wm),
					  y2tPr( Fym=Fym, Wm=Wm),
					  y3Pr(  Fym=Fym, Wm=Wm),
					  y3IPr( Fym=Fym, Wm=Wm)), digits=8)
		# 5) Draw random frequencies in adults
		Fx.gen      <-  as.vector(rmultinom(1, N/2, FxPr)/(N/2))
		Fy.gen      <-  as.vector(rmultinom(1, N/2, FyPr)/(N/2))
		# Realized frequencies
		InvFreq    <-  sum(Fy.gen[c(3,6)])
		E.InvFreq  <-  sum(FyPr[c(3,6)])
		# next gen
		gen        <-  gen + 1
	}

	# Has the inversion reached threshold frequency for establishment (pcrit)? 
	# When did it first reach pcrit?
	if(InvFreq >= pcrit) {
		invEst      <-  1
		invEstTime  <-  gen
	} else {	
		invEst      <-  0
		invEstTime  <-  NA
	}
	res  <-  list(
				"FxPr"        =  FxPr,
				"FyPr"        =  FyPr,
				"Fx.gen"      =  Fx.gen,
				"Fy.gen"      =  Fy.gen,
				"aFreq"       =  aFreq.eq,
				"InvFreq"     =  InvFreq, 
				"E.InvFreq"   =  E.InvFreq,
				"InvEst"      =  invEst,
				"InvEstTime"  =  invEstTime,
				"nGen"        =  gen
			 	)
	return(res)
}




#' Wrapper function to run replicate forward simulations to validate
#' using (lambda - 1) to approximate  the invasion fitness of a Y-linked
#' inversion in a Wright-Fisher population.
#'
#' @title Wright-Fisher forward simulation of genotypic frequencies (default parameter values in parentheses)
#' @param nReps  Number of replicate simulations. With no deleterious mutations, and introducing 
#' 				 a single copy of the inversion.
#' @param N      Effective population size
#' @param hf     Dominance coefficient for a allele at selected locus in females
#' @param sf     Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02) in females
#' @param hm     Dominance coefficient for a allele at selected locus in males
#' @param sm     Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02) in males
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1)
#' @param n      Number of loci at which deleterious mutations may occur
#' @param u      Mutation rate (default value of u = 1e-6)
#' @param h.del  Dominance of deleterious mutations (default value of h = 0).
#' @param noDel  Omit fitness effects of deleterious mutations? Default value of FALSE assumes that 
#' 				       selection against deleterious mutations is twice as strong as selection favouring
#' 				       the locally adaptive alleles. Setting noDel = TRUE runs the W-F simulations as if
#' 				       there were no delterious mutations segregating in the population that are linked
#' 				       to the loci involved in local adaptation. 
#' @param fastSim    Logical. Use threshold frequency for establishment of inversion? 
#' @param newMutant  Switch to choose whether to specify new mutant genotypes, or if they are 
#'                   chosen randomly, given initial genotypic frequencies (Fii.f.init and Fii.m.init). 
#'                   A 2 positions vector ("m"/"f"/"random", "numeric"/"random"). the first position specify wether the inversions come in males, 
#'                   in females, or randomnly. The second select either a numerical position in the haplotype vector for the inversion genotype to be created
#'                   or wether it is random.
#' @seealso `invFwdSimYInversion`
#' @export
#' @author Colin Olito
runWFSimsYinversionValidatePinvApprox  <-  function(par.list) { 

	reps  <-  4*par.list$N

	if(length(par.list$sf) > 1) {
		sfs          <-  par.list$sf
		pInv         <-  c()
		pInvApprox   <-  c()
		pInvApprox2  <-  c()
		for(i in 1:length(sfs)) {
			par.list$sf  <-  sfs[i]
			initEq       <-  getEqFreqs(par.list=par.list, eq.threshold=1e-6)
			aFreq.eq     <-  initEq$aFreq.eq
			Fx.eq        <-  initEq$Fx.eq
			Fy.eq        <-  initEq$Fy.eq
			lYInvWkSel   <-  lambdaYWkSel(Fx=Fx.eq, Fy=Fy.eq, r=par.list$r, hm=par.list$hm, sm=par.list$sm)
			lYInvWkSel2  <-  lambdaYWkSel2(Xf=aFreq.eq[1], Y=aFreq.eq[3], hm=par.list$hm, sm=par.list$sm)

			inv  <-  c()
			for(j in 1:reps) {
				test    <-  invFwdSimYInversion(par.list=par.list, Fx.init=Fx.eq, Fy.init=Fy.eq)
				inv[j]  <-  test$InvEst	
			}
			print(i/length(sfs))

			pInv[i]         <-  sum(inv)/reps
			pInvApprox[i]   <-  2*(lYInvWkSel - 1)
			pInvApprox2[i]  <-  2*(lYInvWkSel2 - 1)
		}

		# compile results as data frame
		results.df  <-  as.data.frame(cbind(
											rep(par.list$N,  length(sfs)),
											rep(par.list$hf, length(sfs)),
											sfs,
											rep(par.list$hm, length(sfs)),
											rep(par.list$sm, length(sfs)),
											rep(par.list$r,  length(sfs)),
											pInv,
											pInvApprox,
											pInvApprox2
									  )     )
		colnames(results.df)  <-  c("N",
									"hf",
									"sf",
									"hm",
									"sm",
									"r",
									"pInv",
									"pInvApprox",
									"pInvApprox2"
									)
		# export data as .csv to ./output/data
		filename <-  paste("./output/data/simResults/SA/SA_Y_InvvalidatePinv_sfGrad", "_N", par.list$N, "_hf", par.list$hf, "_hm", par.list$hm, "_sm", par.list$sm, "_r", par.list$r, ".csv", sep="")
		write.csv(results.df, file=filename, row.names = FALSE)
	}


	if(length(par.list$r) > 1) {
		rs           <-  par.list$r
		pInv         <-  c()
		pInvApprox   <-  c()
		pInvApprox2  <-  c()
		for(i in 1:length(rs)) {
			par.list$r   <-  rs[i]
			initEq       <-  getEqFreqs(par.list=par.list, eq.threshold=1e-6)
			aFreq.eq     <-  initEq$aFreq.eq
			Fx.eq        <-  initEq$Fx.eq
			Fy.eq        <-  initEq$Fy.eq
			lYInvWkSel   <-  lambdaYWkSel(Fx=Fx.eq, Fy=Fy.eq, r=par.list$r, hm=par.list$hm, sm=par.list$sm)
			lYInvWkSel2  <-  lambdaYWkSel2(Xf=aFreq.eq[1], Y=aFreq.eq[3], hm=par.list$hm, sm=par.list$sm)

			inv  <-  c()
			for(j in 1:reps) {
				test    <-  invFwdSimYInversion(par.list=par.list, Fx.init=Fx.eq, Fy.init=Fy.eq)
				inv[j]  <-  test$InvEst	
			}
			print(i/length(rs))

			pInv[i]         <-  sum(inv)/reps
			pInvApprox[i]   <-  2*(lYInvWkSel - 1)
			pInvApprox2[i]  <-  2*(lYInvWkSel2 - 1)
		}

		# compile results as data frame
		results.df  <-  as.data.frame(cbind(
											rep(par.list$N,  length(rs)),
											rep(par.list$hf, length(rs)),
											rep(par.list$sf, length(rs)),
											rep(par.list$hm, length(rs)),
											rep(par.list$sm, length(rs)),
											rs,
											pInv,
											pInvApprox,
											pInvApprox2
									  )     )
		colnames(results.df)  <-  c("N",
									"hf",
									"sf",
									"hm",
									"sm",
									"r",
									"pInv",
									"pInvApprox",
									"pInvApprox2"
									)
		# export data as .csv to ./output/data
		filename <-  paste("./output/data/simResults/SA/SA_Y_InvvalidatePinv_rGrad", "_N", par.list$N, "_hf", par.list$hf, "_sf", par.list$sf, "_hm", par.list$hm, "_sm", par.list$sm, ".csv", sep="")
		write.csv(results.df, file=filename, row.names = FALSE)
	}
}
