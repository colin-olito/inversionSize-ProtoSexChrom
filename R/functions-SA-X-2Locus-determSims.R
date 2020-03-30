##########################################################
#  X-linked inversion spanning the SDR and a second locus
#
#  Necessary functions for deterministic simulation
#  of adult genotypic frequency recursions for the
#  SEXUALLY ANTAGONISTIC SELECTION MODEL of the 
#  evolution of new inversions capturing a recessive   
#  female-determining allele at a sex-determining 
#  locus and a female-beneficial allele at a second
#  (SA) locus with partial linkage (within the PAR)
#
#
#  Author: Colin Olito
#
#  NOTES:  
#          

##########################
##  Key to adult genotypes
#
# Females
#  x1    =  XA  | XA   =  Fx[1]
#  x1I   =  XA  | XAI  =  Fx[2]
#  x1II  =  XAI | XAI  =  Fx[3]
#  x2    =  XA  | Xa   =  Fx[4]
#  x2I   =  XAI | Xa   =  Fx[5]
#  x3    =  Xa  | Xa   =  Fx[6]
#
# Males
#  y1    =  XA  | YA  =  Fy[1]
#  y1I   =  XAI | YA  =  Fy[2]
#  y2c   =  XA  | Ya  =  Fy[3]
#  y2cI  =  XAI | Ya  =  Fy[4]
#  y2t   =  Xa  | YA  =  Fy[5]
#  y3    =  Xa  | Ya  =  Fy[6]

# 'c' denotes 'cis-', and indicates that the 'A' allele is on the X in the heterozygote
# 't' denotes 'trans-', and indicates that the 'a' allele is on the X in the heterozygote
# 'I' indicates an inverted genotype, 'II' indicates an inversion homozgote


#######################################################
## Necessary Functions
#######################################################

# meiosis & random mating
# Females
x1m    <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2 + Fx[4]/2) * (Fy[1] + Fy[3]*(1 - r) + Fy[4] + Fy[5]*r)
}
x1Im   <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2 + Fx[4]/2) * (Fy[2] + Fy[4]) + 
	(Fx[2]/2 + Fx[3] + Fx[5]/2) * (Fy[1] + Fy[3]*(1 - r) + Fy[5]*r)
}
x1IIm  <-  function(Fx, Fy, r) {
	(Fx[2]/2 + Fx[3] + Fx[5]/2) * (Fy[2] + Fy[4])
}
x2m     <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2 + Fx[4]/2) * (Fy[3]*r + Fy[5]*(1 - r) + Fy[6]) +
	(Fx[2]/2 + Fx[4]/2 + Fx[6]) * (Fy[1] + Fy[3]*(1 - r) + Fy[5]*r)
}
x2Im   <-  function(Fx, Fy, r) {
	(Fx[2]/2 + Fx[3] + Fx[5]/2) * (Fy[3]*r + Fy[5]*(1 - r) + Fy[6]) +
	(Fx[4]/2 + Fx[5]/2 + Fx[6]) * (Fy[2] + Fy[4])
}
x3m    <-  function(Fx, Fy, r) {
	(Fx[4]/2 + Fx[5]/2 + Fx[6]) * (Fy[3]*r + Fy[5]*(1 - r) + Fy[6])
}

# Males
y1m    <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2 + Fx[4]/2) * (Fy[1] + Fy[2] + Fy[3]*r + Fy[5]*(1 - r))
}
y1Im   <-  function(Fx, Fy, r) {
	(Fx[2]/2 + Fx[3] + Fx[5]/2) * (Fy[1] + Fy[2] + Fy[3]*r + Fy[5]*(1 - r))
}
y2cm   <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2 + Fx[4]/2) * (Fy[3]*(1 - r) + Fy[4] + Fy[5]*r + Fy[6])
}
y2cIm  <-  function(Fx, Fy, r) {
	(Fx[2]/2 + Fx[3] + Fx[5]/2) * (Fy[3]*(1 - r) + Fy[4] + Fy[5]*r + Fy[6])
}
y2tm   <-  function(Fx, Fy, r) {
	(Fx[4]/2 + Fx[5]/2 + Fx[6]) * (Fy[1] + Fy[2] + Fy[3]*r + Fy[5]*(1 - r))
}
y3m    <-  function(Fx, Fy, r) {
	(Fx[4]/2 + Fx[5]/2 + Fx[6]) * (Fy[3]*(1 - r) + Fy[4] + Fy[5]*r + Fy[6])
}

# Selection
# Females
wfBar  <-  function(Fxm, Wf) {
	Fxm[1]*Wf[1] + Fxm[2]*Wf[1] + Fxm[3]*Wf[1] + Fxm[4]*Wf[2] + Fxm[5]*Wf[2] + Fxm[6]*Wf[3]
}
x1Pr  <-  function(Fxm, Wf) {
	Fxm[1]*Wf[1]/wfBar(Fxm=Fxm, Wf=Wf)
}
x1IPr  <-  function(Fxm, Wf) {
	Fxm[2]*Wf[1]/wfBar(Fxm=Fxm, Wf=Wf)
}
x1IIPr  <-  function(Fxm, Wf) {
	Fxm[3]*Wf[1]/wfBar(Fxm=Fxm, Wf=Wf)
}
x2Pr  <-  function(Fxm, Wf) {
	Fxm[4]*Wf[2]/wfBar(Fxm=Fxm, Wf=Wf)
}
x2IPr  <-  function(Fxm, Wf) {
	Fxm[5]*Wf[2]/wfBar(Fxm=Fxm, Wf=Wf)
}
x3Pr  <-  function(Fxm, Wf) {
	Fxm[6]*Wf[3]/wfBar(Fxm=Fxm, Wf=Wf)
}

# Males
wmBar  <-  function(Fym, Wm) {
	Fym[1]*Wm[1] + Fym[2]*Wm[1] + Fym[3]*Wm[2] + Fym[4]*Wm[2] + Fym[5]*Wm[2] + Fym[6]*Wm[3];
}
y1Pr  <-  function(Fym, Wm) {
	Fym[1]*Wm[1]/wmBar(Fym=Fym, Wm=Wm)
}
y1IPr  <-  function(Fym, Wm) {
	Fym[2]*Wm[1]/wmBar(Fym=Fym, Wm=Wm)
}
y2cPr  <-  function(Fym, Wm) {
	Fym[3]*Wm[2]/wmBar(Fym=Fym, Wm=Wm)
}
y2cIPr  <-  function(Fym, Wm) {
	Fym[4]*Wm[2]/wmBar(Fym=Fym, Wm=Wm)
}
y2tPr  <-  function(Fym, Wm) {
	Fym[5]*Wm[2]/wmBar(Fym=Fym, Wm=Wm)
}
y3Pr  <-  function(Fym, Wm) {
	Fym[6]*Wm[3]/wmBar(Fym=Fym, Wm=Wm)
}


###############################
##  Functions to translate to alternate coordinate system, 
##  where we track frequency of the 'a' allele at the selected 
##  locus among gametes in three classes of chromosomes: 
##  X's from females
##  X's from  males
##  Y's from males

Xf  <-  function(Fx) {
	Fx[1] + Fx[2] + Fx[3] + Fx[4]/2 + Fx[5]/2
}
Xm  <-  function(Fy, r) {
	Fy[1] + Fy[2] + Fy[3]*(1 - r) + Fy[4] + Fy[5]*r
}
Y   <-  function(Fy, r) {
	Fy[1] + Fy[2] + Fy[3]*r + Fy[5]*(1 - r)
}


###############################
##  Eigenvalue from Analytic results
##  indicating whether Inversion genotype
##  should be able to invade
lambdaY  <-  function(Wf, Xf, Xm) {
	(Wf[2] + Xf*(Wf[1] - Wf[2])) / (Wf[3]*(1 - Xf)*(1 - Xm) + Wf[1]*Xf*Xm + Wf[2]*(Xf + Xm*(1 - 2*Xf)))
}
lambdaYWkSel  <-  function(Xf, Xm, hf, sf) {
	1 + sf*(1 - Xm)*(1 - Xf - hf*(1 - 2*Xf)) 
}


#################################
##  Simulation functions
#################################



	
#' Find deterministic equilibrium frequencies in the absence of inversion
#'
#' @title Forward deterministic simulation of adult genotypic recursions for 2-locus model
#' 		  of PAR dynamics with an inversion genotype spanning the male SDR and the second locus
#' 
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   selType  =  "SA",
#' 				   selPars  =  c(hf, sf, hm, sm),
#'				   gen      =  5000,
#'				   r        =  1/2
#'				   )
#' @return	Returns a list with deterministic equilibrium frequencies 
#' 			for each genotype, and for the frequency of the A allele
#' 			on each of the three chromosome classes among gametes
#' @seealso 
#' @export
#' @author Colin Olito.
#' detSimYInversion(par.list, Fx.init, Fy.init, eqStart = TRUE, threshold = 1e-6) 
#'				   selType  =  "OD",
#' 				   selPars  =  c(sf1, sf2, sm1, sm2)
getEqFreqsX  <-  function(par.list, eq.threshold=1e-6, ...) {

	# unpack par.list
	selType  =  par.list$selType
	selPars  =  par.list$selPars
	gen      =  par.list$gen
	r        =  par.list$r

	# Sexually antagonistic fitness expressions
	if(selType == "SA") {
		Wf  <-  c(1,              1 - selPars[1]*selPars[2], 1 - selPars[2])
		Wm  <-  c(1 - selPars[4], 1 - selPars[3]*selPars[4], 1)
	}
	if(selType == "OD") {
		Wf  <-  c(1 - selPars[1], 1, 1 - selPars[2])
		Wm  <-  c(1 - selPars[3], 1, 1 - selPars[4])
	}

	# Arbitrary polymorphic initial frequencies
	Fx.eq     <-  c(1/3, 0,   0, 1/3,   0, 1/3)
	Fy.eq     <-  c(1/4, 0, 1/4,   0, 1/4, 1/4)
	AFreq.eq  <-  c(0,0,0)

	##  Generation Loop
	# Start simulation
	diffs  <-  rep(1,2)
	i=1
	while(i < (gen + 1) & any(abs(diffs[abs(diffs) > 0]) > eq.threshold)) {
		AFreq.eq  <-  c(Xf(Fx=Fx.eq), Xm(Fy=Fy.eq, r=par.list$r), Y(Fy=Fy.eq, r=r))
		Fxm  <-  round(c(x1m(  Fx=Fx.eq, Fy=Fy.eq, r=r),
						 x1Im( Fx=Fx.eq, Fy=Fy.eq, r=r),
						 x1IIm(Fx=Fx.eq, Fy=Fy.eq, r=r),
						 x2m(  Fx=Fx.eq, Fy=Fy.eq, r=r),
						 x2Im( Fx=Fx.eq, Fy=Fy.eq, r=r),
						 x3m(  Fx=Fx.eq, Fy=Fy.eq, r=r)), digits=8)
		Fym  <-  round(c(y1m(Fx=Fx.eq, Fy=Fy.eq, r=r),
						 y1Im(Fx=Fx.eq, Fy=Fy.eq, r=r),
						 y2cm(Fx=Fx.eq, Fy=Fy.eq, r=r),
						 y2cIm(Fx=Fx.eq, Fy=Fy.eq, r=r),
						 y2tm(Fx=Fx.eq, Fy=Fy.eq, r=r),
						 y3m(Fx=Fx.eq, Fy=Fy.eq, r=r)), digits=8)
		FxPr  <-  round(c(x1Pr(Fxm=Fxm, Wf=Wf),
						  x1IPr(Fxm=Fxm, Wf=Wf),
						  x1IIPr(Fxm=Fxm, Wf=Wf),
						  x2Pr(Fxm=Fxm, Wf=Wf),
						  x2IPr(Fxm=Fxm, Wf=Wf),
						  x3Pr(Fxm=Fxm, Wf=Wf)), digits=8)
		FyPr  <-  round(c(y1Pr(Fym=Fym, Wm=Wm),
						  y1IPr(Fym=Fym, Wm=Wm),
						  y2cPr(Fym=Fym, Wm=Wm),
						  y2cIPr(Fym=Fym, Wm=Wm),
						  y2tPr(Fym=Fym, Wm=Wm),
						  y3Pr(Fym=Fym, Wm=Wm)), digits=8)
		diffs  <-  abs(c(FxPr - Fx.eq, FyPr - Fy.eq))
		Fx.eq  <-  FxPr
		Fy.eq  <-  FyPr 
		i=i+1
	}

	res  <-  list(
				  "AFreq.eq"  =  AFreq.eq,
				  "Fx.eq"     =  Fx.eq,
				  "Fy.eq"     =  Fy.eq
				  )
	return(res)
}



#' Forward deterministic simulation of genotypic recursions for invasion of a new X-linked
#' inversion spanning the sex-determining locus an 
#'
#' @title Forward deterministic simulation of adult genotypic recursions for 2-locus model
#' 		  of PAR dynamics with an inversion genotype spanning the female SDR and a 
#' 		  female-benefit allele at the selected locus
#' 
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   selType  =  "SA",
#' 				   selPars  =  c(hf, sf, hm, sm),
#'				   gen      =  5000,
#'				   r        =  1/2
#'				   )
#' @param Fx.init A vector of length 3 containing initial genotypic frequencies for females
#' 				  Fx = c(x1, x2, x3)
#' @param Fy.init A vector of length 6 containing initial genotypic frequencies for males
#' 				  Fy = c(y1, y2c, y2cI, y2t, y3, y3I)
#' @return Returns a list with timeseries for each genotype, 
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#'				   selType  =  "OD",
#' 				   selPars  =  c(sf1, sf2, sm1, sm2)
detSimXInversion  <-  function(par.list, 
							   Fx.init = c(1/3, 0,   0, 1/3,   0, 1/3), 
							   Fy.init = c(1/4, 0, 1/4,   0, 1/4, 1/4), 
							   eqStart = TRUE, eq.threshold=1e-6) {




	if(any(c(round(sum(Fx.init), digits=4), round(sum(Fy.init), digits=4)) != 1))
		stop('Initial frequencies must sum to 1')


	# unpack par.list
	selType  =  par.list$selType
	selPars  =  par.list$selPars
	gen      =  par.list$gen
	r        =  par.list$r

	# Sexually antagonistic fitness expressions
	if(selType == "SA") {
		Wf  <-  c(1,              1 - selPars[1]*selPars[2], 1 - selPars[2])
		Wm  <-  c(1 - selPars[4], 1 - selPars[3]*selPars[4], 1)
	}
	if(selType == "OD") {
		Wf  <-  c(1 - selPars[1], 1, 1 - selPars[2])
		Wm  <-  c(1 - selPars[3], 1, 1 - selPars[4])
	}

	##  Initilize data storage structures
	Fx.gen  <-  matrix(0, ncol=6, nrow=par.list$gen)
	colnames(Fx.gen)  <-  c('x1','x1I','x1II','x2','x2I','x3')
	Fy.gen  <-  matrix(0, ncol=6, nrow=par.list$gen)
	colnames(Fy.gen)  <-  c('y1m','y1Im','y2cm','y2cIm','y2tm','y3')
	
	# Find equilibrium frequencies before introducing the inversion?
	if(eqStart) {
		initEq    <-  getEqFreqsX(par.list=par.list, eq.threshold=1e-6)
		AFreq.eq  <-  initEq$AFreq.eq
		Fx.eq     <-  initEq$Fx.eq
		Fy.eq     <-  initEq$Fy.eq
	} 

	# If starting at arbitrary initial equilibrium, 
	# use *.init frequencies
	if(!eqStart) {
		Fx.eq  <-  Fx.init
		Fy.eq  <-  Fy.init
	} 

	# Evaluate Eigenvalue associated with invasion 
	# of the inversion for initial equilibrium
	lambdaInv  <-  lambdaY(Wf=Wf, Xf=AFreq.eq[1], Xm=AFreq.eq[2])
	
	# Introduce inversion at low frequency
	mutProb  <-  c(Fx.eq[1], Fx.eq[4]/2, Fy.eq[1]/2, Fy.eq[3]/2)/sum(c(Fx.eq[1], Fx.eq[4]/2, Fy.eq[1]/2, Fy.eq[3]/2))
	invMut   <-  c(1,4,7,9)[as.vector(rmultinom(1,prob=mutProb, size=1)) == 1]
	if(invMut < 7) {
		Fx.eq[invMut]      <-  Fx.eq[invMut] - 1e-5
		Fx.eq[invMut + 1]  <-  Fx.eq[invMut + 1] + 1e-5
	}
	if(invMut > 6) {
		Fy.eq[(invMut - 6)]      <-  Fy.eq[(invMut - 6)] - 1e-5
		Fy.eq[(invMut - 6) + 1]  <-  Fy.eq[(invMut - 6) + 1] + 1e-5
	}

	# Initialize .gen storage
	Fx.gen[1,]  <-  Fx.eq
	Fy.gen[1,]  <-  Fy.eq
	
	# save initial equilibrium frequencies
	AFreq.init  <-  AFreq.eq

	# Generation loop
	diffs  <-  rep(1,12)
	i=2
	while(i < (gen + 1) & any(abs(diffs[abs(diffs) > 0]) > eq.threshold)) {
		AFreq.eq  <-  c(Xf(Fx=Fx.gen[i-1,]), Xm(Fy=Fy.gen[i-1,], r=r), Y(Fy=Fy.gen[i-1,], r=r))
		Fxm  <-  round(c(
						 x1m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 x1Im( Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 x1IIm(Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 x2m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 x2Im( Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 x3m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r)), digits=8)
		Fym  <-  round(c(
						 y1m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 y1Im( Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 y2cm( Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 y2cIm(Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 y2tm( Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
						 y3m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r)), digits=8)
		FxPr  <-  round(c(
						  x1Pr(  Fxm=Fxm, Wf=Wf),
						  x1IPr( Fxm=Fxm, Wf=Wf),
						  x1IIPr(Fxm=Fxm, Wf=Wf),
						  x2Pr(  Fxm=Fxm, Wf=Wf),
						  x2IPr( Fxm=Fxm, Wf=Wf),
						  x3Pr(  Fxm=Fxm, Wf=Wf)), digits=8)
		FyPr  <-  round(c(
						  y1Pr(  Fym=Fym, Wm=Wm),
						  y1IPr( Fym=Fym, Wm=Wm),
						  y2cPr( Fym=Fym, Wm=Wm),
						  y2cIPr(Fym=Fym, Wm=Wm),
						  y2tPr( Fym=Fym, Wm=Wm),
						  y3Pr(  Fym=Fym, Wm=Wm)), digits=8)
		FxPr        <-  FxPr/sum(FxPr)
		FyPr        <-  FyPr/sum(FyPr)
		diffs       <-  abs(c(FxPr - Fx.gen[i-1,], FyPr - Fx.gen[i-1,]))
		Fx.gen[i,]  <-  FxPr
		Fy.gen[i,]  <-  FyPr 
		i=i+1
	}	

	# Trim time series if necessary
	Fx.gen  <-  Fx.gen[rowSums(Fx.gen) > 0,]
	Fy.gen  <-  Fy.gen[rowSums(Fy.gen) > 0,]

	# Did inversion invade? 
	# Did the prediction based on the eigenvalue agree with the simulation results?
	lastGen      <-  nrow(Fx.gen)
	initInvFreq  <-  (2*((Fx.gen[1,      2]*1/2) + (Fx.gen[1      ,3]) + (Fx.gen[1      ,5]*1/2)) + (Fy.gen[1      ,2] + Fy.gen[1      ,4]))/3
	eqInvFreq    <-  (2*((Fx.gen[lastGen,2]*1/2) + (Fx.gen[lastGen,3]) + (Fx.gen[lastGen,5]*1/2)) + (Fy.gen[lastGen,2] + Fy.gen[lastGen,4]))/3
	simInvasion  <-  round(eqInvFreq, digits=8) > round(initInvFreq,digits=8)
	eigInvasion  <-  lambdaInv > 1

	#put names on AFreq.eq
	names(AFreq.eq)  <-  c("Xf","Xm","Y")

	# Put results in a list
	res  <-  list ("maxRuntime"   =  nrow(Fx.gen) >= gen,
				   "Fx.gen"       =  Fx.gen,
				   "Fy.gen"       =  Fy.gen,
				   "invAgree"     =  simInvasion == eigInvasion,
				   "eigInvasion"  =  eigInvasion,
				   "simInvasion"  =  simInvasion,
				   "lambdaInv"    =  lambdaInv,
				   "AFreq.init"   =  AFreq.init,
				   "AFreq.eq"     =  AFreq.eq,
				   "eqInvFreq"    =  eqInvFreq,
				   "gen"          =  nrow(Fx.gen)
				   )

	# Return results
	return(res)
}





#' Wrapper function to find equilibrium frequencies 
#' of X inversions across selection parameter space
detXSAInversionEQFreqs  <-  function(selRange = c(0,1), by=0.1, generations = 10000,
								   dom = c(1/2, 1/2), r = 1/2, eq.threshold = 1e-6) {

	# define parameter space being explored
	sfs      <-  seq(from=selRange[1], to=selRange[2], by=by)
	len.s    <-  length(sfs)
	sms      <-  sfs

	# data storage
	Xf.init  <-  rep(0, length=(len.s^2))
	Xm.init  <-  rep(0, length=(len.s^2))
	Y.init   <-  rep(0, length=(len.s^2))
	Xf.eq    <-  rep(0, length=(len.s^2))
	Xm.eq    <-  rep(0, length=(len.s^2))
	Y.eq     <-  rep(0, length=(len.s^2))
	inv.eq   <-  rep(0, length=(len.s^2))
	eigInv   <-  rep(0, length=(len.s^2))
	simInv   <-  rep(0, length=(len.s^2))

	# Simulation Loop
	cat('\r', paste("Progress", 0, "/ ", len.s))

	for(i in 1:len.s) {
		for(j in 1:length(sms)) {
			selPars  <-  c(dom[1], sfs[i], dom[2], sms[j])
			par.list  <-  list(
								selType  =  "SA",
								selPars  =  selPars,
								gen      =  generations,
								r        =  r
							   )

			# Run Deterministic Fwd Simulation
			Res  <-  detSimXInversion(par.list=par.list, eq.threshold=eq.threshold)

			# Store results
			Xf.init[(i-1)*len.s + j]  <-  Res$AFreq.init[1]
			Xm.init[(i-1)*len.s + j]  <-  Res$AFreq.init[2]
			Y.init[(i-1)*len.s + j]   <-  Res$AFreq.init[3]
			Xf.eq[(i-1)*len.s + j]    <-  Res$AFreq.eq[1]
			Xm.eq[(i-1)*len.s + j]    <-  Res$AFreq.eq[2]
			Y.eq[(i-1)*len.s + j]     <-  Res$AFreq.eq[3]
			inv.eq[(i-1)*len.s + j]   <-  Res$eqInvFreq
			eigInv[(i-1)*len.s + j]   <-  Res$eigInvasion
			simInv[(i-1)*len.s + j]   <-  Res$simInvasion

		}
		cat('\r', paste("Progress", i, "/ ", len.s))

	}

	# Vectors for data.frame
	recomb  <-  rep(r, times=len.s^2)
	hfs     <-  rep(selPars[1], times=len.s^2)
	hms     <-  rep(selPars[3], times=len.s^2)
	Sfs     <-  c()
	for(i in 1:len.s) {
		Sfs  <-  c(Sfs, rep(sfs[i], times = len.s))
	}
	Sms     <-  rep(sms, times=len.s)

	# Compile results in a data.frame
	results.df  <-  data.frame(
						   "hf"       =  hfs,
						   "sf"       =  Sfs,
						   "hm"       =  hms,
						   "sm"       =  Sms,
						   "r"        =  recomb,
						   "Xf.init"  =  Xf.init,
						   "Xm.init"  =  Xm.init,
						   "Y.init"   =  Y.init,
						   "Xf.eq"    =  Xf.eq,
						   "Xm.eq"    =  Xm.eq,
						   "Y.eq"     =  Y.eq,
						   "inv.eq"   =  inv.eq,
						   "eigInv"   =  eigInv,
						   "simInv"   =  simInv
							)

		# export data as .csv to ./output/data
		filename <-  paste("./output/data/simResults/deterministicSims-X", "_sMax", selRange[2], "_r", r, "_hf", selPars[1], "_hm", selPars[3], ".csv", sep="")
		write.csv(results.df, file=filename, row.names = FALSE)

}