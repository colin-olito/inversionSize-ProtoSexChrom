##########################################################
#  Y-linked inversion spanning the SDR and a second locus
#
#  Necessary functions for deterministic simulation
#  of adult genotypic frequency recursions for the
#  SEXUALLY ANTAGONISTIC SELECTION MODEL for the 
#  evolution of new inversions capturing the male
#  dominant sex-determining allele at a sex- 
#  determining locus and a male-beneficial 
#  allele at a second (SA) locus with partial 
#  linkage (within the PAR)
#
#  Author: Colin Olito
#
#  NOTES:  
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

# meiosis & random mating
# Females
x1m  <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2) * (Fy[1] + Fy[2]*(1 - r) + Fy[3] + Fy[4]*r)
}
x2m  <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2) * (Fy[2]*r + Fy[4]*(1 - r) + Fy[5] + Fy[6]) + 
	(Fx[3] + Fx[2]/2) * (Fy[1] + Fy[2]*(1 - r) + Fy[3] + Fy[4]*r)
}
x3m  <-  function(Fx, Fy, r) {
	(Fx[3] + Fx[2]/2) * (Fy[2]*r + Fy[4]*(1 - r) + Fy[5] + Fy[6])
}

# Males
y1m  <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2) * (Fy[1] + Fy[2]*r + Fy[4]*(1 - r))
}
y2cm  <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2) * (Fy[2]*(1 - r) + Fy[4]*r + Fy[5])
}
y2cIm  <-  function(Fx, Fy, r) {
	(Fx[1] + Fx[2]/2) * (Fy[3] + Fy[6])
}
y2tm  <-  function(Fx, Fy, r) {
	(Fx[3] + Fx[2]/2) * (Fy[1] + Fy[2]*r + Fy[4]*(1 - r))
}
y3m  <-  function(Fx, Fy, r) {
	(Fx[3] + Fx[2]/2) * (Fy[2]*(1 - r) + Fy[4]*r + Fy[5])
}
y3Im  <-  function(Fx, Fy, r) {
	(Fx[3] + Fx[2]/2) * (Fy[3] + Fy[6])
}

# Selection
# Females
wfBar  <-  function(Fxm, Wf) {
	Fxm[1]*Wf[1] + Fxm[2]*Wf[2] + Fxm[3]*Wf[3]
}
x1Pr  <-  function(Fxm, Wf) {
	Fxm[1]*Wf[1]/wfBar(Fxm=Fxm, Wf=Wf)
}
x2Pr  <-  function(Fxm, Wf) {
	Fxm[2]*Wf[2]/wfBar(Fxm=Fxm, Wf=Wf)
}
x3Pr  <-  function(Fxm, Wf) {
	Fxm[3]*Wf[3]/wfBar(Fxm=Fxm, Wf=Wf)
}
# Males
wmBar  <-  function(Fym, Wm) {
	Fym[1]*Wm[1] + Fym[2]*Wm[2] + Fym[3]*Wm[2] + Fym[4]*Wm[2] + Fym[5]*Wm[3] + Fym[6]*Wm[3];
}
y1Pr  <-  function(Fym, Wm) {
	Fym[1]*Wm[1]/wmBar(Fym=Fym, Wm=Wm)
}
y2cPr  <-  function(Fym, Wm) {
	Fym[2]*Wm[2]/wmBar(Fym=Fym, Wm=Wm)
}
y2cIPr  <-  function(Fym, Wm) {
	Fym[3]*Wm[2]/wmBar(Fym=Fym, Wm=Wm)
}
y2tPr  <-  function(Fym, Wm) {
	Fym[4]*Wm[2]/wmBar(Fym=Fym, Wm=Wm)
}
y3Pr  <-  function(Fym, Wm) {
	Fym[5]*Wm[3]/wmBar(Fym=Fym, Wm=Wm)
}
y3IPr  <-  function(Fym, Wm) {
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
	Fx[3] + Fx[2]/2	
}
Xm  <-  function(Fy, r) {
	Fy[2]*r + Fy[4]*(1 - r) + Fy[5] + Fy[6]
}
Y   <-  function(Fy, r) {
	Fy[2]*(1 - r) + Fy[3] + Fy[4]*r + Fy[5] + Fy[6]
}

###############################
##  Eigenvalue from Analytic results
##  indicating whether Inversion genotype
##  should be able to invade
lambdaY  <-  function(Wm, Xf, Y) {
	(Wm[2]*(1 - Xf) + Wm[3]*Xf)/(Wm[1]*(1 - Xf)*(1 - Y) + Wm[3]*Xf*Y + Wm[2]*(Xf + Y - 2*Xf*Y))
}

lambdaYFgen  <-  function(Fx, Fy, r, Wm) {
	(Wm[3]*(2 - 2*Fx[1] - Fx[2]) + Wm[2]*(2*Fx[1] + Fx[2]))/(Wm[3]*(2 - 2*Fx[1] - Fx[2])*(1 - Fy[1] - r*(Fy[2] - Fy[4]) - Fy[4]) + 
 Wm[1]*(2*Fx[1] + Fx[2])*(Fy[1] + r*(Fy[2] - Fy[4]) + Fy[4]) + Wm[2]*(2*(Fy[1] + r*Fy[2] + Fy[4] - r*Fy[4]) + 
    Fx[2]*(1 - 2*Fy[1] - 2*r*Fy[2] - 2*Fy[4] + 2*r*Fy[4]) + Fx[1]*(2 - 4*Fy[1] - 4*r*Fy[2] - 4*Fy[4] + 4*r*Fy[4])))
}

lambdaYWkSel  <-  function(Fx, Fy, r, hm, sm) {
	1 + (sm*(2*Fx[1] + Fx[2] + 2*hm*(1 - 2*Fx[1] - Fx[2]))*(Fy[1] + r*(Fy[2] - Fy[4]) + Fy[4]))/2
}

lambdaYWkSel2  <-  function(Xf, Y, hm, sm) {
	1 + sm*(1 - Y)*(1 - Xf - hm*(1 - 2*Xf))
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
#'				   gen    =  5000,
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
#' @return Returns a list with timeseries for each genotype, 
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' detSimYInversion(par.list, Fx.init, Fy.init, eqStart = TRUE, threshold = 1e-6) 

getEqFreqsSAY  <-  function(par.list, eq.threshold=1e-6, ...) {

	# unpack par.list
	gen  =  par.list$gen
	sm   =  par.list$sm
	sf   =  par.list$sf
	hm   =  par.list$hm
	hf   =  par.list$hf
	r    =  par.list$r

	# Sexually antagonistic fitness expressions
	 Wf  <-  c(1,      1 - hf*sf, 1 - sf)
	 Wm  <-  c(1 - sm, 1 - hm*sm, 1)

	# Arbitrary polymorphic initial frequencies
	Fx.eq     <-  c(0.25, 0.5, 0.25)
	Fy.eq     <-  c(0.25, 0.25, 0, 0.25, 0.25, 0)
	aFreq.eq  <-  c(0,0,0)

	##  Generation Loop
	# Start simulation
	diffs  <-  rep(1,9)
	i=1
	while(i < (gen + 1) & any(abs(diffs[abs(diffs) > 0]) > eq.threshold)) {
		aFreq.eq  <-  c(Xf(Fx=Fx.eq), Xm(Fy=Fy.eq, r=par.list$r), Y(Fy=Fy.eq, r=r))
		Fxm  <-  round(c(x1m(  Fx=Fx.eq, Fy=Fy.eq, r=r),
				   x2m(  Fx=Fx.eq, Fy=Fy.eq, r=r),
				   x3m(  Fx=Fx.eq, Fy=Fy.eq, r=r)), digits=8)
		Fym  <-  round(c(y1m(  Fx=Fx.eq, Fy=Fy.eq, r=r),
				   y2cm( Fx=Fx.eq, Fy=Fy.eq, r=r),
				   y2cIm(Fx=Fx.eq, Fy=Fy.eq, r=r),
				   y2tm( Fx=Fx.eq, Fy=Fy.eq, r=r),
				   y3m(  Fx=Fx.eq, Fy=Fy.eq, r=r),
				   y3Im( Fx=Fx.eq, Fy=Fy.eq, r=r)), digits=8)
		FxPr  <-  round(c(x1Pr( Fxm=Fxm, Wf=Wf),
				   x2Pr(  Fxm=Fxm, Wf=Wf),
				   x3Pr(  Fxm=Fxm, Wf=Wf)), digits=8)
		FyPr  <-  round(c(y1Pr( Fym=Fym, Wm=Wm),
				   y2cPr( Fym=Fym, Wm=Wm),
				   y2cIPr(Fym=Fym, Wm=Wm),
				   y2tPr( Fym=Fym, Wm=Wm),
				   y3Pr(  Fym=Fym, Wm=Wm),
				   y3IPr( Fym=Fym, Wm=Wm)), digits=8)
		diffs  <-  abs(c(FxPr - Fx.eq, FyPr - Fy.eq))
		Fx.eq  <-  FxPr
		Fy.eq  <-  FyPr 
		i=i+1
	}

	res  <-  list(
				  "aFreq.eq"  =  aFreq.eq,
				  "Fx.eq"     =  Fx.eq,
				  "Fy.eq"     =  Fy.eq
				  )
	return(res)
}



#' Forward deterministic simulation of genotypic recursions for invasion of a Y-linked inversion
#'
#' @title Forward deterministic simulation of adult genotypic recursions for 2-locus model
#' 		  of PAR dynamics with an inversion genotype spanning the male SDR and the second locus
#' 
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   gen    =  5000,
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
#' @return Returns a list with timeseries for each genotype, 
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' detSimYInversion(par.list, Fx.init, Fy.init, eqStart = TRUE, threshold = 1e-6) 

detSimYInversion  <-  function(par.list, 
							   Fx.init = c(0.25, 0.5, 0.25), 
							   Fy.init = c(0.25, 0.25, 0, 0.25, 0.25, 0), 
							   eqStart = TRUE, eq.threshold=1e-6) {
	##  Warnings
#	if(any(par.list[2:6] < 0) | any(par.list[2:6] > 1) | par.list$r > 0.5)
#		stop('The chosen parameter values fall outside of biologically plausible bounds')

	if(any(c(round(sum(Fx.init), digits=4), round(sum(Fy.init), digits=4)) != 1))
		stop('Initial frequencies must sum to 1')

	# unpack par.list
	gen  =  par.list$gen
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

	# Find equilibrium frequencies before introducing the inversion?
	if(eqStart) {
		initEq    <-  getEqFreqsSAY(par.list=par.list, eq.threshold=1e-6)
		aFreq.eq  <-  initEq$aFreq.eq
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
#	lambdaInv  <-  round(lambdaY(Wm=Wm, Xf=aFreq.eq[1], Y=aFreq.eq[3]), digits=5)
#	lambdaInv  <-  round(lambdaYWkSel2(Xf=aFreq.eq[1], Y=aFreq.eq[3], hm=hm, sm=sm), digits=6)
#	lambdaInv  <-  lambdaYWkSel2(Xf=aFreq.eq[1], Y=aFreq.eq[3], hm=hm, sm=sm)
	lambdaInv  <-  lambdaYWkSel(Fx=Fx.eq, Fy=Fy.eq, r=r, hm=hm, sm=sm)

	# Introduce inversion at low frequency
	invMut  <-  c(2,5)[as.vector(rmultinom(1,prob=Fy.eq[c(2,5)], size=1)) == 1]
	Fy.eq[invMut]  <-  Fy.eq[invMut] - 1e-5
	Fy.eq[invMut + 1]  <-  Fy.eq[invMut + 1] + 1e-5

	# Initialize .gen storage
	Fx.gen[1,]  <-  Fx.eq
	Fy.gen[1,]  <-  Fy.eq
	
	# save initial equilibrium frequencies
	aFreq.init  <-  aFreq.eq

	# Generation loop
	diffs  <-  rep(1,9)
	i=2
	while(i < (gen + 1) & any(abs(diffs[abs(diffs) > 0]) > eq.threshold)) {
		aFreq.eq  <-  c(Xf(Fx=Fx.gen[i-1,]), Xm(Fy=Fy.gen[i-1,], r=r), Y(Fy=Fy.gen[i-1,], r=r))
		Fxm  <-  round(c(x1m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
				   x2m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
				   x3m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r)), digits=8)
		Fym  <-  round(c(y1m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
				   y2cm( Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
				   y2cIm(Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
				   y2tm( Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
				   y3m(  Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r),
				   y3Im( Fx=Fx.gen[i-1,], Fy=Fy.gen[i-1,], r=r)), digits=8)
		FxPr  <-  round(c(x1Pr(  Fxm=Fxm, Wf=Wf),
					x2Pr(  Fxm=Fxm, Wf=Wf),
					x3Pr(  Fxm=Fxm, Wf=Wf)), digits=8)
		FyPr  <-  round(c(y1Pr(  Fym=Fym, Wm=Wm),
					y2cPr( Fym=Fym, Wm=Wm),
					y2cIPr(Fym=Fym, Wm=Wm),
					y2tPr( Fym=Fym, Wm=Wm),
					y3Pr(  Fym=Fym, Wm=Wm),
					y3IPr( Fym=Fym, Wm=Wm)), digits=8)
		Fxm         <-  Fxm/sum(Fxm)
		Fym         <-  Fym/sum(Fym)
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
	# Did the prediction based on the eigenvalue 
	# agree with the simulation results?
#browser()	
#	simInvasion  <-  sum(Fy.gen[nrow(Fy.gen),c(3,6)]) > 1e-05
	simInvasion  <-  sum(round(Fy.gen[nrow(Fy.gen),c(3,6)], digits=8)) > sum(round(Fy.gen[1,c(3,6)],digits=8))
	eigInvasion  <-  lambdaInv > 1

	#put names on aFreq.eq
	names(aFreq.eq)  <-  c("Xf","Xm","Y")

	# Put results in a list
	res  <-  list ("maxRuntime"   =  nrow(Fx.gen) >= gen,
				   "Fx.gen"       =  Fx.gen,
				   "Fy.gen"       =  Fy.gen,
				   "invAgree"     =  simInvasion == eigInvasion,
				   "eigInvasion"  =  eigInvasion,
				   "simInvasion"  =  simInvasion,
				   "lambdaInv"    =  lambdaInv,
				   "aFreq.init"   =  aFreq.init,
				   "aFreq.eq"     =  aFreq.eq,
				   "gen"          =  nrow(Fx.gen)
				   )

	# Return results
	return(res)
}




