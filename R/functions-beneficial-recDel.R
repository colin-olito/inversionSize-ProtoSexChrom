################################################################
#'  "Beneficial"
#' 
#'  Functions to perform deterministic and W-F simulations 
#'  for UNCONDITIONALLY BENEFICIAL Y-LINKED inversions 
#'  capturing the sex-determining locus, and n selected
#'  loci with segregating partially recessive deleterious
#'  mutations.
#'
#'
#'  Author: Colin Olito
#'
#'  NOTES:  






# Beneficial inversion expanding SLR w/ partially recessive deleterious mutations

##  EXACT RECURSIONS
wBarY  <-  function(sI, n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del) {
	pt.I.wt    <-  1 - qt.I.wt
	pt.Y.del   <-  1 - qt.Y.del
	pt.Y.wt    <-  1 - qt.Y.wt
		  YI.t*(1 + sI)*(( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )^r) * (1 - s*(h*(pt.I.wt*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))^(n-r) +
	(1 - YI.t)*( (1 - s*(h*((1 - XaOv.t.del)*qt.Y.del + XaOv.t.del*pt.Y.del) + XaOv.t.del*qt.Y.del)^r) * ((1 - s*(h*((1 - XaOv.t.wt)*qt.Y.wt + XaOv.t.wt*pt.Y.wt) + XaOv.t.wt*qt.Y.wt))^(n - r)) )
}

YI.prime  <-  function(sI, n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del, ...) {
	(YI.t*(1 + sI)*(( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )^r) * (1 - s*(h*((1 - qt.I.wt)*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))^(n-r)) / wBarY(sI=sI, n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qt.I.wt, qt.Y.wt=qt.Y.wt, qt.Y.del=qt.Y.del, XaOv.t.wt=XaOv.t.wt, XaOv.t.del=XaOv.t.del)
}

qI.wt.prime  <-  function(n, r, u, s, h, YI.t, qt.I.wt, XaOv.t.wt, XaOv.t.del) {
		(qt.I.wt*(1 - u - s*(XaOv.t.wt + u*(1 - 2*XaOv.t.wt) + h*(1 - XaOv.t.wt)*(1 - 2*u))) + u*(1 - s*(XaOv.t.wt + h*(1 - XaOv.t.wt)))) / 
			(1 - s*XaOv.t.wt*u - qt.I.wt*s*(XaOv.t.wt + u*(1 - 2*XaOv.t.wt)) - h*s*(qt.I.wt + XaOv.t.wt*(1 - 2*qt.I.wt) + u*(2 - 3*XaOv.t.wt - qt.I.wt*(3 - 4*XaOv.t.wt))))
}

qY.wt.prime  <-  function(u, s, h, qt.Y.wt, XaOv.t.wt) {
	(qt.Y.wt*(1 - u - s*(h + 2*XaOv.t.wt*(1 - h) + u*(2 - 3*h - 4*XaOv.t.wt*(1 - h)))) + XaOv.t.wt*(1 - u - s*(h + u*(2 - 3*h))) + 2*u*(1 - h*s)) / 
		(2*(1 - s*XaOv.t.wt*u - qt.Y.wt*s*(XaOv.t.wt + u*(1 - 2*XaOv.t.wt)) - h*s*(qt.Y.wt + XaOv.t.wt*(1 - 2*qt.Y.wt) + u*(2 - 3*XaOv.t.wt - qt.Y.wt*(3 - 4*XaOv.t.wt)))))
}

qY.del.prime  <-  function(u, s, h, qt.Y.del, XaOv.t.del) {
	(qt.Y.del*(1 - u - s*(h + 2*XaOv.t.del*(1 - h ) + u*(2 - 3*h - 4*XaOv.t.del*(1 - h)))) + XaOv.t.del*(1 - u - s*(h + u*(2 - 3*h))) + 2*u*(1 - h*s)) / 
		(2*(1 - s*XaOv.t.del*u - qt.Y.del*s*(XaOv.t.del + u*(1 - 2*XaOv.t.del)) - h*s*(qt.Y.del + XaOv.t.del*(1 - 2*qt.Y.del) + u*(2 - 3*XaOv.t.del - qt.Y.del*(3 - 4*XaOv.t.del)))))
}

XaOv.wt.prime  <-  function(u, s, h, XaOv.t.wt, XaSp.t.wt) {
	(2*u*(1 - h*s) + XaSp.t.wt*(1 - u - s*(h + u*(2 - 3*h))) + XaOv.t.wt*(1 - u - s*(h + 2*XaSp.t.wt*(1 - h) + u*(2 - 3*h - 4*XaSp.t.wt*(1 - h))))) / 
		(2*(1 - s*(XaSp.t.wt*u + XaOv.t.wt*(XaSp.t.wt + u*(1 - 2*XaSp.t.wt))) - h*s*(XaOv.t.wt + XaSp.t.wt*(1 - 2*XaOv.t.wt) + u*(2 - 3*XaSp.t.wt - XaOv.t.wt*(3 - 4*XaSp.t.wt)))))
}

XaOv.del.prime  <-  function(u, s, h, XaOv.t.del, XaSp.t.del) {
	(XaOv.t.del*(1 - u - s*(h + 2*XaSp.t.del*(1 - h) + u*(2 - 3*h - 4*XaSp.t.del*(1 - h)))) + XaSp.t.del*(1 - u - s*(2*u + h*(1 - 3*u))) + 2*u*(1 - h*s)) / 
		(2*(1 - s*XaSp.t.del*u - s*XaOv.t.del*(XaSp.t.del + u*(1 - 2*XaSp.t.del)) - h*s*(XaOv.t.del + XaSp.t.del*(1 - 2*XaOv.t.del) + u*(2 - 3*XaSp.t.del - XaOv.t.del*(3 - 4*XaSp.t.del)))))
}

XaSp.wt.prime  <-  function(u, s, h, YI.t, qt.I.wt, qt.Y.wt, XaOv.t.wt) {
(XaOv.t.wt*(1 + YI.t*(1 - 2*qt.I.wt*s) - h*s*(1 + YI.t*(1 - 2*qt.I.wt))) + u*(2 - 2*h*s - XaOv.t.wt*(1 + 2*s - 3*h*s) - YI.t*(XaOv.t.wt*(1 - h*s) + 2*(1 - h)*qt.I.wt*s*(1 - 2*XaOv.t.wt))) + qt.Y.wt*(1 - YI.t)*(1 - u - s*(h + XaOv.t.wt*(2 - 2*h) + u*(2 - 3*h - 4*XaOv.t.wt*(1 - h))))) / 
	(2 - s*XaOv.t.wt*(2*qt.Y.wt*(1 - YI.t) + 2*qt.I.wt*YI.t) - 2*h*s*(XaOv.t.wt + qt.Y.wt*(1 - 2*XaOv.t.wt) + YI.t*(qt.I.wt - qt.Y.wt)*(1 - 2*XaOv.t.wt)) - 2*s*u*(2*h + qt.Y.wt*(1 - 3*h) + XaOv.t.wt*(1 - 3*h - qt.Y.wt*(2 - 4*h)) + YI.t*(qt.I.wt - qt.Y.wt)*(1 - 2*XaOv.t.wt - h*(3 - 4*XaOv.t.wt))))
}

XaSp.del.prime  <-  function(u, s, h, YI.t, qt.Y.del, XaOv.t.del) {
(XaOv.t.del*(1 - h*s + qt.Y.del*YI.t*(1 - s*(2 - h))) + u*(2 - XaOv.t.del*(1 + qt.Y.del*YI.t) - h*s*(2 - 3*XaOv.t.del)*(1 - qt.Y.del*YI.t) - 2*s*(XaOv.t.del + qt.Y.del*YI.t*(1 - 2*XaOv.t.del))) + qt.Y.del*(1 - YI.t)*(1 - u - s*(h + 2*XaOv.t.del*(1 - h) + u*(2 - 3*h - 4*XaOv.t.del*(1 - h))))) / 
	(2 - 2*qt.Y.del*s*XaOv.t.del*(1 - YI.t) - 2*s*(qt.Y.del*XaOv.t.del*YI.t - h*(qt.Y.del + XaOv.t.del*(1 - 2*qt.Y.del) + YI.t*(qt.Y.del - qt.Y.del)*(1 - 2*XaOv.t.del))) - 2*s*u*(qt.Y.del + h*(2 - 3*qt.Y.del) + XaOv.t.del*(1 - 3*h - 2*qt.Y.del*(1 - 2*h)) + YI.t*(qt.Y.del - qt.Y.del)*(1 - 2*XaOv.t.del - h*(3 - 4*XaOv.t.del))))
}

eucDist <- function(x1, x2) {
  sqrt(sum((x1 - x2)^2))
} 

#################
#' Function to generate fig.1 data
findEq  <- function(sI, r, h, s, n, u, qHat, qHatDel) {

	# Empty Frequency Vectors
	YI.t        <- c()
	XaOv.wt.t   <- c()
	XaSp.wt.t   <- c()
	qI.wt.t     <- c()
	qY.wt.t     <- c()
	XaOv.del.t  <- c()
	XaSp.del.t  <- c()
	qY.del.t    <- c()

	# Find equilibrium prior to inversion 
	# First generation (all loci at equilibrium frequencies when inversion arises)
	YI.t[1]        <-  round(YI.prime(sI=sI, n=n, r=r, s=s, h=h, YI.t=0, qt.I.wt=0, qt.Y.wt=qHat, qt.Y.del=qHatDel, XaOv.t.wt=qHat, XaOv.t.del=qHatDel), digits=9)
	XaOv.wt.t[1]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=qHat, XaSp.t.wt=qHat), digits=9)
	XaSp.wt.t[1]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=0, qt.I.wt=0, qt.Y.wt=qHat, XaOv.t.wt=qHat), digits=9)
	qI.wt.t[1]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=0, qt.I.wt=0, XaOv.t.wt=qHat, XaOv.t.del=qHatDel), digits=9)
	qY.wt.t[1]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qHat, XaOv.t.wt=qHat), digits=9)
	XaOv.del.t[1]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=qHatDel, XaSp.t.del=qHatDel), digits=9)
	XaSp.del.t[1]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=0, qt.Y.del=qHatDel, XaOv.t.del=qHatDel), digits=9)
	qY.del.t[1]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qHatDel, XaOv.t.del=qHatDel), digits=9)

	# Subsequent generations
	i=2
	diff  <-  1
	while(diff > 1e-8) {
		YI.t[i]        <-  round(YI.prime(sI=sI, n=n, r=r, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		XaOv.wt.t[i]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t[i-1], XaSp.t.wt=XaSp.wt.t[i-1]), digits=9)
		XaSp.wt.t[i]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		XaOv.del.t[i]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t[i-1], XaSp.t.del=XaSp.del.t[i-1]), digits=9)
		XaSp.del.t[i]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qI.wt.t[i]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qY.wt.t[i]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		qY.del.t[i]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		diff  <-  eucDist(c(XaOv.wt.t[i-1],XaSp.wt.t[i-1],XaOv.del.t[i-1],XaSp.del.t[i-1],qY.wt.t[i-1],qY.del.t[i-1]),
						  c(XaOv.wt.t[i],XaSp.wt.t[i],XaOv.del.t[i],XaSp.del.t[i],qY.wt.t[i],qY.del.t[i]))
		i  <-  i + 1
	}
	res  <-  list(
				  "XaOv.wt.init" = XaOv.wt.t[i-1],
				  "XaSp.wt.init" = XaSp.wt.t[i-1],
				  "XaOv.del.init" = XaOv.del.t[i-1],
				  "XaSp.del.init" = XaSp.del.t[i-1],
				  "qI.wt.init" = qI.wt.t[i-1],
				  "qY.wt.init" = qY.wt.t[i-1],
				  "qY.del.init" = qY.del.t[i-1]
				  )
	return(res)
}






makeDeterministicFigSimData  <-  function(sI = 0.01, r = 0, x = 0.2, h = 0.1, s = 0.01, Ufactor = 2, generations = 10^4, ...) {
	# Parameters
	U     <-  Ufactor*s
	nTot  <-  10^4
	u     <-  U/nTot
	qHat  <-  (U/(nTot*h*s))
	if(r == 0) {
		qHatDel  <-  0
	} else {qHatDel  <-  qHat}
	n     <-  nTot*x
	N     <-  5000
	YI.0  <-  2/N

	# Empty Frequency Vectors
	YI.t        <- c()
	XaOv.wt.t   <- c()
	XaSp.wt.t   <- c()
	qI.wt.t     <- c()
	qY.wt.t     <- c()
	XaOv.del.t  <- c()
	XaSp.del.t  <- c()
	qY.del.t    <- c()
	wbarYI.t    <- c()

	# Find equilibrium prior to inversion 
	eqs  <-  findEq(sI=sI, r=r, h=h, s=s, n=n, u=u, qHat=qHat, qHatDel=qHatDel)
	
	# First generation (all loci at equilibrium frequencies when inversion arises)
	YI.t[1]        <-  round(YI.prime(sI=sI, n=n, r=r, s=s, h=h, YI.t=YI.0, qt.I.wt=0, qt.Y.wt=eqs$qY.wt.init, qt.Y.del=eqs$qY.del.init, XaOv.t.wt=eqs$XaOv.wt.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
	XaOv.wt.t[1]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=eqs$XaOv.wt.init, XaSp.t.wt=eqs$XaSp.wt.init), digits=9)
	XaSp.wt.t[1]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.0, qt.I.wt=0, qt.Y.wt=eqs$qY.wt.init, XaOv.t.wt=eqs$XaOv.wt.init), digits=9)
	qI.wt.t[1]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.0, qt.I.wt=0, XaOv.t.wt=eqs$XaOv.wt.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
	qY.wt.t[1]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=eqs$qY.wt.init, XaOv.t.wt=eqs$XaOv.wt.init), digits=9)
	XaOv.del.t[1]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=eqs$XaOv.del.init, XaSp.t.del=eqs$XaSp.del.init), digits=9)
	XaSp.del.t[1]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.0, qt.Y.del=eqs$qY.del.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
	qY.del.t[1]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=eqs$qY.del.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
	wbarYI.t[1]    <-  round((YI.t[1]/YI.0), digits=9)
	# Subsequent generations
	i=2
	while(i < generations+1) {
		YI.t[i]        <-  round(YI.prime(sI=sI, n=n, r=r, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		XaOv.wt.t[i]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t[i-1], XaSp.t.wt=XaSp.wt.t[i-1]), digits=9)
		XaSp.wt.t[i]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		XaOv.del.t[i]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t[i-1], XaSp.t.del=XaSp.del.t[i-1]), digits=9)
		XaSp.del.t[i]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qI.wt.t[i]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qY.wt.t[i]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		qY.del.t[i]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		wbarYI.t[i]    <-  round((YI.t[i]/ YI.t[i-1]), digits=9)
		i  <-  i + 1
	}
	# Tack on initial freq.
	YI.t        <- c(YI.0,YI.t)
	XaOv.wt.t   <-  c(eqs$XaOv.wt.init, XaOv.wt.t)
	XaSp.wt.t   <-  c(eqs$XaSp.wt.init, XaSp.wt.t)
	XaOv.del.t  <-  c(eqs$XaOv.del.init, XaOv.del.t)
	XaSp.del.t  <-  c(eqs$XaSp.del.init, XaSp.del.t)
	qI.wt.t     <-  c(0, qI.wt.t)
	qY.wt.t     <-  c(eqs$qY.wt.init, qY.wt.t)
	qY.del.t    <-  c(eqs$qY.del.init, qY.del.t)
	wbarYI.t    <-  c(NA, wbarYI.t)

	# Return results as df
	results  <-  data.frame(
							"YI.t"        =  YI.t,
							"XaOv.wt.t"   =  XaOv.wt.t,
							"XaSp.wt.t"   =  XaSp.wt.t,
							"XaOv.del.t"  =  XaOv.del.t,
							"XaSp.del.t"  =  XaSp.del.t,
							"qI.wt.t"     =  qI.wt.t,
							"qY.wt.t"     =  qY.wt.t,
							"qY.del.t"    =  qY.del.t,
							"wbar.YI.t"   =  wbarYI.t
							)
	return(results)
}






########################################
## Multilocus Wright-Fisher Simulations
## Using exact recursions derived in Mathematica using Xf, Xm, Y
## inluding sampling variance for XaOv, XaSp, etc.
## 
wBarY.multi  <-  function(sI, n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del) {
	pt.I.wt    <-  1 - qt.I.wt
	pt.Y.del   <-  1 - qt.Y.del
	pt.Y.wt    <-  1 - qt.Y.wt
		  YI.t*(1 + sI)*(prod( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )) * prod(1 - s*(h*(pt.I.wt*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt)) +
	(1 - YI.t)*( prod(1 - s*(h*((1 - XaOv.t.del)*qt.Y.del + XaOv.t.del*pt.Y.del) + XaOv.t.del*qt.Y.del)) * prod((1 - s*(h*((1 - XaOv.t.wt)*qt.Y.wt + XaOv.t.wt*pt.Y.wt) + XaOv.t.wt*qt.Y.wt))) )
}

YI.multi.prime  <-  function(sI, n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del, ...) {
	(YI.t*(1 + sI)*(prod( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del))) * prod(1 - s*(h*((1 - qt.I.wt)*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))) / wBarY.multi(sI=sI, n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qt.I.wt, qt.Y.wt=qt.Y.wt, qt.Y.del=qt.Y.del, XaOv.t.wt=XaOv.t.wt, XaOv.t.del=XaOv.t.del)
}

# Function to estimate Pr(fix | x) for different inversion sizes (x)
# Looping over different Population size and U/s ratios
makeDataPrFixInvSize  <-  function(sI = 0.01, h = 0.1, s = 0.01, Us.factor.vals = c(2, 5, 10),
																	 nTot = 10^4, N.vals = c(10^3, 10^4), Nfname = "") {

	# Containers
	PrFix       <-  c()
	YI.t        <-  c()
	XaOv.wt.t   <-  c()
	XaSp.wt.t   <-  c()
	qI.wt.t     <-  c()
	qY.wt.t     <-  c()
	XaOv.del.t  <-  c()
	XaSp.del.t  <-  c()
	qY.del.t    <-  c()

	# inversion sizes
	invSize  <-  c(0.5,1:9)/10

	# Loop over population size
	for(i in 1:length(N.vals)) {
		N  <-  N.vals[i]
		YI.0   <-  2/N

		#number of simulations 
		sims  = 100*N/2

		# Loop over Us factor values
		for(j in 1:length(Us.factor.vals)) {
			U     <-  Us.factor.vals[j]*s
			u     <-  U/nTot
			qHat  <-  (U/(nTot*h*s))

			# Loop over inversion size
			for(k in 1:length(invSize)) {
				# # loci captured by inversion
				n   <-  nTot*invSize[k]

				# counter for fixations
				fix   = 0

				# Loop over replicate simulations
				for(l in 1:sims){


					# Draw random value for # del. mutations captured by inversion (r) given x
					r   <-  rpois(1, lambda=(U*invSize[k]/(s*h)))

					# Assign frequencies
					# Note implicit assumption of equal initial
					# frequencies in XOv, XSp, Y chromosomes
					YI.t        <- YI.0
					XaOv.wt.t   <- qHat
					XaSp.wt.t   <- qHat
					qI.wt.t     <- 0
					qY.wt.t     <- qHat
					XaOv.del.t  <- qHat
					XaSp.del.t  <- qHat
					qY.del.t    <- qHat
# t=1
# plot(NA, ylim=c(0,1.001), xlim=c(0,N))
# abline(h=1)
# text(x=(N/2), y=0.9, labels=paste(l, sep=''))
					# Run forward simulation
					while(YI.t*(1 - YI.t) > 0) {
						#expected frequency after selection
						XaOv.wt.t   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t), digits=9)
						XaOv.del.t  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t), digits=9)
						XaSp.wt.t   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						XaSp.del.t  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qY.wt.t     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
						qY.del.t    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
						qI.wt.t     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
						YI.sel      <-  round(YI.prime(sI=sI, n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, qt.Y.del=qY.del.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
#points(YI.sel/YI.t ~ t, col=2)

						#binomial sampling
						YI.t        <-  rbinom(1, (N/2), YI.sel)/(N/2)

#points(YI.t ~ t)
#points(qI.wt.t ~ t, col=2)
#t=t+1
					}
					fix  <-  fix + YI.t
				}
			PrFix  <-  c(PrFix, (fix/sims))
			cat('\r', paste("N: ", i, "/", length(N.vals), 
										", U: ", j, "/", length(Us.factor.vals), 
										", x: ", round(100*(k/length(invSize))), "% complete", sep=""))
			}

		}

	}

	# Index variables
	Ns        <-  rep(N.vals, each=(length(Us.factor.vals)*length(invSize)))
	Us        <-  rep((Us.factor.vals*s), each=(length(invSize)), times=length(N.vals))
	Ufac      <-  rep((Us.factor.vals), each=(length(invSize)), times=length(N.vals))
	invSizes  <-  rep(invSize, times=length(N.vals)*length(Us.factor.vals))
	
	# Export Results Dataframe
	filename  <-  paste("./output/data/simResults/beneficial_PrFixFig_sI", sI,"_h", h, "_s", s, Nfname, ".csv", sep="")
	d  <-  data.frame(
										"sI"     =  rep(sI, times=length(PrFix)),
										"h"      =  rep(h, times=length(PrFix)),
										"s"      =  rep(s, times=length(PrFix)),
										"N"      =  Ns,
										"U"      =  Us,
										"Ufac"   =  Ufac,
										"x"      =  invSizes,
										"PrFix"  =  PrFix
										)
	write.csv(d, file=filename, row.names=FALSE)

}



