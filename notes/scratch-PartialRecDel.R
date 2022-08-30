qHat  <-  function(n, U, h, s) {
	U/(n*h*(s))
}

WBar.r  <-  function(n, r, h, s, q) {
	((1 - h*s*(1 - q) - s*q)^r)*(1 - h*s*q)^(n-r)
}

WBar.r.approx  <-  function(n, r, h, s, q) {
	exp(-r*s*(h*(1 - q) + q) - (n-r)*h*s*q)
}

n  <-  1000
r  <-  seq(0,200)
h  <-  0.25
s  <-  0.01
U  <-  0.1
q  <-  qHat(n=1000, U=U, h=h, s=s)

real  <-  WBar.r(n=n, r=r, h=h, s=s, q=q)
approx  <-  WBar.r.approx(n=n, r=r, h=h, s=s, q=q)

plot(real ~ seq_along(real))
points(approx ~ seq_along(approx), pch=9)





##################################
# Approx. Deterministic trajectory 
U    <-  0.1
h    <-  0.25
s    <-  0.1
n    <-  1000
r    <-  0
t    <-  1000
x    <-  0.1
YInext  <-  function(YIt, U, h, s, n, r, t, x) {
	qHat  <-  (U/(n*h*s))
	qt    <-  ((1 - exp((h*(-1 + qHat) - qHat)*s*t))*U)/(n*(h + qHat - h*qHat)*s)
	YIt*exp(-s*(r*(h*(1 - qHat) - qHat) + (n*x - r)*(qt*(h + qHat) + h*qHat*(1 - 2*qt)))) / ((1-YIt)*exp(-U*x*(2 + U/(n*(h^2)*s))) + YIt*exp(-s*(r*(h*(1 - qHat) - qHat) + (n*x - r)*(qt*(h + qHat) + h*qHat*(1 - 2*qt)))))
}


# Parameters
YIt   <-  c()
Ud    <-  0.1
hd    <-  0.1
sd    <-  0.01
n     <-  2*10^6
r     <-  5
qHat  <-  (Ud/(n*hd*sd))
n*qHat
x     <-  0.5
N     <-  1000
t     <-  c(1:9999)
YI0   <-  2/N

YIt  <- c()
for(i in 1:length(t)) {
	if(i == 1) {
		YIt[i]  <-  YInext(YIt=YI0, U=Ud, h=hd, s=sd, n=n, r=r, t=t[i], x=x)
	}else{YIt[i]  <-  YInext(YIt = YIt[i-1], U=Ud, h=hd, s=sd, n=n, r=r, t=t[i], x=x) }
	if(YIt[i] == 1) {
		YIt[i]  <-  YIt[i] - 0.00001
	}
}
YIt  <- c(YI0,YIt)
plot(YIt ~ seq_along(YIt), ylim=c(0,1), type='l', lwd=2)









##################################
# Pseudo Deterministic W-F Simulations 
rbernoulli <- function(n,p) runif(n) < p


n   <-  2000000
mu  <-  10^-5
h   <-  0.1
qHat <-  c()
mut  <-  c()
for(i in 1:n) {
	s  <-  rnorm(n=1,mean=0.01, sd=0.001)
	qHat[i]  <-  (mu/(h*s))
	mut[i]   <-  rbernoulli(n=1,p=qHat[i])
}

sum(mut)
sum(qHat)
n*(mu/h*0.01)
sum(rpois(n=2000000,lambda=(mu/(h*0.01))))


# Parameters
U    <-  0.1
h    <-  0.1
s    <-  0.01
n     <-  2*10^6
r     <-  10
qHat  <-  (U/(n*h*s))
n*qHat
x     <-  0.1
N     <-  1000
YI0   <-  2/N

#number of simulations 
sims  = 100*N/2
fix   = 0
for(i in 1:sims){
	YI  <-  YI0
	t   <-  1
	while(YI*(1 - YI) > 0) {

		#expected frequency after selection          
		qt    <-  ((1 - exp((h*(-1 + qHat) - qHat)*s*t))*U)/(n*(h + qHat - h*qHat)*s)
		F.sel = YI*exp(-s*(r*(h*(1 - qHat) - qHat) + (n*x - r)*(qt*(h + qHat) + h*qHat*(1 - 2*qt)))) / ((1-YI)*exp(-U*x*(2 + U/(n*(h^2)*s))) + YI*exp(-s*(r*(h*(1 - qHat) - qHat) + (n*x - r)*(qt*(h + qHat) + h*qHat*(1 - 2*qt)))))

        #binomial sampling
        YI = rbinom(1, (N/2), F.sel)/(N/2)
        t = t + 1
    }
    if(YI == 1){
    	fix  <-  fix + 1
    }
    cat('\r', paste("Progress:", round(100*(i/sims)), "% complete"))
}    

fix/sims





##################################
# Pseudo Deterministic W-F Simulations 
# loop over r to get a profile of Pr(fix)



# Parameters
n     <-  10*10^6
mu    <-  10^-8
U     <-  n*mu
h     <-  0.01
s     <-  0.01
NrR   <-  rpois(n*x,lambda=U*x/(h*s))
r     <-  c(range(NrR)[1]:range(NrR)[2])
qHat  <-  mu/(h*s)
x     <-  0.2
n*x*qHat
N     <-  1000
YI0   <-  2/N
PrFix <-  c()

#number of simulations 
sims  = 100*N/2

# Randomly drawn number of mutations

for(j in 1:length(r)) {
fix   = 0
	for(i in 1:sims){
		YI  <-  YI0
		t   <-  0
		while(YI*(1 - YI) > 0) {

			#expected frequency after selection          
			qt    <-  ((1 - exp((h*(-1 + qHat) - qHat)*s*t))*U)/(n*(h + qHat - h*qHat)*s)
			F.sel = YI*exp(-s*(r[j]*(h*(1 - qHat) - qHat) + (n*x - r[j])*(qt*(h + qHat) + h*qHat*(1 - 2*qt)))) / ((1-YI)*exp(-U*x*(2 + U/(n*(h^2)*s))) + YI*exp(-s*(r[j]*(h*(1 - qHat) - qHat) + (n*x - r[j])*(qt*(h + qHat) + h*qHat*(1 - 2*qt)))))

        	#binomial sampling
        	YI = rbinom(1, (N/2), F.sel)/(N/2)
        	t = t + 1
    	}
    	if(YI == 1){
    	fix  <-  fix + 1
    	}
	}    
	PrFix[j]  <-  fix/sims
	cat('\r', paste("Progress:", round(100*(j/length(r))), "% complete"))
}

PrR  <-  dpois(r,lambda=U*x/(h*s))
par(mfrow=c(2,2))
plot(PrFix ~ r)
plot(PrR ~ r)
plot((PrFix*PrR) ~ r)
plot(PrR*10000)

sum(PrR*10000)
sum(PrR*10000*PrFix)


PrFix[r==n*x*qHat*0.8]*10000
10000*PrFix[r == (nTot*x*qHat*0.8)]


rm(list=ls())

##########################
##  MORE EXACT RECURSIONS
wBar.i.wt  <-  function(s, h, qNt.wt, qIt.wt, YIt) {
	pNt.wt  <-  1  - qNt.wt
	pIt.wt  <-  YIt - qIt.wt
	1 - qNt.wt*s*(h*(pNt.wt + ((1 - YIt)*pIt.wt + YIt*pNt.wt)/2) + (qNt.wt + (1 - YIt)*qIt.wt + YIt*qNt.wt)/2) - pNt.wt*s*h*((qNt.wt + (1 - YIt)*qNt.wt + YIt*qIt.wt)/2)
}
wBar.i.del  <-  function(s, h, qNt.del, YIt) {
	pNt.del  <-  1  - qNt.del
	1 - qNt.del*s*(h*((pNt.del + (1 - YIt)*pNt.del)/2) + (qNt.del + (1 - YIt)*pNt.del + YIt)/2) - pNt.del*s*h*((qNt.del + (1 - YIt)*qNt.del + YIt)/2)
}
wBarI  <-  function(n, r, s, h, YIt, qNt.del, qNt.wt, qIt.wt) {
	QNt.del  <-  qNt.del#/(1 - YIt)
	QNt.wt   <-  qNt.wt#/(1 - YIt)
	QIt.wt   <-  qIt.wt/YIt
		  YIt*(( 1 - s*(h*(1 - QNt.del) + QNt.del) )^r) * (1 - s*(h*((1 - QIt.wt)*QNt.wt + QIt.wt*(1 - QNt.wt)) + QIt.wt*QNt.wt))^(n-r) +
	(1 - YIt)*( ((1 - s*(2*h*QNt.del*(1 - QNt.del) + QNt.del^2))^r) * ((1 - s*(2*h*QNt.wt*(1 - QNt.wt) + QNt.wt^2))^(n - r)) )
}

YI.next  <-  function(n, r, s, h, YIt, qNt.del, qNt.wt, qIt.wt, ...) {
	QNt.del  <-  qNt.del#/(1 - YIt)
	QNt.wt   <-  qNt.wt#/(1 - YIt)
	QIt.wt   <-  qIt.wt/YIt
	(YIt*(( 1 - s*(h*(1 - QNt.del) + QNt.del) )^r) * (1 - s*(h*((1 - QIt.wt)*QNt.wt + QIt.wt*(1 - QNt.wt)) + QIt.wt*QNt.wt))^(n-r)) / wBarI(n=n, r=r, s=s, h=h, YIt=YIt, qNt.del=qNt.del, qNt.wt=qNt.wt, qIt.wt=qIt.wt)
}

qI.wt.next  <-  function(n, r, u, s, h, YIt, qNt.del, qNt.wt, qIt.wt) {
	pNt.wt  <-  1 - qNt.wt
	pIt.wt  <-  YIt - qIt.wt
	YI.next(n=n, r=r, s=s, h=h, YIt=YIt, qNt.del=qNt.del, qIt.wt=qIt.wt, qNt.wt=qNt.wt) *
		( (qIt.wt*(1 - s*(h*pNt.wt + qNt.wt)) + u*pIt.wt*wBar.i.wt(s=s, h=h, qNt.wt=qNt.wt, qIt.wt=qIt.wt, YIt=YIt)) / ( YIt - qIt.wt*s*(h*pNt.wt + qNt.wt) ) )
}

qN.wt.next  <-  function(n, r, u, s, h, YIt, qNt.del, qNt.wt, qIt.wt) {
	pNt.wt  <-  1 -  qNt.wt
	pIt.wt  <-  YIt - qIt.wt
#	(1 - YI.next(n=n, r=r, s=s, h=h, YIt=YIt, qNt.del=qNt.del, qIt.wt=qIt.wt, qNt.wt=qNt.wt)/4)*
		( qNt.wt*(1 - s*(h*((1 - YIt)*pNt.wt + YIt*pIt.wt)/2) + (qNt.wt + (1 - YIt)*qNt.wt + YIt*(qIt.wt))/2) + u*pNt.wt*wBar.i.wt(s=s, h=h, qNt.wt=qNt.wt, qIt.wt=qIt.wt, YIt=YIt) ) / 
			( 1  - qNt.wt*s*(h*((1 - YIt)*pNt.wt + YIt*pIt.wt)/2) + (qNt.wt + (1 - YIt)*qNt.wt + YIt*(qIt.wt))/2) 
}

qN.del.next  <-  function(n, r, u, s, h, YIt, qNt.del, qNt.wt, qIt.wt) {
	pNt.del  <-  1 - qNt.del
#	(1 - YI.next(n=n, r=r, s=s, h=h, YIt=YIt, qNt.del=qNt.del, qIt.wt=qIt.wt, qNt.wt=qNt.wt)/4)*
		(qNt.del*(1 - s*(h*(pNt.del + 2*(1 - YIt)*pNt.del/2) + (qNt.del + 2*(1 - YIt)*qNt.del + YIt)/2)) + ((u*pNt.del*2*(1 - YIt)*pNt.del)/2)*wBar.i.del(s=s, h=h, qNt.del=qNt.del, YIt=YIt)) / 
			(1  - qNt.del*s*(h*(pNt.del + 2*(1 - YIt)*pNt.del/2) + (qNt.del + 2*(1 - YIt)*qNt.del + YIt)/2))
}



# Parameters
h    <-  0.1
s    <-  0.02
U    <-  10*s
nTot  <-  10000
r     <-  0
u     <-  U/nTot
qHat  <-  (U/(nTot*h*s))
x     <-  0.2
n      <-  nTot*x
N     <-  5000
YI0   <-  2/N

YI.t      <- c()
qN.del.t  <- c()
qN.wt.t   <- c()
qI.wt.t   <- c()

# First generation
	YI.t[1]      <-  round(YI.next(n=n, r=r, s=s, h=h, YIt=YI0, qNt.del=qHat, qNt.wt=qHat, qIt.wt=0), digits=9)
	qN.del.t[1]  <-  round(qN.del.next(n=n, r=r, u=u, s=s, h=h, YIt=YI0, qNt.del=qHat, qNt.wt=qHat, qIt.wt=0), digits=9)
	qN.wt.t[1]   <-  round(qN.wt.next(n=n, r=r, u=u, s=s, h=h, YIt=YI0, qNt.del=qHat, qNt.wt=qHat, qIt.wt=0), digits=9)
	qI.wt.t[1]   <-  round(qI.wt.next(n=n, r=r, u=u, s=s, h=h, YIt=YI0, qNt.del=qHat, qNt.wt=qHat, qIt.wt=0), digits=9)
# Subsequent generations
i=2
#while(YI.t[i-1]*(1 - YI.t[i-1]) > 0) {
while(i < 5001) {
		YI.t[i]      <-  round(YI.next(n=n, r=r, s=s, h=h, YIt=YI.t[i-1], qNt.del=qN.del.t[i-1], qNt.wt=qN.wt.t[i-1], qIt.wt=qI.wt.t[i-1]), digits=9)
		qN.del.t[i]  <-  round(qN.del.next(n=n, r=r, u=u, s=s, h=h, YIt=YI.t[i-1], qNt.del=qN.del.t[i-1], qNt.wt=qN.wt.t[i-1], qIt.wt=qI.wt.t[i-1]), digits=9)
		qN.wt.t[i]   <-  round(qN.wt.next(n=n, r=r, u=u, s=s, h=h, YIt=YI.t[i-1], qNt.del=qN.del.t[i-1], qNt.wt=qN.wt.t[i-1], qIt.wt=qI.wt.t[i-1]), digits=9)
		qI.wt.t[i]   <-  round(qI.wt.next(n=n, r=r, u=u, s=s, h=h, YIt=YI.t[i-1], qNt.del=qN.del.t[i-1], qNt.wt=qN.wt.t[i-1], qIt.wt=qI.wt.t[i-1]), digits=9)
		i  <-  i + 1
}
YI.t  <- c(YI0,YI.t)
par(mfrow=c(1,2))
plot(YI.t ~ seq_along(YI.t), ylim=c(0,1), type='l', lwd=2, col=2)
lines(qI.wt.t, lwd=2, lty=1, col=2)
lines(qN.wt.t, lwd=2, lty=1)
if(r > 0){lines(qN.del.t, lwd=2, lty=2)}
legend("topleft",
		legend = c(expression(Y[t]^I), 
				   expression(q[I.wt]), 
				   expression(q[N.wt]), 
				   expression(q[N.del])),
		col = c(2, 2, 1, 1),
		lty = c(1,1,1,2),
		lwd=2
		)
plot(qI.wt.t ~ seq_along(qI.wt.t), ylim=c(0,max(c(qI.wt.t,qN.wt.t,qN.del.t))), type='l', lty=1, lwd=2, col=2)
lines(qN.wt.t, lwd=2, lty=1)
if(r > 0){lines(qN.del.t, lwd=2, lty=2)}
legend("right",
		legend = c(expression(q[I.wt]), 
				   expression(q[N.wt]), 
				   expression(q[N.del])),
		col = c(2, 1, 1),
		lty = c(1,1,2),
		lwd=2
		)
nTot*x
nTot*x*qHat
qHat
r








############
# Trying again with recursions derived in Mathematica using Xf, Xm, Y



rm(list=ls())

##########################
##  MORE EXACT RECURSIONS
wBarY  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del) {
	pt.I.wt    <-  1 - qt.I.wt
	pt.Y.del   <-  1 - qt.Y.del
	pt.Y.wt    <-  1 - qt.Y.wt
		  YI.t*(( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )^r) * (1 - s*(h*(pt.I.wt*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))^(n-r) +
	(1 - YI.t)*( (1 - s*(h*((1 - XaOv.t.del)*qt.Y.del + XaOv.t.del*pt.Y.del) + XaOv.t.del*qt.Y.del)^r) * ((1 - s*(h*((1 - XaOv.t.wt)*qt.Y.wt + XaOv.t.wt*pt.Y.wt) + XaOv.t.wt*qt.Y.wt))^(n - r)) )
}

YI.prime  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del, ...) {
	(YI.t*(( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )^r) * (1 - s*(h*((1 - qt.I.wt)*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))^(n-r)) / wBarY(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qt.I.wt, qt.Y.wt=qt.Y.wt, qt.Y.del=qt.Y.del, XaOv.t.wt=XaOv.t.wt, XaOv.t.del=XaOv.t.del)
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


# Parameters
h    <-  0.25
s    <-  0.01
U    <-  10*s
nTot  <-  10000
r     <-  2
u     <-  U/nTot
qHat  <-  (U/(nTot*h*s))
if(r == 0) {
	qHatDel  <-  0
} else {qHatDel  <-  qHat}
x     <-  0.7
n      <-  nTot*x
N     <-  5000
YI.0   <-  2/N

YI.t        <- c()
XaOv.wt.t   <- c()
XaSp.wt.t   <- c()
qI.wt.t     <- c()
qY.wt.t     <- c()
XaOv.del.t  <- c()
XaSp.del.t  <- c()
qY.del.t    <- c()

# First generation
	YI.t[1]        <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=YI.0, qt.I.wt=0, qt.Y.wt=qHat, qt.Y.del=qHatDel, XaOv.t.wt=qHat, XaOv.t.del=qHatDel), digits=9)
	XaOv.wt.t[1]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=qHat, XaSp.t.wt=qHat), digits=9)
	XaSp.wt.t[1]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.0, qt.I.wt=0, qt.Y.wt=qHat, XaOv.t.wt=qHat), digits=9)
	qI.wt.t[1]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.0, qt.I.wt=0, XaOv.t.wt=qHat, XaOv.t.del=qHatDel), digits=9)
	qY.wt.t[1]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qHat, XaOv.t.wt=qHat), digits=9)
	XaOv.del.t[1]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=qHatDel, XaSp.t.del=qHatDel), digits=9)
	XaSp.del.t[1]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.0, qt.Y.del=qHatDel, XaOv.t.del=qHatDel), digits=9)
	qY.del.t[1]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qHatDel, XaOv.t.del=qHatDel), digits=9)
# Subsequent generations
i=2
#while(YI.t[i-1]*(1 - YI.t[i-1]) > 0) {
while(i < 10001) {
		YI.t[i]        <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		XaOv.wt.t[i]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t[i-1], XaSp.t.wt=XaSp.wt.t[i-1]), digits=9)
		XaSp.wt.t[i]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		XaOv.del.t[i]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t[i-1], XaSp.t.del=XaSp.del.t[i-1]), digits=9)
		XaSp.del.t[i]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qI.wt.t[i]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
		qY.wt.t[i]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
		qY.del.t[i]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
	i  <-  i + 1
}
YI.t  <- c(YI.0,YI.t)
par(mfrow=c(1,2))
plot(YI.t ~ seq_along(YI.t), ylim=c(0,1), main="Inversion captures 0 del. mut. (r = 0)", ylab="Frequency", xlab="Generation", type='l', lwd=2, col=2)
legend("topright",
		legend = c(expression(Y[t]^I)),
		col = c(2),
		lty = c(1),
		lwd=2
		)
plot(qI.wt.t ~ seq_along(qI.wt.t), ylim=c(0,max(c(qI.wt.t,XaOv.wt.t,XaOv.del.t))), ylab="Frequency", xlab="Generation", type='l', lty=1, lwd=2, col=2)
lines(qY.wt.t, lwd=2, lty=1, col=1)
lines(XaOv.wt.t, lwd=2, lty=1, col=2)
lines(XaSp.wt.t, lwd=2, lty=1, col=4)
if(r > 0){
	lines(qY.del.t, lwd=2, lty=1, col=1)
	lines(XaOv.del.t, lwd=2, lty=2)
	lines(XaSp.del.t, lwd=2, lty=2, col=4)
}
legend("bottomright",
		legend = c(expression(q[I.wt]), 
				   expression(q[Y.wt]), 
				   expression(Xa[Ov.wt]), 
				   expression(Xa[Sp.wt]), 
				   expression(q[Y.del]), 
				   expression(Xa[Ov.del]), 
				   expression(Xa[Sp.del])),
		col = c(2,1,2,4,1,2,4),
		lty = c(1,1,1,1,2,2,2),
		lwd=2
		)

min(XaOv.del.t)
u/s

sqrt(u/s)
nTot*x
nTot*x*qHat
qHat
r



########################################
## Wright-Fisher using Exact recursions derived in Mathematica using Xf, Xm, Y
# Parameters
h    <-  0.1
s    <-  0.01
U    <-  10*s
nTot  <-  10000
NrR   <-  rpois(n*x,lambda=U*x/(h*s))
#r     <-  c(range(NrR)[1]:range(NrR)[2])
r     <-  c(0:range(NrR)[2])
#r     <-  c(0:5)
u     <-  U/nTot
qHat  <-  (U/(nTot*h*s))
x     <-  0.7
n     <-  nTot*x
N     <-  1000
YI.0   <-  2/N
PrFix <-  c()


YI.t        <- c()
XaOv.wt.t   <- c()
XaSp.wt.t   <- c()
qI.wt.t     <- c()
qY.wt.t     <- c()
XaOv.del.t  <- c()
XaSp.del.t  <- c()
qY.del.t    <- c()


#number of simulations 
sims  = 100*N/2

# Loop over # deleterious mutations initially captured by inversion
for(j in 1:length(r)) {
fix   = 0
	
	# Loop over replicate simulations
	for(i in 1:sims){
		YI.t        <- YI.0
		XaOv.wt.t   <- qHat
		XaSp.wt.t   <- qHat
		qI.wt.t     <- 0
		qY.wt.t     <- qHat
		XaOv.del.t  <- qHat
		XaSp.del.t  <- qHat
		qY.del.t    <- qHat

		# Run forward simulation
		while(YI.t*(1 - YI.t) > 0) {

			#expected frequency after selection
			XaOv.wt.t   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t), digits=9)
			XaSp.wt.t   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
			XaOv.del.t  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t), digits=9)
			XaSp.del.t  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
			qI.wt.t     <-  round(qI.wt.prime(n=n, r=r[j], u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
			qY.wt.t     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
			qY.del.t    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
			YI.sel      <-  round(YI.prime(n=n, r=r[j], s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, qt.Y.del=qY.del.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)

			#binomial sampling
			YI.t      <-  rbinom(1, (N/2), YI.sel)/(N/2)
#			qN.del  <-  rbinom(1, (2*N - (YI*(N/2))), qN.del.sel)/(2*N - (YI*(N/2)))
#			qN.wt   <-  rbinom(1, (2*N - (YI*(N/2))), qN.wt.sel)/(2*N - (YI*(N/2)))
#			qI.wt   <-  rbinom(1, (YI*(N/2)), qI.wt.sel)/(N/2)
		}
		if(YI.t == 1){
		fix  <-  fix + 1
		}
	}
	PrFix[j]  <-  fix/sims
	cat('\r', paste("Progress:", round(100*(j/length(r))), "% complete"))
}

PrR  <-  dpois(r,lambda=U*x/(h*s))
par(mfrow=c(2,3))
plot(PrFix ~ r)
plot(PrR ~ r)
plot((PrFix*PrR) ~ r)
plot(PrR*nTot)
plot((PrR*nTot*PrFix) ~ r)


sum(PrR*nTot)
sum(PrR*nTot*PrFix)

nTot*x*qHat*0.8
nTot*PrFix[r == (nTot*x*qHat*0.8)]




########################################
## Wright-Fisher using Exact recursions derived in Mathematica using Xf, Xm, Y
## inluding sampling variance for XaOv, XaSp, etc.
# Parameters
h    <-  0.1
s    <-  0.01
U    <-  10*s
x     <-  0.2
nTot  <-  10000
NrR   <-  rpois(500,lambda=nTot*x*qHat)
r     <-  seq(0,range(NrR)[2], by=4)
u     <-  U/nTot
qHat  <-  (U/(nTot*h*s))
n     <-  nTot*x
N     <-  1000
YI.0   <-  2/N
PrFix <-  c()


YI.t        <- c()
XaOv.wt.t   <- c()
XaSp.wt.t   <- c()
qI.wt.t     <- c()
qY.wt.t     <- c()
XaOv.del.t  <- c()
XaSp.del.t  <- c()
qY.del.t    <- c()


#number of simulations 
sims  = 100*N/2

# Loop over # deleterious mutations initially captured by inversion
for(j in 1:length(r)) {
fix   = 0
	
	# Loop over replicate simulations
	for(i in 1:sims){
		YI.t        <- YI.0
		XaOv.wt.t   <- qHat
		XaSp.wt.t   <- qHat
		qI.wt.t     <- 0
		qY.wt.t     <- qHat
		XaOv.del.t  <- qHat
		XaSp.del.t  <- qHat
		qY.del.t    <- qHat

		# Run forward simulation
		while(YI.t*(1 - YI.t) > 0) {

			#expected frequency after selection
			XaOv.wt.sel   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t), digits=9)
			XaOv.del.sel  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t), digits=9)
			XaSp.wt.sel   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
			XaSp.del.sel  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
			qY.wt.sel     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
			qY.del.sel    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
			qI.wt.sel     <-  round(qI.wt.prime(n=n, r=r[j], u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
			YI.sel        <-  round(YI.prime(n=n, r=r[j], s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, qt.Y.del=qY.del.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)

			#binomial sampling
			XaOv.wt.t   <-  rbinom(1, N, XaOv.wt.sel)/N
			XaOv.del.t  <-  rbinom(1, N, XaOv.del.sel)/N
			YI.t        <-  rbinom(1, (N/2), YI.sel)/(N/2)
			qY.wt.t     <-  rbinom(1, (round((N/2)*(1 - YI.t))), qY.wt.sel)/(round((N/2)*(1 - YI.t)))
			qY.del.t    <-  rbinom(1, (round((N/2)*(1 - YI.t))), qY.del.sel)/(round((N/2)*(1 - YI.t)))
			qI.wt.t     <-  rbinom(1, (round((N/2)*YI.t)), qI.wt.sel)/(round((N/2)*YI.t))
			XaSp.wt.t   <-  rbinom(1, (N/2), XaSp.wt.sel)/(N/2)
			XaSp.del.t  <-  rbinom(1, (N/2), XaSp.del.sel)/(N/2)

		}
		if(YI.t == 1){
		fix  <-  fix + 1
		}
	}
	PrFix[j]  <-  fix/sims
	cat('\r', paste("Progress:", round(100*(j/length(r))), "% complete"))
}

PrR  <-  dpois(r,lambda=n*qHat)
dev.new()
par(mfrow=c(2,3))
plot(PrFix ~ r, main="Pr(Fix | r)")
plot(PrR ~ r, main="Prob. capturing r del. alleles")
plot((PrFix*PrR) ~ r, main="Pr(Fix | r) * Pr(r)")
plot((PrFix*sims) ~ r, main="# fixations")
plot((PrR*sims)  ~ r, main="# inv. per r bin")
plot((PrR*sims*PrFix) ~ r, , main="# fixations")

nTot*x*qHat*0.8
nTot*PrFix[r == (nTot*x*qHat*0.8)]




########################################
## Wright-Fisher using Exact recursions derived in Mathematica using Xf, Xm, Y
## inluding sampling variance for XaOv, XaSp, etc.
# Gradient across inversion size (x)

# Parameters
h    <-  0.1
s    <-  0.01
U    <-  10*s
nTot  <-  10000
u     <-  U/nTot
qHat  <-  (U/(nTot*h*s))
N     <-  1000
YI.0   <-  2/N

# Containers
PrFix <-  c()
YI.t        <- c()
XaOv.wt.t   <- c()
XaSp.wt.t   <- c()
qI.wt.t     <- c()
qY.wt.t     <- c()
XaOv.del.t  <- c()
XaSp.del.t  <- c()
qY.del.t    <- c()

#number of simulations 
sims  = 100*N/2

# inversion sizes
invSize  <-  c(0.5,1:9)/10

# Loop over inversion size
for(j in 1:length(invSize)) {
fix   = 0
	
	# Loop over replicate simulations
	for(i in 1:sims){

		# initial frequencies
		YI.t        <- YI.0
		XaOv.wt.t   <- qHat
		XaSp.wt.t   <- qHat
		qI.wt.t     <- 0
		qY.wt.t     <- qHat
		XaOv.del.t  <- qHat
		XaSp.del.t  <- qHat
		qY.del.t    <- qHat

		# Draw random value for # loci captured by inversion
		n     <-  rpois(1,lambda=nTot*invSize[j])

		# Draw random value for # del. mutations captured by inversion (r) given x
		r   <-  rbinom(1, size=n, prob=qHat)

		# Run forward simulation
		while(YI.t*(1 - YI.t) > 0) {

			#expected frequency after selection
			XaOv.wt.sel   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t), digits=9)
			XaOv.del.sel  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t), digits=9)
			XaSp.wt.sel   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
			XaSp.del.sel  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
			qY.wt.sel     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
			qY.del.sel    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
			qI.wt.sel     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
			YI.sel        <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, qt.Y.del=qY.del.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)

			#binomial sampling
			XaOv.wt.t   <-  rbinom(1, N, XaOv.wt.sel)/N
			XaOv.del.t  <-  rbinom(1, N, XaOv.del.sel)/N
			YI.t        <-  rbinom(1, (N/2), YI.sel)/(N/2)
			qY.wt.t     <-  rbinom(1, (round((N/2)*(1 - YI.t))), qY.wt.sel)/(round((N/2)*(1 - YI.t)))
			qY.del.t    <-  rbinom(1, (round((N/2)*(1 - YI.t))), qY.del.sel)/(round((N/2)*(1 - YI.t)))
			qI.wt.t     <-  rbinom(1, (round((N/2)*YI.t)), qI.wt.sel)/(round((N/2)*YI.t))
			XaSp.wt.t   <-  rbinom(1, (N/2), XaSp.wt.sel)/(N/2)
			XaSp.del.t  <-  rbinom(1, (N/2), XaSp.del.sel)/(N/2)

		}
		if(YI.t == 1){
		fix  <-  fix + 1
		}
	}
	PrFix[j]  <-  fix/sims
	cat('\r', paste("Progress:", round(100*(j/length(invSize))), "% complete"))
}

par(mfrow=c(1,1))
plot(PrFix ~ invSize, xlim=c(0,max(invSize)), ylim=c(0,max(PrFix)), main="Pr(Fix | x)")
abline(h=2/N, lwd=2)
text(x=0.9, y=(2/N)+0.00075, labels=c("2/N"))







########################################
## Multilocus Wright-Fisher Simulations
## Using exact recursions derived in Mathematica using Xf, Xm, Y
## inluding sampling variance for XaOv, XaSp, etc.
## 
# Gradient across inversion size (x)

wBarY.multi  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del) {
	pt.I.wt    <-  1 - qt.I.wt
	pt.Y.del   <-  1 - qt.Y.del
	pt.Y.wt    <-  1 - qt.Y.wt
		  YI.t*(prod( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del) )) * prod(1 - s*(h*(pt.I.wt*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt)) +
	(1 - YI.t)*( prod(1 - s*(h*((1 - XaOv.t.del)*qt.Y.del + XaOv.t.del*pt.Y.del) + XaOv.t.del*qt.Y.del)) * prod((1 - s*(h*((1 - XaOv.t.wt)*qt.Y.wt + XaOv.t.wt*pt.Y.wt) + XaOv.t.wt*qt.Y.wt))) )
}

YI.multi.prime  <-  function(n, r, s, h, YI.t, qt.I.wt, qt.Y.wt, qt.Y.del, XaOv.t.wt, XaOv.t.del, ...) {
	(YI.t*(prod( 1 - s*(h*(1 - XaOv.t.del) + XaOv.t.del))) * prod(1 - s*(h*((1 - qt.I.wt)*XaOv.t.wt + qt.I.wt*(1 - XaOv.t.wt)) + qt.I.wt*XaOv.t.wt))) / wBarY.multi(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qt.I.wt, qt.Y.wt=qt.Y.wt, qt.Y.del=qt.Y.del, XaOv.t.wt=XaOv.t.wt, XaOv.t.del=XaOv.t.del)
}






# Gradient across inversion size (x)

# Parameters
h    <-  0.1
s    <-  0.01
U    <-  5*s
nTot  <-  1000
u     <-  U/nTot
qHat  <-  (U/(nTot*h*s))
N     <-  1000
YI.0   <-  2/N

# Containers
PrFix <-  c()
YI.t        <- c()
XaOv.wt.t   <- c()
XaSp.wt.t   <- c()
qI.wt.t     <- c()
qY.wt.t     <- c()
XaOv.del.t  <- c()
XaSp.del.t  <- c()
qY.del.t    <- c()

#number of simulations 
sims  = 100*N/2

# inversion sizes
invSize  <-  c(0.5,1:9)/10

# Loop over inversion size
for(j in 1:length(invSize)) {
fix   = 0
	
	# Loop over replicate simulations
	for(i in 1:sims){

		# initial frequencies

		# Draw random value for # loci captured by inversion
		n   <-  rpois(1,lambda=nTot*invSize[j])

		# Draw random initial frequencies for del. alleles at each locus
		qi  <-  rbinom(n, size = 2*N, prob=qHat)/(2*N)

		# Draw random value for # del. mutations captured by inversion (r) given x
		ri  <-  rbinom(n, size=1, prob=qi)
		r   <-  sum(ri)

		# Assign 
		YI.t        <- YI.0
		XaOv.wt.t   <- qi[ri == 0]
		XaSp.wt.t   <- qi[ri == 0]
		qI.wt.t     <- 0
		qY.wt.t     <- qi[ri == 0]
		XaOv.del.t  <- qi[ri == 1]
		XaSp.del.t  <- qi[ri == 1]
		qY.del.t    <- qi[ri == 1]


		# Run forward simulation
		while(YI.t*(1 - YI.t) > 0) {

			#expected frequency after selection
			XaOv.wt.sel   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t, XaSp.t.wt=XaSp.wt.t), digits=9)
			XaOv.del.sel  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t, XaSp.t.del=XaSp.del.t), digits=9)
			XaSp.wt.sel   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
			XaSp.del.sel  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
			qY.wt.sel     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t, XaOv.t.wt=XaOv.wt.t), digits=9)
			qY.del.sel    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t, XaOv.t.del=XaOv.del.t), digits=9)
			qI.wt.sel     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)
			YI.sel        <-  round(YI.multi.prime(n=n, r=r, s=s, h=h, YI.t=YI.t, qt.I.wt=qI.wt.t, qt.Y.wt=qY.wt.t, qt.Y.del=qY.del.t, XaOv.t.wt=XaOv.wt.t, XaOv.t.del=XaOv.del.t), digits=9)

			#binomial sampling
			XaOv.wt.t   <-  rbinom(n, N, prob=XaOv.wt.sel)/N
			XaOv.del.t  <-  rbinom(r, N, XaOv.del.sel)/N
			YI.t        <-  rbinom(1, (N/2), YI.sel)/(N/2)
			qY.wt.t     <-  rbinom(n, (round((N/2)*(1 - YI.t))), qY.wt.sel)/(round((N/2)*(1 - YI.t)))
			qY.del.t    <-  rbinom(r, (round((N/2)*(1 - YI.t))), qY.del.sel)/(round((N/2)*(1 - YI.t)))
			qI.wt.t     <-  rbinom(n, (round((N/2)*YI.t)), qI.wt.sel)/(round((N/2)*YI.t))
			XaSp.wt.t   <-  rbinom(n, (N/2), XaSp.wt.sel)/(N/2)
			XaSp.del.t  <-  rbinom(r, (N/2), XaSp.del.sel)/(N/2)

		}
		if(YI.t == 1){
		fix  <-  fix + 1
		}
	}
	PrFix[j]  <-  fix/sims
	cat('\r', paste("Progress:", round(100*(j/length(invSize))), "% complete"))
}

par(mfrow=c(1,1))
plot(PrFix ~ invSize, xlim=c(0,max(invSize)), ylim=c(0,max(PrFix)), main="Pr(Fix | x)")
abline(h=2/N, lwd=2)
text(x=0.9, y=(2/N)+0.00075, labels=c("2/N"))


d  <-  data.frame("invSize" = invSize, "PrFix" = PrFix)
write.csv(d, file="./output/data/simResults/neutralPartialRecessive_Pfix_Ud5s_N1k.csv", row.names=FALSE)

getwd()







###################################################
###################################################
## Older Code trying to follow Nei et al. (1967)


########################################
## Wright-Fisher using Exact Recursions
rm(list=ls())
# Parameters
h    <-  0.1
s    <-  0.01
U    <-  10*s
nTot  <-  10000
NrR   <-  rpois(n*x,lambda=U*x/(h*s))
#r     <-  c(range(NrR)[1]:range(NrR)[2])
r     <-  c(0:range(NrR)[2])
#r     <-  c(0:5)
u     <-  U/nTot
qHat  <-  (U/(nTot*h*s))
x     <-  0.2
n     <-  nTot*x
N     <-  1000
YI0   <-  2/N
PrFix <-  c()


YI.t      <- c()
qN.del.t  <- c()
qN.wt.t   <- c()
qI.wt.t   <- c()


#number of simulations 
sims  = 100*N/2

# Loop over # deleterious mutations initially captured by inversion
for(j in 1:length(r)) {
fix   = 0
	
	# Loop over replicate simulations
	for(i in 1:sims){
		YI      <-  YI0
		qN.del  <-  qHat
		qN.wt   <-  qHat
		qI.wt   <-  0

		# Run forward simulation
		while(YI*(1 - YI) > 0) {

			#expected frequency after selection
			F.I.sel     <-  YI.next(n=n, r=r[j], s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt)
			if(r[j] == 0) {
				qN.del.sel  <-  round(qN.del.next(n=n, r=r[j], u=u, s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt), digits=9)
			}
			if(r[j] > 0) {
				qN.del.sel  <-  round(qN.del.next(n=n, r=r[j], u=u, s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt), digits=9)
			}
			qN.wt.sel   <-  round(qN.wt.next(n=n, r=r[j], u=u, s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt), digits=9)
			qI.wt.sel   <-  round(qI.wt.next(n=n, r=r[j], u=u, s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt), digits=9)

			#binomial sampling
			YI      <-  rbinom(1, (N/2), F.I.sel)/(N/2)
#			qN.del  <-  rbinom(1, (2*N - (YI*(N/2))), qN.del.sel)/(2*N - (YI*(N/2)))
#			qN.wt   <-  rbinom(1, (2*N - (YI*(N/2))), qN.wt.sel)/(2*N - (YI*(N/2)))
#			qI.wt   <-  rbinom(1, (YI*(N/2)), qI.wt.sel)/(N/2)
		}
		if(YI == 1){
		fix  <-  fix + 1
		}
	}
	PrFix[j]  <-  fix/sims
	cat('\r', paste("Progress:", round(100*(j/length(r))), "% complete"))
}

PrR  <-  dpois(r,lambda=U*x/(h*s))
par(mfrow=c(2,2))
plot(PrFix ~ r)
plot(PrR ~ r)
plot((PrFix*PrR) ~ r)
plot(PrR*10000)

plot((PrR*10000*PrFix) ~ r)


sum(PrR*10000)
sum(PrR*10000*PrFix)

nTot*x*qHat*0.8
10000*PrFix[r == (nTot*x*qHat*0.8)]









#######################################
## Wright-Fisher w/ multi-locus Exact Recursions

# Parameters
h     <-  0.1
s     <-  0.01
U     <-  10*s
nTot  <-  10000
NrR   <-  rpois(n*x,lambda=U*x/(h*s))
#r     <-  c(range(NrR)[1]:range(NrR)[2])
r     <-  c(0:range(NrR)[2])
#r     <-  c(0:5)
u     <-  U/nTot
qHat  <-  (U/(nTot*h*s))
x     <-  0.2
n     <-  nTot*x
N     <-  1000
YI0   <-  2/N
PrFix <-  c()

	QN.del  <-   (2*N - (YI*(N/2)))*qN.del/(2*N - (YI*(N/2))) # qN.del/(1 - YI)
	QN.wt   <-  (2*N - (YI*(N/2)))*qN.wt / (2*N - (YI*(N/2))) #qN.wt/(1 - YI)
	QI.wt   <-  qI.wt/YI
		  YI*(prod( 1 - s*(h*(1 - QN.del) + QN.del) ) * prod(1 - s*(h*((1 - QI.wt)*QN.wt + QI.wt*(1 - QN.wt)) + QI.wt*QN.wt))) +
	(1 - YI)*( prod(1 - s*(2*h*QN.del*(1 - QN.del) + QN.del^2)) * prod(1 - s*(2*h*QN.wt*(1 - QN.wt) + QN.wt^2)))



wBarI.multi  <-  function(n, r, s, h, YIt, qNt.del, qNt.wt, qIt.wt) {
	QNt.del  <-  qNt.del/(1 - YIt)
	QNt.wt   <-  qNt.wt/(1 - YIt)
	QIt.wt   <-  qIt.wt/YIt
		  YIt*(prod( 1 - s*(h*(1 - QNt.del) + QNt.del) ) * prod(1 - s*(h*((1 - QIt.wt)*QNt.wt + QIt.wt*(1 - QNt.wt)) + QIt.wt*QNt.wt))) +
	(1 - YIt)*( prod(1 - s*(2*h*QNt.del*(1 - QNt.del) + QNt.del^2)) * prod(1 - s*(2*h*QNt.wt*(1 - QNt.wt) + QNt.wt^2)))
}

YI.multi.next  <-  function(n, r, s, h, YIt, qNt.del, qNt.wt, qIt.wt, ...) {
	QNt.del  <-  qNt.del/(1 - YIt)
	QNt.wt   <-  qNt.wt/(1 - YIt)
	QIt.wt   <-  qIt.wt/YIt
	(YIt*( prod( 1 - s*(h*(1 - QNt.del) + QNt.del) )) * prod(1 - s*(h*((1 - QIt.wt)*QNt.wt + QIt.wt*(1 - QNt.wt)) + QIt.wt*QNt.wt))) / wBarI.multi(n=n, r=r, s=s, h=h, YIt=YIt, qNt.del=qNt.del, qNt.wt=qNt.wt, qIt.wt=qIt.wt)
}

qI.wt.multi.next  <-  function(n, r, u, s, h, YIt, qNt.del, qNt.wt, qIt.wt) {
	pNt.wt  <-  1 - YIt - qNt.wt
	pIt.wt  <-  YIt - qIt.wt
	( (qIt.wt*(1 - s*(h*pNt.wt + qNt.wt)) + u*pIt.wt*wBar.i.wt(s=s, h=h, qNt.wt=qNt.wt, qIt.wt=qIt.wt, YIt=YIt)) / ( YIt - qIt.wt*s*(h*pNt.wt + qNt.wt) ) ) * YI.multi.next(n=n, r=r, s=s, h=h, YIt=YIt, qNt.del=qNt.del, qIt.wt=qIt.wt, qNt.wt=qNt.wt)
}

qN.wt.multi.next  <-  function(n, r, u, s, h, YIt, qNt.del, qNt.wt, qIt.wt) {
	pNt.wt  <-  1 - YIt - qNt.wt
	pIt.wt  <-  YIt - qIt.wt
	(3/4)*(1 - YI.multi.next(n=n, r=r, s=s, h=h, YIt=YIt, qNt.del=qNt.del, qIt.wt=qIt.wt, qNt.wt=qNt.wt))*
		( qNt.wt*(1 - s*(h*(pNt.wt + YIt*(pIt.wt - pNt.wt)/2) + qNt.wt + YIt*(qIt.wt - qNt.wt)/2)) + u*pNt.wt*wBar.i.wt(s=s, h=h, qNt.wt=qNt.wt, qIt.wt=qIt.wt, YIt=YIt) ) / 
			( 1 - (YIt/4) - qNt.wt*s*(h*(pNt.wt + YIt*(pIt.wt - pNt.wt)/2) + qNt.wt + YIt*(qIt.wt - qNt.wt)/2) )
}

qN.del.multi.next  <-  function(n, r, u, s, h, YIt, qNt.del, qNt.wt, qIt.wt) {
	pNt.del  <-  1 - YIt - qNt.del
	(3/4)*(1 - YI.multi.next(n=n, r=r, s=s, h=h, YIt=YIt, qNt.del=qNt.del, qIt.wt=qIt.wt, qNt.wt=qNt.wt))*
		(qNt.del*(1 - s*(h*(pNt.del - YIt*pNt.del/2) + qNt.del + YIt*pNt.del/2)) + (u*pNt.del*(1 - YIt/2))*wBar.i.del(s=s, h=h, qNt.del=qNt.del, YIt=YIt)) / 
			(1 - (YIt/4) - qNt.del*s*(h*(pNt.del - YIt*pNt.del/2) + qNt.del + YIt*pNt.del/2))
}

#number of simulations 
sims  = 100*N/2

# Loop over # deleterious mutations initially captured by inversion
for(j in 1:length(r)) {
fix   = 0
	
	# Loop over replicate simulations
	for(i in 1:sims){
		YI      <-  YI0
		qN.del  <-  rbinom(r[j], (2*N - (YI*(N/2))), qHat)/(2*N - (YI*(N/2)))
		qN.wt   <-  rbinom((n - r[j]), (2*N - (YI*(N/2))), qHat)/(2*N - (YI*(N/2)))
		qI.wt   <-  rep(0, times = (n-r[j]))

		# Run forward simulation
#		plot(YI ~ 0, xlim=c(0,5000), ylim=c(0,1))
#		t  <-  1
		while(YI*(1 - YI) > 0) {

			#expected frequency after selection
			if(YI.multi.next(n=n, r=r[j], s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt) > 1) {
				browser()
			}
			F.I.sel     <-  YI.multi.next(n=n, r=r[j], s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt)
			qN.del.sel  <-  qN.del.multi.next(n=n, r=r[j], u=u, s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt)
			qN.wt.sel   <-  qN.wt.multi.next(n=n, r=r[j], u=u, s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt)
			qI.wt.sel   <-  qI.wt.multi.next(n=n, r=r[j], u=u, s=s, h=h, YIt=YI, qNt.del=qN.del, qNt.wt=qN.wt, qIt.wt=qI.wt)

			#binomial sampling
			YI      <-  rbinom(1, (N/2), F.I.sel)/(N/2)
			qN.del  <-  apply(as.matrix(qN.del.sel), 1, function(x) rbinom(1, (2*N - (YI*(N/2))), x)/(2*N - (YI*(N/2))))
			qN.wt   <-  apply(as.matrix(qN.wt.sel), 1, function(x) rbinom(1, (2*N - (YI*(N/2))), x)/(2*N - (YI*(N/2))))
			qI.wt   <-  apply(as.matrix(qI.wt.sel), 1, function(x) rbinom(1, (YI*(N/2)), x)/(N/2))

#points(YI ~ t, cex=0.5)
#t  <-  t + 1
		}
		if(YI == 1){
		fix  <-  fix + 1
		}
	cat('\r', paste("Progress:", round(100*(i/sims)), "% complete"))
	}
	PrFix[j]  <-  fix/sims
	cat('\r', paste("Progress:", round(100*(j/length(r))), "% complete"))
}

PrR  <-  dpois(r,lambda=U*x/(h*s))
par(mfrow=c(2,2))
plot(PrFix ~ r)
plot(PrR ~ r)
plot((PrFix*PrR) ~ r)
plot(PrR*10000)

plot((PrR*10000*PrFix) ~ r)


sum(PrR*10000)
sum(PrR*10000*PrFix)

nTot*x*qHat*0.8
10000*PrFix[r == (nTot*x*qHat*0.8)]
