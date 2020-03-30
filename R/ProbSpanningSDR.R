# Simulate a bunch of inversions with 
# -- a gradient of locations for the SDR and inversion 
# -- a gradient of inversion sizes 
reps  <-  10000
SDRLoc  <- c(1:20)/20
xs      <- c(1:11)/11
pCatch  <-  c()
for(i in 1:length(SDRLoc)) {	
	for(j in 1:length(xs)) {
		LHS  <-  runif(n=reps, min=0, max=(1 - xs[j]))
		pCatch[((i-1)*length(xs)) + j]   <-  sum(LHS < SDRLoc[i] & SDRLoc[i] < (LHS + xs[j])) / reps
	}
}

# Make a data frame for plotting
res  <-  data.frame("SDRLoc"   =  rep(SDRLoc, each=length(xs)),
					"xs"       = rep(xs,times=length(SDRLoc)),
					"pCatch"   =  pCatch)

# What are y1 & y2 (the proportions of the chromosome arm to the LEFT & RIGHT of the SDR)
y1  <-  res$SDRLoc
y2  <-  1 - y1


# Analytic expressions for Pr(span SDR | x)
PrCatch  <-  rep(NA,times=nrow(res))
PrCatch[res$xs < y1 & res$xs < y2]  <- ((res$xs)/(1-res$xs))[res$xs < y1 & res$xs < y2]
PrCatch[res$xs > y1 & res$xs < y2]  <-  (y1/(1-res$xs))[res$xs > y1 & res$xs < y2]
PrCatch[res$xs < y1 & res$xs > y2]  <-  (y2/(1-res$xs))[res$xs < y1 & res$xs > y2]
PrCatch[res$xs > y1 & res$xs > y2]  <-  1

res$PrCatch  <-  PrCatch
head(res)

par(mfrow=c(1,2))

# Pr(catch SDR) as a function of the location of the SDR
# different lines are plotted for each inversion size (x)
plot(pCatch ~ SDRLoc, ylim=c(0,1), xlim=c(0,1), pch=NA, data=res)
for(i in 1:length(xs)) {
		lines(pCatch[xs==xs[i]] ~ SDRLoc[xs==xs[i]], lwd=2, col=i, data=res)
		lines(PrCatch[xs==xs[i]] ~ SDRLoc[xs==xs[i]], lwd=2, lty=2, col=1, data=res)
}

# Pr(catch SDR) as a function of inversion size, now with
# different lines plotted for each SDR location
plot(pCatch ~ xs, ylim=c(0,1), xlim=c(0,1), pch=NA, data=res)
SDRs  <-  SDRLoc
for(i in 1:(length(SDRLoc)/2)) {
		lines(pCatch[SDRLoc==SDRs[i]] ~ xs[SDRLoc==SDRs[i]], lwd=2, col=i, data=res)
		lines(PrCatch[SDRLoc==SDRs[i]] ~ xs[SDRLoc==SDRs[i]], lwd=2, lty=2, col=1, data=res)
#i=i+1
}
