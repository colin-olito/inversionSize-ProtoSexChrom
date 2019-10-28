 par.list  <-  list(
					  N    =  1000,		# Population size
					  gen  =  1000,		# Number of generations for determ. sim. of eq. frequencies
					  reps =  100,		# Number of replicates
					  sm   =  0.01,		# selection coefficient in males
					  sf   =  0.01,		# selection coefficient in females
					  hm   =  1/2,		# dominance coefficient in males
					  hf   =  1/2,		# dominance coefficient in females
					  r    =  1/2,		# recombination rate
					  A    =  4,		# frequency of SA loci on the chromosome (Poisson rate parameter=)
					  P    =  0.1,		# Size of sl-PAR (fraction of chromosome)
					  Ud   =  0.1,		# Deleterious mutation rate
					  sd   =  0.05		# Selection coefficient for deleterious mutations.
					)















ss  <-  c()
a.aPAR  <-  c()
a.slPAR  <-  c()
for(i in 1:10000) {
			capturedAlleles.aPAR   <-  sum(rbinom(n = (n - n.slPAR), prob=qHat.aPAR, size=1))
			capturedAlleles.slPAR  <-  sum(rbinom(n = n.slPAR, prob=qHat.slPAR, size=1))
a.aPAR[i]  <-  capturedAlleles.aPAR
a.slPAR[i]  <- capturedAlleles.slPAR 

			ss[i]  <-  ((n-n.slPAR)*sm*(1-qHat.aPAR)*(1 - hm*(1-2*qHat.aPAR)-qHat.aPAR) + 
						n.slPAR*sm*(1-qHat.slPAR)*(1 - hm*(1-2*qHat.slPAR)-qHat.slPAR)) -
					  ((n-n.slPAR-capturedAlleles.aPAR)*sm*(1 - hm*(1-2*qHat.aPAR)-qHat.aPAR) + 
						(n.slPAR-capturedAlleles.slPAR)*sm*(1 - hm*(1-2*qHat.slPAR)-qHat.slPAR))

}
mean(ss)
mean(a.aPAR)
mean(a.slPAR)

hist(ss,breaks=30)

par(mfrow=c(2,2))
N  <-  10^4
 par.list  <-  list(
					  N    =  N,		# Population size
					  gen  =  10^4,		# Number of generations for determ. sim. of eq. frequencies
					  reps =  10^5,		# Number of replicates
					  sm   =  0.01,		# selection coefficient in males
					  sf   =  0.01,		# selection coefficient in females
					  hm   =  1/2,		# dominance coefficient in males
					  hf   =  1/2,		# dominance coefficient in females
					  r    =  1/2,		# recombination rate
					  A    =  2,		# frequency of SA loci on the chromosome (Poisson rate parameter=)
					  P    =  0.2,		# Size of sl-PAR (fraction of chromosome)
					  Ud   =  0.1,		# Deleterious mutation rate
					  sd   =  0.05		# Selection coefficient for deleterious mutations.
					)
res1  <-  simPrInvMultiLocusSA(par.list)

par.list$P  <-  0.001
res2  <-  simPrInvMultiLocusSA(par.list)

par.list$r  <-  0.005
par.list$P  <-  0.2
res3  <-  simPrInvMultiLocusSA(par.list)

par.list$P  <-  0.001
res4  <-  simPrInvMultiLocusSA(par.list)

ymax  <-  max(res1$simPrFix,res2$simPrFix,res3$simPrFix,res4$simPrFix)
plot(simPrFix   ~ invSize, ylim=c(0,ymax*1.05), xlim=c(0,1), data = res1)
points(simPrFix ~ invSize, col=2, data = res2)
points(simPrFix ~ invSize, pch=21, bg=3, data = res3)
points(simPrFix ~ invSize, pch=21, bg=4, data = res4)

ymax  <-  max(res1$simPrFixDel,res2$simPrFixDel,res3$simPrFixDel,res4$simPrFixDel)
plot(simPrFixDel   ~ invSize, ylim=c(0,ymax*1.05), xlim=c(0,1), data = res1)
points(simPrFixDel ~ invSize, col=2, data = res2)
points(simPrFixDel ~ invSize, pch=21, bg=3, data = res3)
points(simPrFixDel ~ invSize, pch=21, bg=4, data = res4)




mean(res3$invSize*(res3$simPrFixDel/sum(res3$simPrFixDel)))

weighted.mean(res3$invSize,res3$simPrFixDel)
weighted.mean(res3$invSize,(res3$simPrFixDel/sum(res3$simPrFixDel)))
	weighted.mean(res4$invSize,res4$simPrFixDel)

lines(approxPrFix ~ invSize, lty=2, lwd=2, col=1, data = res1)
lines(approxPrFix ~ invSize, lty=2, lwd=2, col=2, data = res2)
lines(approxPrFix ~ invSize, lty=2, lwd=2, col=3, data = res3)
lines(approxPrFix ~ invSize, lty=2, lwd=2, col=3, data = res4)

plot(approxPrFixDel ~ invSize, type='l', lty=2, lwd=2, col=1, data = res1)
lines(approxPrFixDel ~ invSize, lty=2, lwd=2, col=2, data = res2)
lines(approxPrFixDel ~ invSize, lty=2, lwd=2, col=3, data = res3)
lines(approxPrFixDel ~ invSize, lty=2, lwd=2, col=4, data = res4)
abline(h=0)


res4
plot(simPrFix ~ invSize, ylim=c(0,0.15), xlim=c(0,1), data = res1)
points(simPrFix ~ invSize, ylim=c(0,0.15), xlim=c(0,1), pch=21, bg=1, data = res3)
points(simPrFix ~ invSize, ylim=c(0,0.15), xlim=c(0,1), col=2, data = res2)
points(simPrFix ~ invSize, ylim=c(0,0.15), xlim=c(0,1), pch=21, bg=2, data = res4)







N  <-  30^4
 par.list  <-  list(
					  N    =  N,		# Population size
					  gen  =  10^4,		# Number of generations for determ. sim. of eq. frequencies
					  reps =  10^5,		# Number of replicates
					  locSA = c(0,2),	# location of SA loci c(slPAR, aPAR)
					  sm   =  0.01,		# selection coefficient in males
					  sf   =  0.01,		# selection coefficient in females
					  hm   =  1/2,		# dominance coefficient in males
					  hf   =  1/2,		# dominance coefficient in females
					  A    =  2,		# frequency of SA loci on the chromosome (Poisson rate parameter=)
					  r    =  1/2,		# recombination rate
					  Ud   =  0.1,		# Deleterious mutation rate
					  sd   =  0.05		# Selection coefficient for deleterious mutations.
					)
res1  <-  simPrInvMultiLocusSATwoLocus(par.list)
par.list$locSA  <-  c(2,0)
res2  <-  simPrInvMultiLocusSATwoLocus(par.list)
par.list$locSA  <-  c(1,1)
res3  <-  simPrInvMultiLocusSATwoLocus(par.list)
par.list$r  <-  0.005
par.list$locSA  <-  c(0,2)
res4  <-  simPrInvMultiLocusSATwoLocus(par.list)
par.list$locSA  <-  c(2,0)
res5  <-  simPrInvMultiLocusSATwoLocus(par.list)
par.list$locSA  <-  c(1,1)
res6  <-  simPrInvMultiLocusSATwoLocus(par.list)


par(mfrow=c(2,2))
ymax  <-  max(res1$simPrFix,res2$simPrFix,res3$simPrFix,res4$simPrFix,res5$simPrFix,res6$simPrFix)
plot(simPrFix   ~ invSize, ylim=c(0,ymax*1.05), xlim=c(0,1), type='l', lwd=2, lty=2, data = res1)
lines(simPrFix ~ invSize, col=2, lwd=2, lty=2, data = res2)
lines(simPrFix ~ invSize, col=4, lwd=2, lty=2, data = res3)
lines(simPrFix ~ invSize, col='#252525', lwd=2, lty=2, data = res4)
lines(simPrFix ~ invSize, col='tomato', lwd=2, lty=2, data = res5)
lines(simPrFix ~ invSize, col='dodgerblue', lwd=2, lty=2, data = res6)

ymax  <-  max(res1$simPrFixDel,res2$simPrFixDel,res3$simPrFixDel,res4$simPrFixDel,res5$simPrFixDel,res6$simPrFixDel)
plot(simPrFixDel   ~ invSize, ylim=c(0,ymax*1.05), xlim=c(0,1), type='l', lwd=2, lty=2, data = res1)
lines(simPrFixDel ~ invSize, col=2, lwd=2, lty=2, data = res2)
lines(simPrFixDel ~ invSize, col=4, lwd=2, lty=2, data = res3)
lines(simPrFixDel ~ invSize, col='#252525', lwd=2, lty=2, data = res4)
lines(simPrFixDel ~ invSize, col='tomato', lwd=2, lty=2, data = res5)
lines(simPrFixDel ~ invSize, col='dodgerblue', lwd=2, lty=2, data = res6)

ymax  <-  max(res1$simPrFixCatch,res2$simPrFixCatch,res3$simPrFixCatch,res4$simPrFixCatch,res5$simPrFixCatch,res6$simPrFixCatch)
plot(simPrFixCatch   ~ invSize, ylim=c(0,ymax*1.05), xlim=c(0,1), type='l', lwd=2, lty=2, data = res1)
lines(simPrFixCatch ~ invSize, col=2, lwd=2, lty=2, data = res2)
lines(simPrFixCatch ~ invSize, col=4, lwd=2, lty=2, data = res3)
lines(simPrFixCatch ~ invSize, col='#252525', lwd=2, lty=2, data = res4)
lines(simPrFixCatch ~ invSize, col='tomato', lwd=2, lty=2, data = res5)
lines(simPrFixCatch ~ invSize, col='dodgerblue', lwd=2, lty=2, data = res6)

ymax  <-  max(res1$simPrFixDelCatch,res2$simPrFixDelCatch,res3$simPrFixDelCatch,res4$simPrFixDelCatch,res5$simPrFixDelCatch,res6$simPrFixDelCatch)
plot(simPrFixDelCatch   ~ invSize, ylim=c(0,ymax*1.05), xlim=c(0,1), type='l', lwd=2, lty=2, data = res1)
lines(simPrFixDelCatch ~ invSize, col=2, lwd=2, lty=2, data = res2)
lines(simPrFixDelCatch ~ invSize, col=4, lwd=2, lty=2, data = res3)
lines(simPrFixDelCatch ~ invSize, col='#252525', lwd=2, lty=2, data = res4)
lines(simPrFixDelCatch ~ invSize, col='tomato', lwd=2, lty=2, data = res5)
lines(simPrFixDelCatch ~ invSize, col='dodgerblue', lwd=2, lty=2, data = res6)
