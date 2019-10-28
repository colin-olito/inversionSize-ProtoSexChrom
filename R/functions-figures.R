###############
# DEPENDENCIES
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)
library(plyr)

#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('output/figures', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}

##############################################################
##  Analytic Approximations

pFixNeutralY  <-  function(x, Ny, Ud, sd) {
    1/(Ny*exp(Ud*x/sd))
}

pFixBeneficialY  <-  function(x, sI, Ud, sd) {
    2*sI*(1 + ( (Ud*x)/(1 - (1 - sI)*exp(-sd))))*exp((-2*Ud*x)/sd)
}
logisticEq  <-  function(x, a, b) {
    exp(a+b*x)/(1 + exp(a+b*x))
}
qHat_SAAdd  <-  function(sm, sf) {
    (sm - sf + sm*sf) / (2*sm*sf)
}
sI_SAAdd  <-  function(sm, qHat) {
    (sm*qHat*(1 - qHat))/2
}
pFix_SAUnlinkedAdd  <-  function(x, sm, qHat, A, Ud, sd) {
    sm*qHat*(1 - qHat)*(x^2)*A*exp(-x*A)*( (Ud*x)/(1 - (1 - (sm*qHat*(1 - qHat))/2)*exp(-sd)) )*exp((-2*Ud*x)/sd)
}
pFix_SALinkedAdd  <-  function(x, sI, A, P, Ud, sd) {
    2*sI*x*A*exp(-A*P)*( (Ud*x)/(1 - (1 - sI)*exp(-sd)) )*exp((-2*Ud*x)/sd)
}

meanpFixSAUnlinked  <-  function(A, Ud, sd) {
    (3*sd)/(A*sd + 2*Ud)
}
sciNotation <- function(x, digits = 1) {
    if (length(x) > 1) {
        return(append(sciNotation(x[1]), sciNotation(x[-1])))
    }
    if (!x) return(0)
    exponent <- floor(log10(x))
    base <- round(x / 10^exponent, digits)
    as.expression(substitute(base %*% 10^exponent, 
            list(base = base, exponent = exponent)))
}
##############################################################
##############################################################
##  Final figures for paper

fixationProbabilityFigure  <-  function() {

    ## import data.frames
    # Neutral sims
    neutralNy10     <-  read.csv(file = "./output/data/simResults/neutral-Y-pFix_Ny10_Ud0.1_sd0.05.csv")
    neutralNy100    <-  read.csv(file = "./output/data/simResults/neutral-Y-pFix_Ny100_Ud0.1_sd0.05.csv")
    neutralNy1000   <-  read.csv(file = "./output/data/simResults/neutral-Y-pFix_Ny1000_Ud0.1_sd0.05.csv")

    # beneficial sims
    benNy10     <-  read.csv(file = "./output/data/simResults/beneficial-Y-pFix_Ny10_sI0.02_Ud0.1_sd0.05.csv")
    benNy50     <-  read.csv(file = "./output/data/simResults/beneficial-Y-pFix_Ny50_sI0.02_Ud0.1_sd0.05.csv")
    benNy100    <-  read.csv(file = "./output/data/simResults/beneficial-Y-pFix_Ny100_sI0.02_Ud0.1_sd0.05.csv")
    benNy1000   <-  read.csv(file = "./output/data/simResults/beneficial-Y-pFix_Ny1000_sI0.02_Ud0.1_sd0.05.csv")

    # sexually antagonistic sims
    SAaddLinked     <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI_sf0.02_sm0.02_hf0.5_hm0.5.csv")

    # set constants
    x     <-  seq(0,1,by=0.005)
    Ny    <-  1*10^3
    Ud    <-  0.1
    sd    <-  0.05
    sb    <-  0.02
    sf    <-  0.02
    sm    <-  0.02
    r_SA  <-  SAaddLinked$r[2]
    A     <-  1
    P     <-  0.05

    # set plot layout
    layout.mat <- matrix(c(1:3), nrow=1, ncol=3, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A: Neutral Inversions
    pFixNeutral  <-  pFixNeutralY(x=x, Ny=Ny, Ud=Ud, sd=sd)
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,5,5,2), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(pFixNeutral*Ny ~ x, lty=1, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points(PrFixSimDel*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.1), data=neutralNy10)
        points(PrFixSimDel*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.3), data=neutralNy100)
        points(PrFixSimDel*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.7), data=neutralNy1000)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(0, 1/4, 1/2, 3/4, 1), labels = NA)
        labelPos  <-  axTicks(2)/0.91 + 0.035
        proportionalLabel(-0.13,  labelPos[5],   expression(paste("1/", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.13,  labelPos[4],   expression(paste("3/4", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.13,  labelPos[3],   expression(paste("1/2", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.13,  labelPos[2],   expression(paste("1/4", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.13,  labelPos[1],   expression(paste("0.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        # Plot labels etc.
        proportionalLabel( 0.025,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.2,   expression(paste("Neutral")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Pr(fix | ", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(x))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
#               x       =  usr[2]*0.955,
               x       =  usr[2]*0.995,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Eq(1a)"))),
               lty     =  1,
               lwd     =  2,
               col      =  c(transparentColor('#252525', opacity=1.0)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.93,
               legend  =  c(
                            expression(paste(italic(N[Y]), " = ",10^3)),
                            expression(paste(italic(N[Y]), " = ",10^2)),
                            expression(paste(italic(N[Y]), " = ",10))),
               pch     =  c(21,21,21),
               pt.bg      =  c(transparentColor('#252525', opacity=0.7),
                               transparentColor('#252525', opacity=0.3),
                               transparentColor('#252525', opacity=0.1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



## Panel B: Beneficial Inverions
    pFixBeneficial  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.1*max(pFixBeneficial/sb)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(pFixBeneficial/sb ~ x, lty=1, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points(PrFixSimDel/sb ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.1), data=benNy50)
        points(PrFixSimDel/sb ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.3), data=benNy100)
        points(PrFixSimDel/sb ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.7), data=benNy1000)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        labelPos  <-  axTicks(2)/2.4 + 0.035
#        proportionalLabel(-0.13,  labelPos[6],   expression(paste("5/4", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.12,  labelPos[5],   expression(paste("2",italic(s[I]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.12,  labelPos[4],   expression(paste("3", italic(s[I]),"/4")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.12,  labelPos[3],   expression(paste(italic(s[I]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.12,  labelPos[2],   expression(paste(italic(s[I]),"/2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.12,  labelPos[1],   expression(paste("0.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        # Annnotations
        proportionalLabel( 0.025,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.2,   expression(paste("Beneficial")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.2,  0.5,   expression(paste("Pr(fix | ", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(x))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Plot annotations
#        proportionalLabel( 0.5,   0.25,  expression(paste("SA poly.")), cex=0.8, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2]*0.995,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Eq(4)"))),
               lty     =  1,
               lwd     =  2,
               col      =  c(transparentColor('#252525', opacity=1.0)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel C: Sexual Antagonism (2-locus models)
    qHatSA             <-  qHat_SAAdd(sf=sf, sm=sm)
    sISAAdd            <-  sI_SAAdd(sm=sm, qHat=qHatSA)
    pFixSAunlinkedAdd  <-  pFix_SAUnlinkedAdd(x=x, sm=sm, qHat=qHatSA, A=A, Ud=Ud, sd=sd)
    meanUnlinked       <-  meanpFixSAUnlinked(A=A, Ud=Ud, sd=sd)
    sISAaddLinked      <-  SAaddLinked$sI[SAaddLinked$r == r_SA]
    pFixSAlinked       <-  pFix_SALinkedAdd(x=x, sI=sISAaddLinked, A=A, P=P, Ud=Ud, sd=sd)
    RelpFixSAlinked    <-  ((pFixSAlinked/sISAaddLinked)/max(pFixSAlinked/sISAaddLinked))*max(pFixSAunlinkedAdd/sISAAdd)
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.2*max(c((pFixSAunlinkedAdd),(pFixSAlinked)))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(pFixSAlinked ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1))
        # Left axes
        axis(2, las=1, at=axTicks(2), labels=sciNotation(axTicks(2),1))
        par(new=TRUE)
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.2*max(c((pFixSAunlinkedAdd)))), ylab='', xlab='', cex.lab=1.2)
        lines(pFixSAunlinkedAdd ~ x, lty=1, lwd=2, col=transparentColor('#252525', opacity=1))
        abline(v=c((sd/Ud),meanUnlinked), lwd=2, lty=c(3,3), col='black')
        # Right Axes
        axis(4, las=1, at=axTicks(2), labels=sciNotation(axTicks(2),1))
        # x Axis
        axis(1, las=1)
        # Plot labels etc.
        proportionalLabel( 0.025,  1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.2,   expression(paste("Sex Antag. (2-locus)")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 1.4,  0.5,   expression(paste("Pr(fix | ", italic(x), ", ",italic(r), " = 1/2)")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(x))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Annotations
        proportionalLabel( 0.85,  0.3,  expression(paste(italic(s[d]/U[d]))), cex=0.8, adj=c(0.5, 0.5), xpd=NA)        
        arrows(x0 = usr[1] + 0.73*(usr[2] - usr[1]), 
               y0 = usr[3] + 0.165*(usr[4] - usr[3]), 
               x1 = usr[1] + 0.51*(usr[2] - usr[1]), 
               y1 = usr[3] + 0.165*(usr[4] - usr[3]), angle=30, length=0.05)
        proportionalLabel( 0.85,  0.15,  expression(paste(italic(3*s[d]))), cex=0.8, adj=c(0.5, 0.5), xpd=NA)
        segments(x0 = usr[1] + 0.75*(usr[2] - usr[1]), 
                 y0 = usr[3] + 0.07*(usr[4] - usr[3]), 
                 x1 = usr[1] + 0.95*(usr[2] - usr[1]), 
                 y1 = usr[3] + 0.07*(usr[4] - usr[3]))
        proportionalLabel( 0.85,  0.05,  expression(paste(italic(A*s[d]+2*U[d]))), cex=0.8, adj=c(0.5, 0.5), xpd=NA)        
        arrows(x0 = usr[1] + 0.73*(usr[2] - usr[1]), 
               y0 = usr[3] + 0.07*(usr[4] - usr[3]), 
               x1 = usr[1] + 0.6*(usr[2] - usr[1]), 
               y1 = usr[3] + 0.07*(usr[4] - usr[3]), angle=30, length=0.05)
        # Legend
        legend(
               x       =  usr[2]*0.45,
               y       =  usr[4]*0.5,
               legend  =  c(
                            expression(paste("r = 0.002")),
                            expression(paste("r = 1/2"))),
               lty     =  c(2,1),
               lwd     =  2,
               col      =  c(transparentColor('#252525', opacity=1.0)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}





##############################################################
##############################################################
##  Exploratory figures



#' Creates figure to validate approximations for invasion probability for 
#' a Y-linked inversion capturing the sex-determining locus and a male- 
#' beneficial allele at a second SA locus 
#'
#' @title Plots for SA Y-linked inversion validation of invasion probability
#' @param s.df     filename for selection gradient simulations. 
#'                 Data available in ./output/data/simResults/SA/ 
#' @param r.df     filename for recombination gradient simulations. 
#'                 Data available in ./output/data/simResults/SA/ 
#' @export
SA_Y_validatePinv  <-  function(s.df, r.df) {

    # importa data.frames
    ext  <-  "./output/data/simResults/SA/"
    sDat  <-  read.csv(paste(ext,s.df, sep=""), header=TRUE)
    rDat  <-  read.csv(paste(ext,r.df, sep=""), header=TRUE)

    # Set plot layout
    layout.mat <- matrix(c(1:2), nrow=1, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    SA_polyLims  <-  c(sDat$sm[1]/(1 + sDat$sm[1]),sDat$sm[1]/(1 - sDat$sm[1]))

## Panel A: sel.grad 
    par(omi=c(0.5, 0.5, 0.5, 0.5), mar = c(4,4,0.5,0.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(min(sDat$sf),max(sDat$sf)), ylim = c(0,1.05*max(sDat$pInv, sDat$pInvApprox2)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        polygon(x=c(SA_polyLims, rev(SA_polyLims)), 
            y=c(usr[3],usr[3],usr[4],usr[4]), 
            col=transparentColor('grey80', opacity=0.4), border=NA)
        lines(sDat$pInvApprox2 ~ sDat$sf, lty=1, lwd=2, col=transparentColor('#252525', opacity=0.75))
        points(sDat$pInv ~ sDat$sf, pch=21, bg=transparentColor('dodgerblue', opacity=0.5), col=1)
        abline(v=sDat$sm[1]/(1 + sDat$sm[1]), col=1, lty=1)
        abline(v=sDat$sm[1]/(1 - sDat$sm[1]), col=1, lty=1)
        abline(h=1/sDat$N[1], lwd=2, lty=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.45,  1.12,   expression(paste(italic(N)," = ")), cex=1.0, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.6 ,  1.12,  substitute(N,list(N=sDat$N[1])), cex=1.0, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.35,  1.05,   expression(paste(italic(s[m])," =         ,")), cex=1.0, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.4,   1.055,  substitute(sm,list(sm=sDat$sm[1])), cex=1.0, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.55,  1.05,   expression(paste(italic(r)," =")), cex=1.0, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.65, 1.055,  substitute(r,list(r=sDat$r[1])), cex=1.0, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Pr(", italic(inv.), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Plot annotations
        proportionalLabel( 0.5,   0.25,  expression(paste("SA poly.")), cex=0.8, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.85,   0.075,  expression(paste(1/N)), cex=0.8, adj=c(0.5, 0.5), xpd=NA)

## Panel B: r.grad 
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,max(rDat$r)), ylim = c(0,1.05*max(rDat$pInv, sDat$pInvApprox2)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        # pInv
        lines(rDat$pInvApprox2 ~ rDat$r, lty=1, lwd=2, col=transparentColor('#252525', opacity=0.75))
        points(rDat$pInv ~ rDat$r, pch=21, bg=transparentColor('dodgerblue', opacity=0.6), col=1)
        abline(h=1/rDat$N[1], lwd=2, lty=2, col='black')
        # axes
        axis(1, las=1)
        axis(2, las=1, labels = NA)
        # Plot labels etc.
        proportionalLabel( 0.05,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.45,  1.12,   expression(paste(italic(N)," = ")), cex=1.0, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.6 ,  1.12,  substitute(N,list(N=sDat$N[1])), cex=1.0, adj=c(0.5, 0.5), xpd=NA)

        proportionalLabel( 0.45,  1.05,   expression(paste(italic(s[m])," = ", italic(s[f])," = ")), cex=1.0, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.62,  1.06,  substitute(sm,list(sm=sDat$sm[1])), cex=1.0, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.85,   0.075,  expression(paste(1/N)), cex=0.8, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(r))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Sim.")),
                            expression(paste(2(lambda - 1)))),
               pch     =  c(21,NA),
               pt.bg   =  c(transparentColor('dodgerblue', opacity=0.5),1),
               col     =  c(1,transparentColor('#252525', opacity=0.75)),
               lty     =  c(NA,1),
               lwd     =  c(NA,2),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}
