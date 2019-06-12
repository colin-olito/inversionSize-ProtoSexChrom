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
##############################################################
##  Final figures for paper







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
