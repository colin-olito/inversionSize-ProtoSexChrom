###############
# DEPENDENCIES
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)
library(plyr)
library(lattice)
library(latticeExtra)
library(wesanderson)

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


proportionalArrows <- function(px1, py1, px2, py2, adj=c(0, 1), log=FALSE, length=length, ...) {
    usr  <-  par('usr')
    x.p1  <-  usr[1] + px1*(usr[2] - usr[1])
    y.p1  <-  usr[3] + py1*(usr[4] - usr[3])
    x.p2  <-  usr[1] + px2*(usr[2] - usr[1])
    y.p2  <-  usr[3] + py2*(usr[4] - usr[3])
    if(log=='x') {
        x.p1  <-  10^(x.p1)
        x.p2  <-  10^(x.p2)
    }
    if(log=='y') {
        y.p1  <-  10^(y.p1)
        y.p2  <-  10^(y.p2)
    }
    if(log=='xy') {
        x.p1  <-  10^(x.p1)
        y.p1  <-  10^(y.p1)
        x.p2  <-  10^(x.p2)
        y.p2  <-  10^(y.p2)
    }
    arrows(x0=x.p1, y0=y.p1, x1=x.p2, y1=y.p2, length=length,...)
}


proportionalRect <- function(px1, py1, px2, py2, adj=c(0, 1), log=FALSE, border=TRUE, ...) {
    usr  <-  par('usr')
    x.p1  <-  usr[1] + px1*(usr[2] - usr[1])
    y.p1  <-  usr[3] + py1*(usr[4] - usr[3])
    x.p2  <-  usr[1] + px2*(usr[2] - usr[1])
    y.p2  <-  usr[3] + py2*(usr[4] - usr[3])
    if(log=='x') {
        x.p1  <-  10^(x.p1)
        x.p2  <-  10^(x.p2)
    }
    if(log=='y') {
        y.p1  <-  10^(y.p1)
        y.p2  <-  10^(y.p2)
    }
    if(log=='xy') {
        x.p1  <-  10^(x.p1)
        y.p1  <-  10^(y.p1)
        x.p2  <-  10^(x.p2)
        y.p2  <-  10^(y.p2)
    }
    rect(xleft=x.p1, ybottom=y.p1, xright=x.p2, ytop=y.p2,...)
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

fibonacci.scale  <-  function(n) {
    fibs  <-  c(0,1)
    for(i in 2:n) {
        fibs  <-  c(fibs, (fibs[i] + fibs[i-1]))
    }
    (fibs/max(fibs))[-2]
}

filled.contour3  <-  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            wesCol = 'Zissou1', col = wes_palette(wesCol, length(levels) - 1, 'continuous'), 
            plot.title, plot.axes, key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = FALSE, frame.plot = axes, mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page

  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
 # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
 # on.exit(par(par.orig))
 # w <- (3 + mar.orig[2]) * par("csi") * 2.54
 # par(las = las)
 # mar <- mar.orig
                            #plot.new()
 # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
       col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

filled.legend <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                            length.out = ncol(z)), z, xlim=range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes, ...) 
   {
    # modification of filled.contour by Carey McGilliard and Bridget Ferris
    # designed to just plot the legend
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
      if (axes) 
        axis(4)
    }
    else key.axes
    box()
  }
##############################################################
##  Analytic Approximations

pCatchDel  <-  function(Ud, x, sd) {
    exp(-Ud*x/sd)
}
PrCatchSLR  <-  function(SLR_pos, x) {
    y1   <-  SLR_pos
    y2   <-  1 - y1
    p1   <-  x/(1 - x)
    p2   <-  y1/(1 - x)
    p3   <-  y2/(1 - x)
    p4   <-  rep(1, length(x))
    res  <-  apply(rbind(p1,p2,p3,p4),2,min)
    res
}
pFixNeutralY  <-  function(Ny) {
    Ny^(-1)
}

pFixBeneficialY  <-  function(x, sI, Ud, sd) {
    2*sI*(1 + ((Ud*x)/(1 - (1 - sI)*exp(-sd))))*exp((-Ud*x)/sd)
}
logisticEq  <-  function(x, a, b) {
    exp(a+b*x)/(1 + exp(a+b*x))
}
qHat_SAAdd  <-  function(sm, sf) {
    (sm - sf + sm*sf) / (2*sm*sf)
}
sI_SAAdd  <-  function(sm, qHat) {
    (sm*(1 - qHat))/2
}
sI_SAAdd2  <-  function(sm, qHat) {
    ((sm - qHat*sm)/(2 - 2*(1 - qHat)*sm))
}
#pFix_SAUnlinkedAdd  <-  function(x, sm, qHat, A, Ud, sd) {
#    2*((sm - qHat*sm)/(2 - 2*(1 - qHat)*sm))*(x^2)*A*exp(-x*A)*(1 + (Ud*x)/(1 - (1 - ((sm - qHat*sm)/(2 - 2*(1 - qHat)*sm))*exp(-sd))) )*exp((-Ud*x)/sd)
#}
pFix_SAUnlinkedAdd  <-  function(x, sm, qHat, A, Ud, sd) {
    2*sI_SAAdd(sm=sm, qHat=qHat)*x*A*exp(-x*A)*(1 + ( (Ud*x)/(1 - (1 - sI_SAAdd(sm=sm, qHat=qHat))*exp(-sd))))*exp((-Ud*x)/sd)
}
pFix_SALinkedAdd  <-  function(x, sI, A, P, Ud, sd) {
    2*sI*A*P*exp(-A*P)*(1 + ( (Ud*x)/(1 - (1 - sI)*exp(-sd))))*exp((-Ud*x)/sd)
}
sciNotation <- function(x, digits = 1) {
    if (length(x) > 1) {
        return(append(sciNotation(x[1]), sciNotation(x[-1])))
    }
    if (!x) return(0)
    exponent <- floor(log10(x))
    base     <- round(x / 10^exponent, digits)
    as.expression(substitute(base %*% 10^exponent, 
            list(base = base, exponent = exponent)))
}
distNewInvRBP  <-  function(x) {
    2*(1 - x)
}
distNewInvEXP  <-  function(lambda, x) {
    (lambda*exp(-lambda*x)) / (1 - exp(-lambda))
}

##############################################################
##############################################################
##  Figures - partially recessive del. mut.

pFixFig_Neutral_Ben  <-  function() {

    ## import data.frames
    # Neutral sims
    neutral_N10k  <-  read.csv(file = "./output/data/simResults/shelter_PrFixFig_h0.25_s0.01_N10k_deterministic_q.csv")

    # beneficial sims
    ben_N10k <-  read.csv(file = "./output/data/simResults/beneficial_PrFixFig_sI0.02_h0.25_s0.01_N10k_deterministic_q.csv")


    # set constants
    x     <-  seq(0,1,by=0.005)
    Ny    <-  10^4/2

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)

    # set plot layout
    layout.mat <- matrix(c(1:2), nrow=2, ncol=1, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A: Neutral Inversions
    pFixNeutral  <-  pFixNeutralY(Ny=Ny)
    par(omi=c(0.25, 0.5, 0.25, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^6), 2/(10^2)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # pInv
        abline(h=pFixNeutral, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[1], opacity=0.4), data=neutral_N10k)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[1], opacity=0.4), data=neutral_N10k)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=24, bg=transparentColor(COLS[1], opacity=0.4), data=neutral_N10k)
        # Axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), 
                labels=c(0, expression(2/10^5), expression(2/N), expression(2/10^3), expression(2/10^2), expression(2/10^1)))
        # Plot labels etc.
        proportionalLabel(0.5, 1.06, expression("Neutral") , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')        
        # Legend
        legend(
               x       =  1.02,
               y       =  0.025,
               legend  =  c(
                            expression(paste(italic(U/hs), " = ", 8)),
                            expression(paste(italic(U/hs), " = ", 20)),
                            expression(paste(italic(U/hs), " = ", 40))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg     =  transparentColor(COLS[2], opacity=0.6),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel B: Beneficial Inverions
    sb  <-  0.02
    sd  <-  0.01
    pFixBeneficial  <-  pFixBeneficialY(x=x, sI=sb, Ud=0.02, sd=sd)
    pFixBeneficial2  <-  pFixBeneficialY(x=x, sI=sb, Ud=0.05, sd=sd)
    pFixBeneficial3  <-  pFixBeneficialY(x=x, sI=sb, Ud=0.1, sd=sd)
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,3), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        abline(h=2, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        points(PrFix[U == 0.02]/sI[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[1], opacity=0.4), data=ben_N10k)
        points(PrFix[U == 0.05]/sI[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[1], opacity=0.4), data=ben_N10k)
        points(PrFix[U == 0.1]/sI[U == 0.1] ~ x[U == 0.1], pch=24, bg=transparentColor(COLS[1], opacity=0.4), data=ben_N10k)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(0, 1, 2, 3), 
                labels=c(0, 
                         expression(italic(s[b])), 
                         expression(2*italic(s[b])), 
                         expression(3*italic(s[b]))))
        # Plot labels etc.
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(B)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste("Beneficial")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}



## Fixation probability ~ r 
## equal SA selection, varying U/hs
pFixFig_SA_r  <-  function() {

    ## import data.frames
    SA_add_r1  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.002_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r2  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.008_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r3  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.032_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r4  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.124_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r5  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.498_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")

    SA_DomRev_r1  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.002_sf0.05_sm0.05_hSA0.25_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_DomRev_r2  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.008_sf0.05_sm0.05_hSA0.25_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_DomRev_r3  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.032_sf0.05_sm0.05_hSA0.25_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_DomRev_r4  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.124_sf0.05_sm0.05_hSA0.25_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_DomRev_r5  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.498_sf0.05_sm0.05_hSA0.25_h0.25_s0.01_NumEqFreq_N10k.csv")

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey90'))
    COLS     <-  colfunc(5)

    # set plot layout
    layout.mat <- matrix(c(1:6), nrow=2, ncol=3, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A -- U/hs = 8 -- Additive SA Fitness Effects
    par(omi=c(0.2, 0.5, 0.2, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[1], opacity=0.75), data = SA_add_r1)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[2], opacity=0.75), data = SA_add_r2)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[3], opacity=0.75), data = SA_add_r3)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.75), data = SA_add_r4)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[5], opacity=0.75), data = SA_add_r5)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(A)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste(italic(U)/italic(hs)==8)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.5,  0.5,   expression(paste(italic(h[SA]==1/2))), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        

## Panel B -- U/hs = 20 -- Additive SA Fitness Effects
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[1], opacity=0.75), data = SA_add_r1)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[2], opacity=0.75), data = SA_add_r2)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[3], opacity=0.75), data = SA_add_r3)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[4], opacity=0.75), data = SA_add_r4)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[5], opacity=0.75), data = SA_add_r5)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Annnotations
        proportionalLabel( 0.5,   1.25,   expression(paste("Sexually Antagonistic Selection")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03,  1.075, expression(paste(B)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste(italic(U)/italic(hs)==20)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

## Panel C -- U/hs = 40 -- Additive SA Fitness Effects
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[1], opacity=0.75), data = SA_add_r1)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[2], opacity=0.75), data = SA_add_r2)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[3], opacity=0.75), data = SA_add_r3)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[4], opacity=0.75), data = SA_add_r4)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[5], opacity=0.75), data = SA_add_r5)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(C)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste(italic(U)/italic(hs)==40)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  1.02,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(r), " = ", 0.5)),
                            expression(paste(italic(r), " = ", 0.124)),
                            expression(paste(italic(r), " = ", 0.032)),
                            expression(paste(italic(r), " = ", 0.008)),
                            expression(paste(italic(r), " = ", 0.002))),
               pch     =  21,
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(transparentColor(COLS[5], opacity=0.75),
                            transparentColor(COLS[4], opacity=0.75),
                            transparentColor(COLS[3], opacity=0.75),
                            transparentColor(COLS[2], opacity=0.75),
                            transparentColor(COLS[1], opacity=0.75)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



## Panel D -- U/hs = 8 -- DomRev SA Fitness Effects
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.06), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[1], opacity=0.75), data = SA_DomRev_r1)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[2], opacity=0.75), data = SA_DomRev_r2)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[3], opacity=0.75), data = SA_DomRev_r3)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.75), data = SA_DomRev_r4)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[5], opacity=0.75), data = SA_DomRev_r5)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 0.03, 1.075, expression(paste(D)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.5,  0.5,   expression(paste(italic(h[SA]==1/4))), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5, -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

## Panel E -- U/hs = 20 -- DomRev SA Fitness Effects
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.06), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[1], opacity=0.75), data = SA_DomRev_r1)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[2], opacity=0.75), data = SA_DomRev_r2)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[3], opacity=0.75), data = SA_DomRev_r3)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[4], opacity=0.75), data = SA_DomRev_r4)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[5], opacity=0.75), data = SA_DomRev_r5)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Annnotations
        proportionalLabel( 0.03, 1.075, expression(paste(E)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5, -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

## Panel F -- U/hs = 40 -- DomRev SA Fitness Effects
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.06), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[1], opacity=0.75), data = SA_DomRev_r1)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[2], opacity=0.75), data = SA_DomRev_r2)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[3], opacity=0.75), data = SA_DomRev_r3)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[4], opacity=0.75), data = SA_DomRev_r4)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[5], opacity=0.75), data = SA_DomRev_r5)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(F)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}






## Fixation probability ~ r 
## biased SA selection (sf = sm/(1 - sm)), varying U/hs
pFixFig_SA_biasedSel_r  <-  function() {

    ## import data.frames
    SA_add_biased_r1  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.002_sf0.0526315789473684_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_biased_r2  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.008_sf0.0526315789473684_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_biased_r3  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.032_sf0.0526315789473684_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_biased_r4  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.124_sf0.0526315789473684_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_biased_r5  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.498_sf0.0526315789473684_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey90'))
    COLS     <-  colfunc(5)

    # set plot layout
    layout.mat <- matrix(c(1:6), nrow=1, ncol=3, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A -- U/hs = 8 -- Additive SA Fitness Effects
    par(omi=c(0.2, 0.5, 0.2, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[1], opacity=0.75), data = SA_add_biased_r1)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[2], opacity=0.75), data = SA_add_biased_r2)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[3], opacity=0.75), data = SA_add_biased_r3)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.75), data = SA_add_biased_r4)
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[5], opacity=0.75), data = SA_add_biased_r5)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 0.03, 1.075, expression(paste(A)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste(italic(U)/italic(hs)==8)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

## Panel B -- U/hs = 20 -- DomRev SA Fitness Effects
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[1], opacity=0.75), data = SA_add_biased_r1)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[2], opacity=0.75), data = SA_add_biased_r2)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[3], opacity=0.75), data = SA_add_biased_r3)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[4], opacity=0.75), data = SA_add_biased_r4)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=21, bg=transparentColor(COLS[5], opacity=0.75), data = SA_add_biased_r5)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Annnotations
        proportionalLabel( 0.5,   1.25,   expression(paste("Biased SA Selection: ", italic(s[f])==italic(s[m])/(1-italic(s[m])))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03, 1.075, expression(paste(B)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste(italic(U)/italic(hs)==20)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

## Panel c -- U/hs = 40 -- DomRev SA Fitness Effects
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[1], opacity=0.75), data = SA_add_biased_r1)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[2], opacity=0.75), data = SA_add_biased_r2)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[3], opacity=0.75), data = SA_add_biased_r3)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[4], opacity=0.75), data = SA_add_biased_r4)
        points(PrFix[U == 0.1] ~ x[U == 0.1], pch=21, bg=transparentColor(COLS[5], opacity=0.75), data = SA_add_biased_r5)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(C)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste(italic(U)/italic(hs)==40)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  1.02,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(r), " = ", 0.5)),
                            expression(paste(italic(r), " = ", 0.124)),
                            expression(paste(italic(r), " = ", 0.032)),
                            expression(paste(italic(r), " = ", 0.008)),
                            expression(paste(italic(r), " = ", 0.002))),
               pch     =  21,
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(transparentColor(COLS[5], opacity=0.75),
                            transparentColor(COLS[4], opacity=0.75),
                            transparentColor(COLS[3], opacity=0.75),
                            transparentColor(COLS[2], opacity=0.75),
                            transparentColor(COLS[1], opacity=0.75)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}









sI_SA  <-  function(sm, hm, Yhat, Xfhat) {
    (sm*(1 - Yhat)*(1 - Xfhat - hm*(1 - 2*Xfhat))) / (1 - sm*(1 - Xfhat - Yhat*(1 - hm - Xfhat) + hm*Xfhat*(1 - 2*Yhat)))
}

#' Combined Pr(fix | x) figure
#'  including neutral, beneficial, and SA scenarios
## Fixation probability ~ r 
## equal SA selection, varying U/hs
pFixFig_combined  <-  function() {

    ## import data.frames
    # Neutral sims
    neutral_N10k  <-  read.csv(file = "./output/data/simResults/shelter_PrFixFig_h0.25_s0.01_N10k_deterministic_q.csv")

    # beneficial sims
    ben_N10k <-  read.csv(file = "./output/data/simResults/beneficial_PrFixFig_sI0.02_h0.25_s0.01_N10k_deterministic_q.csv")

    ## SA selection
    SA_add_r1  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.002_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r2  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.008_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r3  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.032_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r4  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.124_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r5  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.498_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")

    # set constants
    x     <-  seq(0,1,by=0.005)
    Ny    <-  10^4/2
    A     <-  1
    P     <-  0.05

    # Calculate deterministic s_I for SA selection
    log(0.002)
    log(0.498)
    logSeq  <-  seq(from=-6.214608, to=-0.6971552, len=5)
    rSeq  <-  round(exp(logSeq), digits=3)
    rSeq[4]  <-  rSeq[4] - 0.001
    eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-EqualSel_sf0.05_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
        # subset sI data to these values
        eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]
        sI.vec  <-  c()
        for(i in 1:nrow(eqFreq.subdat)) {
            sI.vec[i]  <-  sI_SA(sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], Yhat = eqFreq.subdat$Y[i], Xfhat = eqFreq.subdat$Xf[i])
        }
        sI_SA_vals  <-  data.frame(cbind(eqFreq.subdat[,1:5], sI.vec))

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey90'))
    COLS     <-  colfunc(5)

    # set plot layout
    layout.mat <- matrix(c(1:6), nrow=2, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)



## Panel A: Neutral Inversions
    pFixNeutral  <-  pFixNeutralY(Ny=Ny)
    par(omi=c(0.25, 0.5, 0.5, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^6), 2/(10^2)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # pInv
        abline(h=pFixNeutral, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        points(PrFix[U == 0.1]  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        # Axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), 
                labels=c(0, expression(2/10^5), expression(2/N), expression(2/10^3), expression(2/10^2), expression(2/10^1)))
        # Plot labels etc.
        proportionalLabel(0.5, 1.25, expression("Neutral") , cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')        
        # Legend
        legend(
               x       =  1.02,
               y       =  0.025,
               legend  =  c(
                            expression(paste(italic(U/hs), " = ", 8)),
                            expression(paste(italic(U/hs), " = ", 20)),
                            expression(paste(italic(U/hs), " = ", 40))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg     =  transparentColor(COLS[4], opacity=0.8),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel B: Beneficial Inversions
    sb  <-  0.02
    sd  <-  0.01
    pFixBeneficial  <-  pFixBeneficialY(x=x, sI=sb, Ud=0.02, sd=sd)
    pFixBeneficial2  <-  pFixBeneficialY(x=x, sI=sb, Ud=0.05, sd=sd)
    pFixBeneficial3  <-  pFixBeneficialY(x=x, sI=sb, Ud=0.1, sd=sd)
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,3), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        abline(h=2, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        points(PrFix[U == 0.02]/sI[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        points(PrFix[U == 0.05]/sI[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        points(PrFix[U == 0.1]/sI[U == 0.1]   ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(0, 1, 2, 3), 
                labels=c(0, 
                         expression(italic(s[b])), 
                         expression(2*italic(s[b])), 
                         expression(3*italic(s[b]))))
        # Plot labels etc.
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(B)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.25,   expression(paste("Beneficial")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## Panel C -- Additive SA Fitness Effects -- r = 0.002
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        abline(h=2*sI_SA_vals$sI.vec[1], lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        points(PrFix[U == 0.1]  ~ x[U == 0.1], pch=24,  bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 1.15,   1.25,   expression(paste("Sexually Antagonistic Selection")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03,  1.075, expression(paste(C)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.002)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  1.02,
               y       =  0.07,
               legend  =  c(
                            expression(paste(italic(2*s[I])))),
               lty     =  2,
               lwd     =  2,
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



## Panel E -- U/hs = 40 -- Additive SA Fitness Effects
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        abline(h=2*sI_SA_vals$sI.vec[3], lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r3)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r3)
        points(PrFix[U == 0.1]  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r3)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(E)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.032)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
         # Legend
        legend(
               x       =  1.02,
               y       =  0.07,
               legend  =  c(
                            expression(paste(italic(2*s[I])))),
               lty     =  2,
               lwd     =  2,
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



## Panel D -- Additive SA Fitness Effects -- r = 0.002
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        abline(h=2*sI_SA_vals$sI.vec[2], lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r2)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r2)
        points(PrFix[U == 0.1]  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r2)

        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(D)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.008)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  1.02,
               y       =  0.07,
               legend  =  c(
                            expression(paste(italic(2*s[I])))),
               lty     =  2,
               lwd     =  2,
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )




## Panel F -- Additive SA Fitness Effects -- r_SA =0.5
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.07), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        abline(h=2*sI_SA_vals$sI.vec[5], lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        points(PrFix[U == 0.1]  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Annnotations
        proportionalLabel( 0.03, 1.075, expression(paste(F)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.5)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5, -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  1.02,
               y       =  0.07,
               legend  =  c(
                            expression(paste(italic(2*s[I])))),
               lty     =  2,
               lwd     =  2,
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


}









#' Combined Pr(fix | x) figure
#'  including neutral, beneficial, and SA scenarios
## Fixation probability ~ r 
## equal SA selection, varying U/hs
pFix_pCatchSAFig_combined  <-  function() {

    ## import data.frames
    # Neutral sims
    neutral_N10k  <-  read.csv(file = "./output/data/simResults/shelter_PrFixFig_h0.25_s0.01_N10k_deterministic_q.csv")

    # beneficial sims
    ben_N10k <-  read.csv(file = "./output/data/simResults/beneficial_PrFixFig_sI0.02_h0.25_s0.01_N10k_deterministic_q.csv")

    ## SA selection
    SA_add_r1  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.002_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r5  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.498_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")

    # set constants
    x     <-  seq(0,1,by=0.005)
    Ny    <-  10^4/2
    A     <-  1
    P     <-  0.05

    # Calculate deterministic s_I for SA selection
    log(0.002)
    log(0.498)
    logSeq  <-  seq(from=-6.214608, to=-0.6971552, len=5)
    rSeq  <-  round(exp(logSeq), digits=3)
    rSeq[4]  <-  rSeq[4] - 0.001
    eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-EqualSel_sf0.05_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
        # subset sI data to these values
        eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]
        sI.vec  <-  c()
        for(i in 1:nrow(eqFreq.subdat)) {
            sI.vec[i]  <-  sI_SA(sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], Yhat = eqFreq.subdat$Y[i], Xfhat = eqFreq.subdat$Xf[i])
        }
        sI_SA_vals  <-  data.frame(cbind(eqFreq.subdat, sI.vec))

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey90'))
    COLS     <-  colfunc(5)

    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)



## Panel A: Neutral Inversions
    pFixNeutral  <-  pFixNeutralY(Ny=Ny)
    par(omi=c(0.25, 0.25, 0.5, 0.25), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^6), 2/(10^2)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # pInv
        abline(h=pFixNeutral, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        points(PrFix[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        points(PrFix[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        points(PrFix[U == 0.1]  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        # Axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), 
                labels=c(0, expression(2/10^5), expression(2/N), expression(2/10^3), expression(2/10^2), expression(2/10^1)))
        # Plot labels etc.
        proportionalLabel(0.5, 1.2, expression("Neutral") , cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')        


## Panel B: Beneficial Inversions
    sb  <-  0.02
    sd  <-  0.01
    pFixBeneficial   <-  pFixBeneficialY(x=x, sI=sb, Ud=0.02, sd=sd)
    pFixBeneficial2  <-  pFixBeneficialY(x=x, sI=sb, Ud=0.05, sd=sd)
    pFixBeneficial3  <-  pFixBeneficialY(x=x, sI=sb, Ud=0.1, sd=sd)
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,3), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        abline(h=2, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        points(PrFix[U == 0.02]/sI[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        points(PrFix[U == 0.05]/sI[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        points(PrFix[U == 0.1]/sI[U == 0.1]   ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(0, 1, 2, 3), 
                labels=c(0, 
                         expression(italic(s[b])), 
                         expression(2*italic(s[b])), 
                         expression(3*italic(s[b]))))
        # Plot labels etc.
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(B)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.2,   expression(paste("Beneficial")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        


## Panel C -- Additive SA Fitness Effects -- r = 1/2
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.01), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        lines(2*sI_SA_vals$sI.vec[5]*sI_SA_vals$Y[5]*A*x*exp(-A*x) ~ x, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
#        abline(h=2*sI_SA_vals$sI.vec[5], lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02]*A*x[U == 0.02]*exp(-A*x[U == 0.02])*sI_SA_vals$Y[5] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        points(PrFix[U == 0.05]*A*x[U == 0.05]*exp(-A*x[U == 0.05])*sI_SA_vals$Y[5] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        points(PrFix[U == 0.1]*A*x[U == 0.1]*exp(-A*x[U == 0.1])*sI_SA_vals$Y[5]  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 0.5,   1.2,   expression(paste("Sexually Antagonistic Selection")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03, 1.075, expression(paste(C)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.5)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)



## Panel D -- Additive SA Fitness Effects -- r = 0.002
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        abline(h=2*sI_SA_vals$sI.vec[1]*sI_SA_vals$Y[1]*A*P*exp(-A*P), lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02]*sI_SA_vals$Y[1]*A*P*exp(-A*P) ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        points(PrFix[U == 0.05]*sI_SA_vals$Y[1]*A*P*exp(-A*P) ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        points(PrFix[U == 0.1]*sI_SA_vals$Y[1]*A*P*exp(-A*P)  ~ x[U == 0.1], pch=24,  bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(D)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.002)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5, -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*1,
               y       =  usr[4]*0.98,
               legend  =  c(
                            expression(paste(Pr(fix), " ignoring del. mut."))),
               lty     =  2,
               lwd     =  2,
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )        # Legend
        legend(
               x       =  usr[2]*1,
               y       =  usr[4]*0.85,
               legend  =  c(
                            expression(paste(italic(U/hs), " = ", 8)),
                            expression(paste(italic(U/hs), " = ", 20)),
                            expression(paste(italic(U/hs), " = ", 40))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg     =  transparentColor(COLS[4], opacity=0.8),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



}



#' Probability of spanning SLR & fixing
#' 

PrSpanSLRFixFigure  <-  function(SLR_pos = c(1/2, 1/10)) {

    ## import data.frames
    # Neutral sims
    neutral_N10k  <-  read.csv(file = "./output/data/simResults/shelter_PrFixFig_h0.25_s0.01_N10k_deterministic_q.csv")

    # beneficial sims
    ben_N10k <-  read.csv(file = "./output/data/simResults/beneficial_PrFixFig_sI0.02_h0.25_s0.01_N10k_deterministic_q.csv")

    ## SA selection
    SA_add_r1  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.002_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r5  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.498_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey90'))
    COLS     <-  colfunc(5)

    # set constants
    x   <-  seq(0,1,by=0.005)
    Ny  <-  10^4/2
    A   <-  1
    P   <-  0.05
    SLR_positions  <-  c(0.5, 0.1)

    # Calculate deterministic s_I for SA selection
    log(0.002)
    log(0.498)
    logSeq  <-  seq(from=-6.214608, to=-0.6971552, len=5)
    rSeq  <-  round(exp(logSeq), digits=3)
    rSeq[4]  <-  rSeq[4] - 0.001
    eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-EqualSel_sf0.05_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
        # subset sI data to these values
        eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]
        sI.vec  <-  c()
        for(i in 1:nrow(eqFreq.subdat)) {
            sI.vec[i]  <-  sI_SA(sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], Yhat = eqFreq.subdat$Y[i], Xfhat = eqFreq.subdat$Xf[i])
        }
        sI_SA_vals  <-  data.frame(cbind(eqFreq.subdat, sI.vec))


    # set plot layout
    layout.mat <- matrix(c(1:8), nrow=2, ncol=4, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

### ROW 1: SLR_pos = 1/2
    SLR_pos  <-  SLR_positions[1]
## Panel A: Neutral Inversions
    pFixNeutral  <-  pFixNeutralY(Ny=Ny)
    par(omi=c(0.25, 0.5, 0.5, 0.25), mar = c(4.5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/10^6), 4/(10^3)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # pInv
        lines(pFixNeutral*PrCatchSLR(SLR_pos=SLR_pos, x = x) ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFix[U == 0.02]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.02]) ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        points(PrFix[U == 0.05]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.05]) ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        points(PrFix[U == 0.1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.1])  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        # Axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(1/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), 
                labels=c(0, expression(2/10^5), expression(2/N), expression(2/10^3), expression(2/10^2), expression(2/10^1)))
        # Plot labels etc.
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, expression("Neutral") , cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.02, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.02, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.3,  0.5,   expression(paste("Pr(fix | ", italic(x), ", SLR)")), cex=1.2, adj=c(0.5, 0.5), log='y', xpd=NA, srt=90)        


## Panel B: Beneficial Inversions
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,2.2), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(2*PrCatchSLR(SLR_pos=SLR_pos, x = x) ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points(PrFix[U == 0.02]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.02])/sI[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        points(PrFix[U == 0.05]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.05])/sI[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        points(PrFix[U == 0.1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.1])/sI[U == 0.1]   ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(0, 1, 2, 3), 
                labels=c(0, 
                         expression(italic(s[b])), 
                         expression(2*italic(s[b])), 
                         expression(3*italic(s[b]))))
        # Annnotations
        proportionalLabel( 1.2,  1.22,   expression(paste(SLR[pos], " = 1/2")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03, 1.075, expression(paste(B)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.085,   expression(paste("Beneficial")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(x))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## Panel C -- Additive SA Fitness Effects -- r = 0.002
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.00018), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        lines(2*sI_SA_vals$sI.vec[1]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = x) ~ x, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.02]) ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        points(PrFix[U == 0.05]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.05]) ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        points(PrFix[U == 0.1]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.1])  ~ x[U == 0.1], pch=24,  bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 1.2,   1.2,   expression(paste("SA Selection")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03,  1.075, expression(paste(C)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.002)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        


## Panel D -- Additive SA Fitness Effects -- r = 1/2
plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.0095), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        lines(2*sI_SA_vals$sI.vec[5]*A*x*exp(-A*x)*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = x) ~ x, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
#        abline(h=2*sI_SA_vals$sI.vec[5], lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02]*A*x[U == 0.02]*exp(-A*x[U == 0.02])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.02]) ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        points(PrFix[U == 0.05]*A*x[U == 0.05]*exp(-A*x[U == 0.05])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.05]) ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        points(PrFix[U == 0.1]*A*x[U == 0.1]*exp(-A*x[U == 0.1])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.1])  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 0.03, 1.075, expression(paste(D)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.5)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.5, -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        




### ROW 2: SLR_pos = 1/2
    SLR_pos  <-  SLR_positions[2]

## Panel E: Beneficial Inversions
    pFixNeutral  <-  pFixNeutralY(Ny=Ny)
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^7), 4/(10^3)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # pInv
        yvar=pFixNeutral*PrCatchSLR(SLR_pos=SLR_pos, x = x)
        yvar[1]=2/10^7
        lines(yvar ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFix[U == 0.02]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.02]) ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        points(PrFix[U == 0.05]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.05]) ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        points(PrFix[U == 0.1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.1])  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data=neutral_N10k)
        # Axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^7, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), 
                labels=c(0, expression(2/10^5), expression(2/N), expression(2/10^3), expression(2/10^2), expression(2/10^1)))
        # Plot labels etc.
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.02, 0.24, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y', srt=30)
        proportionalLabel(-0.02, 0.22, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y', srt=30)
        proportionalLabel(-0.3,  0.5,   expression(paste("Pr(fix | ", italic(x), ", SLR)")), cex=1.2, adj=c(0.5, 0.5), log='y', xpd=NA, srt=90)        
        proportionalLabel( 0.5, -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), log='y', xpd=NA)


## Panel F: Beneficial Inversions
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,2.2), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(2*PrCatchSLR(SLR_pos=SLR_pos, x = x) ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points(PrFix[U == 0.02]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.02])/sI[U == 0.02] ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        points(PrFix[U == 0.05]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.05])/sI[U == 0.05] ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        points(PrFix[U == 0.1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.1])/sI[U == 0.1]   ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data=ben_N10k)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(0, 1, 2, 3), 
                labels=c(0, 
                         expression(italic(s[b])), 
                         expression(2*italic(s[b])), 
                         expression(3*italic(s[b]))))
        # Annnotations
        proportionalLabel( 1.2,  1.2,   expression(paste(SLR[pos], " = 1/10")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03, 1.075, expression(paste(F)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(x))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        proportionalLabel( 0.5, -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## Panel G -- Additive SA Fitness Effects -- r = 0.002
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.00018), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        lines(2*sI_SA_vals$sI.vec[1]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = x) ~ x, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.02]) ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        points(PrFix[U == 0.05]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.05]) ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        points(PrFix[U == 0.1]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.1])  ~ x[U == 0.1], pch=24,  bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r1)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(G)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.002)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5, -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        


## Panel H -- Additive SA Fitness Effects -- r = 1/2
plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.0095), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # deterministic sI_SA
        lines(2*sI_SA_vals$sI.vec[5]*A*x*exp(-A*x)*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = x) ~ x, lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
#        abline(h=2*sI_SA_vals$sI.vec[5], lwd=2, lty=2, col=transparentColor(COLS[1], opacity=1.0))
        # pInv
        points(PrFix[U == 0.02]*A*x[U == 0.02]*exp(-A*x[U == 0.02])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.02]) ~ x[U == 0.02], pch=21, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        points(PrFix[U == 0.05]*A*x[U == 0.05]*exp(-A*x[U == 0.05])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.05]) ~ x[U == 0.05], pch=22, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        points(PrFix[U == 0.1]*A*x[U == 0.1]*exp(-A*x[U == 0.1])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = x[U == 0.1])  ~ x[U == 0.1],  pch=24, bg=transparentColor(COLS[4], opacity=0.8), data = SA_add_r5)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Annnotations
        proportionalLabel( 0.03, 1.075, expression(paste(H)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,   1.075,   expression(paste(italic(r[SA])==0.5)), cex=1.25, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5, -0.3,   expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.8,
               y       =  usr[4]*0.98,
               legend  =  c(
                            expression("Ignoring del. mut.")),
               lty     =  2,
               lwd     =  2,
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )        # Legend
        legend(
               x       =  usr[2]*0.5,
               y       =  usr[4]*0.875,
               legend  =  c(
                            expression(paste(italic(U/hs), " = ", 8)),
                            expression(paste(italic(U/hs), " = ", 20)),
                            expression(paste(italic(U/hs), " = ", 40))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg     =  transparentColor(COLS[4], opacity=0.8),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


}





#' Expected overall distribution of fixed inversions 
#' parameter Ud can take values of 0.02, 0.05, 0.1
expectedInvDistFig  <-  function(Ud = 0.05) {


    ## import data.frames
    # Neutral sims
    neutral_N10k  <-  read.csv(file = "./output/data/simResults/shelter_PrFixFig_h0.25_s0.01_N10k_deterministic_q.csv")

    # beneficial sims
    ben_N10k <-  read.csv(file = "./output/data/simResults/beneficial_PrFixFig_sI0.02_h0.25_s0.01_N10k_deterministic_q.csv")

    ## SA selection
    SA_add_r1  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.002_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")
    SA_add_r5  <-  read.csv(file = "./output/data/simResults/SA_PrFixFig_rSA0.498_sf0.05_sm0.05_hSA0.5_h0.25_s0.01_NumEqFreq_N10k.csv")

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey90'))
    COLS     <-  colfunc(5)

    # set constants
    x              <-  seq(0,1,by=0.005)
    Ny             <-  10^4/2
    A              <-  1
    P              <-  0.05
    lambda         <-  10
    SLR_positions  <-  c(0.5, 0.1)

    # Calculate deterministic s_I for SA selection
    log(0.002)
    log(0.498)
    logSeq  <-  seq(from=-6.214608, to=-0.6971552, len=5)
    rSeq  <-  round(exp(logSeq), digits=3)
    rSeq[4]  <-  rSeq[4] - 0.001
    eqFreq.dat  <-  read.csv('./output/data/simResults/SA-EqFreqs-EqualSel_sf0.05_sm0.05_hf0.5_hm0.5.csv', header=TRUE)
        # subset sI data to these values
        eqFreq.subdat  <-  eqFreq.dat[eqFreq.dat$r %in% rSeq,]
        sI.vec  <-  c()
        for(i in 1:nrow(eqFreq.subdat)) {
            sI.vec[i]  <-  sI_SA(sm = eqFreq.subdat$sm[i], hm = eqFreq.subdat$hm[i], Yhat = eqFreq.subdat$Y[i], Xfhat = eqFreq.subdat$Xf[i])
        }
        sI_SA_vals  <-  data.frame(cbind(eqFreq.subdat, sI.vec))


    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

### ROW 1: SLR_pos = 1/2
    SLR_pos  <-  SLR_positions[1]

    # Results for Random Breakpoint Model    
    neutralRBP  <-  neutral_N10k$PrFix[neutral_N10k$U == Ud]*PrCatchSLR(SLR_pos=SLR_pos, x = neutral_N10k$x[neutral_N10k$U == Ud])*distNewInvRBP(x=neutral_N10k$x[neutral_N10k$U == Ud])
    neutralRBP  <-  neutralRBP / sum(neutralRBP)
    beneficialRBP  <-  (ben_N10k$PrFix[ben_N10k$U == Ud]*PrCatchSLR(SLR_pos=SLR_pos, x = ben_N10k$x[ben_N10k$U == Ud])/ben_N10k$sI[ben_N10k$U == Ud])*distNewInvRBP(x=ben_N10k$x[ben_N10k$U == Ud])
    beneficialRBP  <-  beneficialRBP / sum(beneficialRBP)
    SAlinked_RBP  <-  SA_add_r1$PrFix[SA_add_r1$U == Ud]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = SA_add_r1$x[SA_add_r1$U == Ud])*distNewInvRBP(x=SA_add_r1$x[SA_add_r1$U == Ud])
    SAlinked_RBP  <-  SAlinked_RBP / sum(SAlinked_RBP)
    SAunlinked_RBP  <-  (SA_add_r5$PrFix[SA_add_r5$U == Ud]*A*SA_add_r5$x[SA_add_r5$U == Ud]*exp(-A*SA_add_r5$x[SA_add_r5$U == Ud])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = SA_add_r5$x[SA_add_r5$U == Ud])*distNewInvRBP(x=SA_add_r5$x[SA_add_r5$U == Ud]))
    SAunlinked_RBP  <-  SAunlinked_RBP / sum(SAunlinked_RBP)

    # Results for Exponential Model    
    neutralEXP  <-  neutral_N10k$PrFix[neutral_N10k$U == Ud]*PrCatchSLR(SLR_pos=SLR_pos, x = neutral_N10k$x[neutral_N10k$U == Ud])*distNewInvEXP(x=neutral_N10k$x[neutral_N10k$U == Ud], lambda=10)
    neutralEXP  <-  neutralEXP / sum(neutralEXP)
    beneficialEXP  <-  (ben_N10k$PrFix[ben_N10k$U == Ud]*PrCatchSLR(SLR_pos=SLR_pos, x = ben_N10k$x[ben_N10k$U == Ud])/ben_N10k$sI[ben_N10k$U == Ud])*distNewInvEXP(x=ben_N10k$x[ben_N10k$U == Ud], lambda=10)
    beneficialEXP  <-  beneficialEXP / sum(beneficialEXP)
    SAlinked_EXP  <-  SA_add_r1$PrFix[SA_add_r1$U == Ud]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = SA_add_r1$x[SA_add_r1$U == Ud])*distNewInvEXP(x=SA_add_r1$x[SA_add_r1$U == Ud], lambda=10)
    SAlinked_EXP  <-  SAlinked_EXP / sum(SAlinked_EXP)
    SAunlinked_EXP  <-  (SA_add_r5$PrFix[SA_add_r5$U == Ud]*A*SA_add_r5$x[SA_add_r5$U == Ud]*exp(-A*SA_add_r5$x[SA_add_r5$U == Ud])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = SA_add_r5$x[SA_add_r5$U == Ud])*distNewInvEXP(x=SA_add_r5$x[SA_add_r5$U == Ud], lambda=10))
    SAunlinked_EXP  <-  SAunlinked_EXP / sum(SAunlinked_EXP)

    x  <-  neutral_N10k$x[neutral_N10k$U == Ud]

## Panel A: Random Breakpoint 
    par(omi=c(0.1, 0.4, 0.1, 0.4), mar = c(2.5,3,3,1.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.05*(max(neutralRBP,beneficialRBP,SAlinked_RBP,SAunlinked_RBP))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Probability Density
        points(neutralRBP ~ x,     pch=21, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=neutralRBP), pch=21, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(beneficialRBP ~ x,  pch=22, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=beneficialRBP), pch=22, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(SAlinked_RBP ~ x,   pch=24, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=SAlinked_RBP), pch=24, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(SAunlinked_RBP ~ x, pch=25, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=SAunlinked_RBP), pch=25, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 1.18,  1.275,  expression(paste(SLR[pos], " = 0.5")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03,  1.075, expression(paste(A)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.1,  expression(paste("Random Breakpoint ", italic(f)(italic(x)))), cex=1.3, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,    expression(paste("Probability density (scaled)")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        

## Panel B: Exponential Model
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.05*(max(neutralEXP,beneficialEXP,SAlinked_EXP,SAunlinked_EXP))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Prob. Density
        points(neutralEXP ~ x,     pch=21, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=neutralEXP), pch=21, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(beneficialEXP ~ x,  pch=22, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=beneficialEXP), pch=22, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(SAlinked_EXP ~ x,   pch=24, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=SAlinked_EXP), pch=24, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(SAunlinked_EXP ~ x, pch=25, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=SAunlinked_EXP), pch=25, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(B)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.1,  expression(paste("Exponential ", italic(f)(italic(x)))), cex=1.3, adj=c(0.5, 0.5), xpd=NA)


### ROW 2: SLR_pos = 1/10
    SLR_pos  <-  SLR_positions[2]

    # Results for Random Breakpoint Model    
    neutralRBP  <-  neutral_N10k$PrFix[neutral_N10k$U == Ud]*PrCatchSLR(SLR_pos=SLR_pos, x = neutral_N10k$x[neutral_N10k$U == Ud])*distNewInvRBP(x=neutral_N10k$x[neutral_N10k$U == Ud])
    neutralRBP  <-  neutralRBP / sum(neutralRBP)
    beneficialRBP  <-  (ben_N10k$PrFix[ben_N10k$U == Ud]*PrCatchSLR(SLR_pos=SLR_pos, x = ben_N10k$x[ben_N10k$U == Ud])/ben_N10k$sI[ben_N10k$U == Ud])*distNewInvRBP(x=ben_N10k$x[ben_N10k$U == Ud])
    beneficialRBP  <-  beneficialRBP / sum(beneficialRBP)
    SAlinked_RBP  <-  SA_add_r1$PrFix[SA_add_r1$U == Ud]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = SA_add_r1$x[SA_add_r1$U == Ud])*distNewInvRBP(x=SA_add_r1$x[SA_add_r1$U == Ud])
    SAlinked_RBP  <-  SAlinked_RBP / sum(SAlinked_RBP)
    SAunlinked_RBP  <-  (SA_add_r5$PrFix[SA_add_r5$U == Ud]*A*SA_add_r5$x[SA_add_r5$U == Ud]*exp(-A*SA_add_r5$x[SA_add_r5$U == Ud])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = SA_add_r5$x[SA_add_r5$U == Ud])*distNewInvRBP(x=SA_add_r5$x[SA_add_r5$U == Ud]))
    SAunlinked_RBP  <-  SAunlinked_RBP / sum(SAunlinked_RBP)

    # Results for Exponential Model    
    neutralEXP  <-  neutral_N10k$PrFix[neutral_N10k$U == Ud]*PrCatchSLR(SLR_pos=SLR_pos, x = neutral_N10k$x[neutral_N10k$U == Ud])*distNewInvEXP(x=neutral_N10k$x[neutral_N10k$U == Ud], lambda=10)
    neutralEXP  <-  neutralEXP / sum(neutralEXP)
    beneficialEXP  <-  (ben_N10k$PrFix[ben_N10k$U == Ud]*PrCatchSLR(SLR_pos=SLR_pos, x = ben_N10k$x[ben_N10k$U == Ud])/ben_N10k$sI[ben_N10k$U == Ud])*distNewInvEXP(x=ben_N10k$x[ben_N10k$U == Ud], lambda=10)
    beneficialEXP  <-  beneficialEXP / sum(beneficialEXP)
    SAlinked_EXP  <-  SA_add_r1$PrFix[SA_add_r1$U == Ud]*A*P*exp(-A*P)*sI_SA_vals$Y[1]*PrCatchSLR(SLR_pos=SLR_pos, x = SA_add_r1$x[SA_add_r1$U == Ud])*distNewInvEXP(x=SA_add_r1$x[SA_add_r1$U == Ud], lambda=10)
    SAlinked_EXP  <-  SAlinked_EXP / sum(SAlinked_EXP)
    SAunlinked_EXP  <-  (SA_add_r5$PrFix[SA_add_r5$U == Ud]*A*SA_add_r5$x[SA_add_r5$U == Ud]*exp(-A*SA_add_r5$x[SA_add_r5$U == Ud])*sI_SA_vals$Y[5]*PrCatchSLR(SLR_pos=SLR_pos, x = SA_add_r5$x[SA_add_r5$U == Ud])*distNewInvEXP(x=SA_add_r5$x[SA_add_r5$U == Ud], lambda=10))
    SAunlinked_EXP  <-  SAunlinked_EXP / sum(SAunlinked_EXP)


## Panel C: Random Breakpoint 
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.05*(max(neutralRBP,beneficialRBP,SAlinked_RBP,SAunlinked_RBP))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Prob Density
        points(neutralRBP ~ x,     pch=21, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=neutralRBP), pch=21, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(beneficialRBP ~ x,  pch=22, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=beneficialRBP), pch=22, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(SAlinked_RBP ~ x,   pch=24, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=SAlinked_RBP), pch=24, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(SAunlinked_RBP ~ x, pch=25, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=SAunlinked_RBP), pch=25, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 1.18,  1.2,   expression(paste(SLR[pos], " = 0.1")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03,  1.075, expression(paste(C)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,    expression(paste("Probability density (scaled)")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

## Panel D: Exponential Model
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.05*(max(neutralEXP,beneficialEXP,SAlinked_EXP,SAunlinked_EXP))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Prob. Density
        points(neutralEXP ~ x,     pch=21, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=neutralEXP), pch=21, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(beneficialEXP ~ x,  pch=22, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=beneficialEXP), pch=22, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(SAlinked_EXP ~ x,   pch=24, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=SAlinked_EXP), pch=24, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        points(SAunlinked_EXP ~ x, pch=25, bg=transparentColor(COLS[4], opacity=0.8), type='b', col=transparentColor('#252525', opacity=1))
        points(0 ~ weighted.mean(x=x,w=SAunlinked_EXP), pch=25, bg=transparentColor(COLS[1], opacity=1), type='b', col=transparentColor('#252525', opacity=1))
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(D)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Neutral")),
                            expression(paste("Beneficial")),
                            expression(paste("SA Sel. ", r[SA]==0.002)),
                            expression(paste("SA Sel. ", r[SA]==0.5)),
                            expression(paste("          Mean lengths"))),
               col     =  c(1,transparentColor('#252525', opacity=1)),
               pch     =  c(21,22,24,25,NA),
               pt.bg   =  transparentColor(COLS[4], opacity=0.8),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        proportionalLabel( px=c(0.43,0.48,0.53,0.58), py=c(0.615,0.615,0.615,0.615), pch=c(21,22,24,25), bg= transparentColor(COLS[1], opacity=1), cex=1, adj=c(0.5, 0.5), xpd=NA, text=FALSE)


}







##############################################################
##############################################################
##  Supplementary Figures 

PrSpanSLR_Fig  <-  function() {

    # Colors
    COLS     <-  c('#252525', 'white')

    # set constants
    x   <-  seq(0,1,by=0.005)
    SLR_positions  <-  c(0.5, 0.1)

    # set plot layout
    layout.mat <- matrix(c(1:2), nrow=1, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # calculate Pr(SLR | x)
    PrSLR_0.5  <-  PrCatchSLR(SLR_pos=SLR_positions[1], x = x)
    PrSLR_0.1  <-  PrCatchSLR(SLR_pos=SLR_positions[2], x = x)

## Panel A: SLR_loc = 1/2
    par(omi=c(0.25, 0.5, 0.5, 0.25), mar = c(4.5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0, 1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Shading
        polygon(x = c(0,SLR_positions[1], SLR_positions[1], 0), y = c(usr[4],usr[4],usr[1],usr[1]), col = transparentColor(COLS[1], opacity=0.2), border=FALSE)
        # pInv
        lines(PrSLR_0.5 ~ x, lty=1, lwd=2, col=transparentColor(COLS[1], opacity=1))
        # Axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.5, 1.3, expression(SLR[pos]==1/2) , cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3, 0.5,  expression(paste("Pr(SLR | ", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5, -0.3, expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
#        proportionalLabel( 0.25, 1.15, expression(italic(y)[1]==0.5) , cex=0.75, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.75, 1.15, expression(italic(y)[3]==0.5) , cex=0.75, adj=c(0.5, 0.5), xpd=NA)
        #Annotations
        proportionalLabel(0.072, 0.75, expression(italic(x)/(1 - italic(x))))
        proportionalRect(0.06, 0.635, 0.4, 0.76)
        proportionalArrows(px1=0.24, py1=0.6325, px2=0.34, py2=0.5, length=0.05, lwd=1.2)

        proportionalLabel(0.75, 0.73, expression(1))
        proportionalRect(0.73, 0.64, 0.81, 0.76)
        proportionalArrows(px1=0.77, py1=0.76, px2=0.77, py2=0.9, length=0.05, lwd=1.2)


## Panel B: SLR_loc = 1/10
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0, 1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Shading
        polygon(x = c(0,SLR_positions[2], SLR_positions[2], 0), y = c(usr[4],usr[4],usr[1],usr[1]), col = transparentColor(COLS[1], opacity=0.2), border=FALSE)
        polygon(x = c(SLR_positions[2], (1-SLR_positions[2]), (1-SLR_positions[2]), SLR_positions[2]), y = c(usr[4],usr[4],usr[1],usr[1]), col = transparentColor(COLS[1], opacity=0.1), border=FALSE)
        # pInv
        lines(PrSLR_0.1 ~ x, lty=1, lwd=2, col=transparentColor(COLS[1], opacity=1))
        # Axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.5, 1.3, expression(SLR[pos]==1/10) , cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5, -0.3, expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
#        proportionalLabel( 0.1, 1.15, expression(italic(y)[1]==0.1) , cex=0.75, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.575, 1.15, expression(italic(y)[3]==0.9) , cex=0.75, adj=c(0.5, 0.5), xpd=NA)
        #Annotations
        proportionalLabel(0.075, 0.45, expression(italic(x)/(1 - italic(x))))
#        proportionalRect(0.06, 0.335, 0.41, 0.46)
        proportionalArrows(px1=0.235, py1=0.335, px2=0.09, py2=0.12, length=0.05, lwd=1.2)

        proportionalLabel(0.4, 0.6, expression(italic(y[1])/(1 - italic(x))))
#        proportionalRect(0.39, 0.485, 0.76, 0.61)
        proportionalArrows(px1=0.575, py1=0.4825, px2=0.71, py2=0.39, length=0.05, lwd=1.2)

        proportionalLabel(0.905, 0.73, expression(1))
#        proportionalRect(0.885, 0.64, 0.965, 0.76)
        proportionalArrows(px1=0.925, py1=0.765, px2=0.925, py2=0.9, length=0.05, lwd=1.2)

}



##############################################################
##############################################################
##  Old Figures for additive del. mutations

fixationProbabilityFigure  <-  function() {

    ## import data.frames
    # Neutral sims
    neutralNy500    <-  read.csv(file = "./output/data/simResults/neutral-Y-pFix_Ny100_Ud0.2_sd0.05.csv")
    neutralNy1000   <-  read.csv(file = "./output/data/simResults/neutral-Y-pFix_Ny1000_Ud0.2_sd0.05.csv")
    neutralNy10000  <-  read.csv(file = "./output/data/simResults/neutral-Y-pFix_Ny10000_Ud0.2_sd0.05.csv")

    # beneficial sims
    benNy1000   <-  read.csv(file = "./output/data/simResults/beneficial-Y-pFix_Ny1000_sI0.02_Ud0.2_sd0.02.csv")
    benNy1000_2 <-  read.csv(file = "./output/data/simResults/beneficial-Y-pFix_Ny1000_sI0.02_Ud0.1_sd0.02.csv")

    # sexually antagonistic sims
    SAaddLinked      <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI_sf0.05_sm0.05_hf0.5_hm0.5.csv")
#    SAaddLinkedWF    <-  read.csv(file = "./output/data/simResults/SA-Y-pFix2_N10000_hf0.5_sf0.05_hm0.5_sm0.05_r0.0101_Ud0.2_sd0.02.csv")
#    SAaddUnLinkedWF  <-  read.csv(file = "./output/data/simResults/SA-Y-pFix2_N10000_hf0.5_sf0.05_hm0.5_sm0.05_r0.5_Ud0.2_sd0.02.csv")
    SAaddLinkedWF    <-  read.csv(file = "./output/data/simResults/SA-Y-pFix2_N10000_hf0.5_sf0.05_hm0.5_sm0.05_r0.0101_Ud0.1_sd0.01.csv")
    SAaddUnLinkedWF  <-  read.csv(file = "./output/data/simResults/SA-Y-pFix2_N10000_hf0.5_sf0.05_hm0.5_sm0.05_r0.5_Ud0.1_sd0.01.csv")


    # set constants
    x     <-  seq(0,1,by=0.005)
    Ny    <-  1*10^3
    sd    <-  0.02
    A     <-  1
    P     <-  0.05

    # set plot layout
    layout.mat <- matrix(c(1:3), nrow=1, ncol=3, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A: Neutral Inversions
    pFixNeutral  <-  pFixNeutralY(Ny=Ny)
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,5,5,2), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        abline(h=pFixNeutral*Ny, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points(PrFixSimDel*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.1), data=neutralNy500)
        points(PrFixSimDel*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.3), data=neutralNy1000)
        points(PrFixSimDel*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.7), data=neutralNy10000)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(0, 1/4, 1/2, 3/4, 1), labels = NA)
        labelPos  <-  axTicks(2) + 0.035
        proportionalLabel(-0.14,  labelPos[5],   expression(paste("1/", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[4],   expression(paste("3/4", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[3],   expression(paste("1/2", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[2],   expression(paste("1/4", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[1],   expression(paste("0.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste("Neutral")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Pr(fix | ", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.37,
               y       =  usr[4]*0.325,
               legend  =  c(
                            expression(paste(italic(N[Y]), " = ",10^4)),
                            expression(paste(italic(N[Y]), " = ",10^3)),
                            expression(paste(italic(N[Y]), " = ",10^2))),
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
    Ud  <-  0.2
    sb  <-  0.02
    sd  <-  0.02
    pFixBeneficial  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)
    Ud  <-  0.1
    pFixBeneficial_2  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)
#    Ud  <-  0.2
#    sd  <-  0.03
#    pFixBeneficial_3  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.1*max(pFixBeneficial/sb)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(pFixBeneficial/sb ~ x, lty=1, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points((PrFixSim*PrMutFree)/sb ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.7), data=benNy1000)
        lines(pFixBeneficial_2/sb ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points((PrFixSim*PrMutFree)/sb ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.1), data=benNy1000_2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        labelPos  <-  axTicks(2)/2.4 + 0.035
        proportionalLabel(-0.12,  labelPos[5],   expression(paste("2",italic(s[I]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.12,  labelPos[4],   expression(paste("3", italic(s[I]),"/4")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.12,  labelPos[3],   expression(paste(italic(s[I]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.12,  labelPos[2],   expression(paste(italic(s[I]),"/2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.12,  labelPos[1],   expression(paste("0.0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        # Annnotations
        proportionalLabel( 0.03,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste("Beneficial")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Plot annotations
        # Legend
        legend(
               x       =  usr[2]*0.76,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("")),
                            expression(paste(""))),
               lty     =  c(2,1),
               seg.len =  2,
               lwd     =  2,
               col     =  c(transparentColor('#252525', opacity=1.0)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2]*0.975,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(~U[d]), " = ",0.1)),
                            expression(paste(italic(~U[d]), " = ",0.2))),
               pch     =  c(21),
               pt.bg   =  c(transparentColor('#252525', opacity=0.1),
                            transparentColor('#252525', opacity=0.7)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel C: Sexual Antagonism (2-locus models)
    Ud    <-  0.2
    sf    <-  0.05
    sm    <-  0.05
    sd    <-  0.02
    r_SA  <-  SAaddLinked$r[6]
    qHatSA             <-  qHat_SAAdd(sf=sf, sm=sm)
    pFixSAunlinkedAdd  <-  pFix_SAUnlinkedAdd(x=x, sm=sm, qHat=qHatSA, A=A, Ud=Ud, sd=sd)
    sISAaddLinked      <-  SAaddLinked$sI[SAaddLinked$r == r_SA]
    pFixSAlinked       <-  pFix_SALinkedAdd(x=x, sI=sISAaddLinked, A=A, P=P, Ud=Ud, sd=sd)
    yMaxLinked         <-  1.1*max(c(pFixSAlinked,SAaddLinkedWF$PrFixSimDel))
    yMaxUnLinked       <-  1.1*max(c(pFixSAunlinkedAdd,SAaddUnLinkedWF$PrFixSimDel))
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,yMaxUnLinked), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv Linked
        lines(pFixSAlinked ~ x, lty=1, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFixSimDel/invSizes ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.75), data=SAaddLinkedWF)
        axis(2, las=1, at=axTicks(2), labels=sciNotation(axTicks(2),1))
        # pInv Unlinked
        lines(pFixSAunlinkedAdd ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFixSimDel/invSizes ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.25), data=SAaddUnLinkedWF)
        # x Axis
        axis(1, las=1)
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste("Sex Antagonism")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Annotations
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2]*0.76,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("")),
                            expression(paste(""))),
               lty     =  c(2,1),
               seg.len =  2,
               lwd     =  2,
               col     =  c(transparentColor('#252525', opacity=1.0)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2]*0.975,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(~~italic(r)," = 1/2")),
                            expression(paste(~~italic(r)," = 0.01"))),
               pch     =  c(21),
               pt.bg   =  c(transparentColor('#252525', opacity=0.1),
                            transparentColor('#252525', opacity=0.7)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


}


PrCatchFixFigure  <-  function(SLR_pos = c(1/2, 1/10)) {

    ## import data.frames
    # Neutral sims
    neutralNy500    <-  read.csv(file = "./output/data/simResults/neutral-Y-pFix_Ny100_Ud0.2_sd0.05.csv")
    neutralNy1000   <-  read.csv(file = "./output/data/simResults/neutral-Y-pFix_Ny1000_Ud0.2_sd0.05.csv")
    neutralNy10000  <-  read.csv(file = "./output/data/simResults/neutral-Y-pFix_Ny10000_Ud0.2_sd0.05.csv")

    # beneficial sims
    benNy1000   <-  read.csv(file = "./output/data/simResults/beneficial-Y-pFix_Ny1000_sI0.02_Ud0.2_sd0.02.csv")
    benNy1000_2 <-  read.csv(file = "./output/data/simResults/beneficial-Y-pFix_Ny1000_sI0.02_Ud0.1_sd0.02.csv")

    # sexually antagonistic sims
    SAaddLinked      <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI_sf0.05_sm0.05_hf0.5_hm0.5.csv")
#    SAaddLinkedWF    <-  read.csv(file = "./output/data/simResults/SA-Y-pFix2_N10000_hf0.5_sf0.05_hm0.5_sm0.05_r0.0101_Ud0.2_sd0.02.csv")
#    SAaddUnLinkedWF  <-  read.csv(file = "./output/data/simResults/SA-Y-pFix2_N10000_hf0.5_sf0.05_hm0.5_sm0.05_r0.5_Ud0.2_sd0.02.csv")
    SAaddLinkedWF    <-  read.csv(file = "./output/data/simResults/SA-Y-pFix2_N10000_hf0.5_sf0.05_hm0.5_sm0.05_r0.0101_Ud0.1_sd0.01.csv")
    SAaddUnLinkedWF  <-  read.csv(file = "./output/data/simResults/SA-Y-pFix2_N10000_hf0.5_sf0.05_hm0.5_sm0.05_r0.5_Ud0.1_sd0.01.csv")


    # set constants
    x       <-  seq(0,1,by=0.005)
    Ny      <-  1*10^3
    sd      <-  0.02
    A       <-  1
    P       <-  0.05
    

    # set plot layout
    layout.mat <- matrix(c(1:6), nrow=2, ncol=3, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

### ROW 1: SLR_pos = 1/2
    SLR_pos  <-  SLR_pos[1]
## Panel A: Neutral Inversions
    pFixNeutral  <-  pFixNeutralY(Ny=Ny)
    par(omi=c(0.1, 0.5, 0.1, 0.5), mar = c(4,5,4,2), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(PrCatchSLR(SLR_pos=SLR_pos, x = x) ~ x, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.1), data=neutralNy500)
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.3), data=neutralNy1000)
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.7), data=neutralNy10000)
        # axes
        axis(1, labels=NA)
        axis(2, las=1, at=c(0, 1/4, 1/2, 3/4, 1), labels = NA)
        labelPos  <-  axTicks(2) + 0.035
        proportionalLabel(-0.14,  labelPos[5],   expression(paste("1/", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[4],   expression(paste("3/4", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[3],   expression(paste("1/2", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[2],   expression(paste("1/4", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[1],   expression(paste("0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste("Neutral")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Pr(fix | ", italic(x), ", SDR)")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(x))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.37,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(N[Y]), " = ",10^4)),
                            expression(paste(italic(N[Y]), " = ",10^3)),
                            expression(paste(italic(N[Y]), " = ",10^2))),
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



## Panel B: Beneficial Inversions
    Ud  <-  0.2
    sb  <-  0.02
    sd  <-  0.02
    pFixBeneficial  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)
    Ud  <-  0.1
    pFixBeneficial_2  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.1*max((benNy1000_2$PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = benNy1000_2$invSizes)/sb),(pFixBeneficial_2*PrCatchSLR(SLR_pos=SLR_pos, x = x)/sb))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(pFixBeneficial*PrCatchSLR(SLR_pos=SLR_pos, x = x)/sb ~ x, lty=1, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)/sb ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.7), data=benNy1000)
        lines(pFixBeneficial_2*PrCatchSLR(SLR_pos=SLR_pos, x = x)/sb ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)/sb ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.1), data=benNy1000_2)
        # axes
        axis(1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        labelPos  <-  axTicks(2)*2 + 0.035
        proportionalLabel(-0.13,  labelPos[5],   expression(paste("2",italic(s[I]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.13,  labelPos[4],   expression(paste("3", italic(s[I]),"/4")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.13,  labelPos[3],   expression(paste(italic(s[I]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.13,  labelPos[2],   expression(paste(italic(s[I]),"/2")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.13,  labelPos[1],   expression(paste("0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        # Annnotations
        proportionalLabel( 0.5,  1.25,   expression(paste(SDR[loc], " = 1/2")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.085,   expression(paste("Beneficial")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(x))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.785,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("")),
                            expression(paste(""))),
               lty     =  c(2,1),
               seg.len =  2,
               lwd     =  2,
               col     =  c(transparentColor('#252525', opacity=1.0)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2]*0.975,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(~U[d]), " = ",0.1)),
                            expression(paste(italic(~U[d]), " = ",0.2))),
               pch     =  c(21),
               pt.bg   =  c(transparentColor('#252525', opacity=0.1),
                            transparentColor('#252525', opacity=0.7)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel C: Sexual Antagonism (2-locus models)
    Ud    <-  0.2
    sf    <-  0.05
    sm    <-  0.05
    sd    <-  0.02
    r_SA  <-  SAaddLinked$r[6]
    qHatSA             <-  qHat_SAAdd(sf=sf, sm=sm)
    pFixSAunlinkedAdd  <-  pFix_SAUnlinkedAdd(x=x, sm=sm, qHat=qHatSA, A=A, Ud=Ud, sd=sd)*PrCatchSLR(SLR_pos=SLR_pos, x = x)
    sISAaddLinked      <-  SAaddLinked$sI[SAaddLinked$r == r_SA]
    pFixSAlinked       <-  pFix_SALinkedAdd(x=x, sI=sISAaddLinked, A=A, P=P, Ud=Ud, sd=sd)*PrCatchSLR(SLR_pos=SLR_pos, x = x)
    yMaxLinked         <-  1.1*max(c(pFixSAlinked,(SAaddLinkedWF$PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = SAaddLinkedWF$invSizes)/SAaddLinkedWF$invSizes)))
    yMaxUnLinked       <-  1.1*max(c(pFixSAunlinkedAdd,(SAaddUnLinkedWF$PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = SAaddUnLinkedWF$invSizes)/SAaddUnLinkedWF$invSizes)))
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,yMaxUnLinked), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv Linked
        lines(pFixSAlinked ~ x, lty=1, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)/invSizes ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.75), data=SAaddLinkedWF)
        axis(2, las=1, at=axTicks(2), labels=sciNotation(axTicks(2),1))
        # pInv Unlinked
        lines(pFixSAunlinkedAdd ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)/invSizes ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.25), data=SAaddUnLinkedWF)
        # x Axis
        axis(1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.085,   expression(paste("Sex Antagonism")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(x))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2]*0.785,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("")),
                            expression(paste(""))),
               lty     =  c(2,1),
               seg.len =  2,
               lwd     =  2,
               col     =  c(transparentColor('#252525', opacity=1.0)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2]*0.975,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(~~italic(r)," = 1/2")),
                            expression(paste(~~italic(r)," = 0.01"))),
               pch     =  c(21),
               pt.bg   =  c(transparentColor('#252525', opacity=0.1),
                            transparentColor('#252525', opacity=0.7)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


### ROW 2: SLR_pos = 1/2
    SLR_pos  <-  SLR_pos[2]

## Panel D: Neutral Inversions
    pFixNeutral  <-  pFixNeutralY(Ny=Ny)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(PrCatchSLR(SLR_pos=SLR_pos, x = x) ~ x, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.1), data=neutralNy500)
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.3), data=neutralNy1000)
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)*Ny ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.7), data=neutralNy10000)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(0, 1/4, 1/2, 3/4, 1), labels = NA)
        labelPos  <-  axTicks(2) + 0.035
        proportionalLabel(-0.14,  labelPos[5],   expression(paste("1/", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[4],   expression(paste("3/4", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[3],   expression(paste("1/2", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[2],   expression(paste("1/4", italic(N[y]))), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.14,  labelPos[1],   expression(paste("0")), cex=1, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Pr(fix | ", italic(x), ", SDR)")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Inversion size (", italic(x), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.37,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(N[Y]), " = ",10^4)),
                            expression(paste(italic(N[Y]), " = ",10^3)),
                            expression(paste(italic(N[Y]), " = ",10^2))),
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



## Panel E: Beneficial Inverions
    Ud  <-  0.2
    sb  <-  0.02
    sd  <-  0.02
    pFixBeneficial  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)
    Ud  <-  0.1
    pFixBeneficial_2  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.1*max((benNy1000_2$PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = benNy1000_2$invSizes)/sb),(pFixBeneficial_2*PrCatchSLR(SLR_pos=SLR_pos, x = x)/sb))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(pFixBeneficial*PrCatchSLR(SLR_pos=SLR_pos, x = x)/sb ~ x, lty=1, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)/sb ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.7), data=benNy1000)
        lines(pFixBeneficial_2*PrCatchSLR(SLR_pos=SLR_pos, x = x)/sb ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1.0))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)/sb ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.1), data=benNy1000_2)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=axTicks(2), labels=c(
                                expression(paste("0")),
                                expression(paste(italic(s[I]),"/2")),
                                expression(paste(italic(s[I]))),
                                expression(paste("3", italic(s[I]),"/4"))))
        # Annnotations
        proportionalLabel( 0.5,  1.25,   expression(paste(SDR[loc], " = 1/10")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03,  1.075, expression(paste(bold(E))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Inversion size (", italic(x), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.785,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("")),
                            expression(paste(""))),
               lty     =  c(2,1),
               seg.len =  2,
               lwd     =  2,
               col     =  c(transparentColor('#252525', opacity=1.0)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2]*0.975,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(~U[d]), " = ",0.1)),
                            expression(paste(italic(~U[d]), " = ",0.2))),
               pch     =  c(21),
               pt.bg   =  c(transparentColor('#252525', opacity=0.1),
                            transparentColor('#252525', opacity=0.7)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
## Panel F: Sexual Antagonism (2-locus models)
    Ud    <-  0.2
    sf    <-  0.05
    sm    <-  0.05
    sd    <-  0.02
    r_SA  <-  SAaddLinked$r[6]
    qHatSA             <-  qHat_SAAdd(sf=sf, sm=sm)
    pFixSAunlinkedAdd  <-  pFix_SAUnlinkedAdd(x=x, sm=sm, qHat=qHatSA, A=A, Ud=Ud, sd=sd)*PrCatchSLR(SLR_pos=SLR_pos, x = x)
    sISAaddLinked      <-  SAaddLinked$sI[SAaddLinked$r == r_SA]
    pFixSAlinked       <-  pFix_SALinkedAdd(x=x, sI=sISAaddLinked, A=A, P=P, Ud=Ud, sd=sd)*PrCatchSLR(SLR_pos=SLR_pos, x = x)
    yMaxLinked         <-  1.1*max(c(pFixSAlinked,(SAaddLinkedWF$PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = SAaddLinkedWF$invSizes)/SAaddLinkedWF$invSizes)))
    yMaxUnLinked       <-  1.1*max(c(pFixSAunlinkedAdd,(SAaddUnLinkedWF$PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = SAaddUnLinkedWF$invSizes)/SAaddUnLinkedWF$invSizes)))
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,yMaxUnLinked), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv Linked
        lines(pFixSAlinked ~ x, lty=1, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)/invSizes ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.75), data=SAaddLinkedWF)
        axis(2, las=1, at=axTicks(2), labels=sciNotation(axTicks(2),1))
        # pInv Unlinked
        lines(pFixSAunlinkedAdd ~ x, lty=2, lwd=2, col=transparentColor('#252525', opacity=1))
        points(PrFixSimDel*PrCatchSLR(SLR_pos=SLR_pos, x = invSizes)/invSizes ~ invSizes, pch=21, bg=transparentColor('#252525', opacity=0.25), data=SAaddUnLinkedWF)
        # x Axis
        axis(1, las=1)
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(F))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Inversion size (", italic(x), ")")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2]*0.785,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("")),
                            expression(paste(""))),
               lty     =  c(2,1),
               seg.len =  2,
               lwd     =  2,
               col     =  c(transparentColor('#252525', opacity=1.0)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2]*0.975,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(~~italic(r)," = 1/2")),
                            expression(paste(~~italic(r)," = 0.01"))),
               pch     =  c(21),
               pt.bg   =  c(transparentColor('#252525', opacity=0.1),
                            transparentColor('#252525', opacity=0.7)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}

determInvFreqPlot  <-  function(wesMovie = 'Zissou1') {

    ## import data.frames
    yInvUnLinked  <-  read.csv(file = "./output/data/simResults/deterministicSims-Y_sMax1_r0.5_hf0.5_hm0.5.csv")
    yInvLinked    <-  read.csv(file = "./output/data/simResults/deterministicSims-Y_sMax1_r0.02_hf0.5_hm0.5.csv")
    xInvUnLinked  <-  read.csv(file = "./output/data/simResults/deterministicSims-X_sMax1_r0.5_hf0.5_hm0.5.csv")
    xInvLinked    <-  read.csv(file = "./output/data/simResults/deterministicSims-X_sMax1_r0.02_hf0.5_hm0.5.csv")

    COLS  <-  wes_palette(wesMovie, 20, type="continuous")

## Panel A: Y-linked inversion Unlinked
plot.new()
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,5,5,2), bty='o', xaxt='s', yaxt='s')    
        usr      <-  par('usr')
        x        <-  yInvUnLinked$sm
        y        <-  yInvUnLinked$sf
        z        <-  round(yInvUnLinked$inv.eq, digits=3)
        xcoords  <-  unique(x)
        ycoords  <-  unique(y)
        surface.matrix = matrix(z,nrow=length(xcoords),ncol=length(ycoords),byrow=FALSE)
par(new = "TRUE",plt = c(0.1,0.4,0.60,0.95),las = 1,cex.axis = 1, bty='o', xaxt='s', yaxt='s',xpd = NA)
filled.contour3(xcoords, ycoords, surface.matrix, xlim = c(min(xcoords),max(xcoords)),ylim = c(min(ycoords),max(ycoords)),zlim = c(min(surface.matrix),max(surface.matrix)),
                col=COLS, axes=FALSE, xlab = "",ylab = "")
axis(2,labels = TRUE)
axis(1,labels = NA)
box()
proportionalLabel( 0.025,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel( 0.5,  1.2,   expression(paste(italic(r), " = 0.5")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel( -0.45,  0.5,   expression(paste("Y inversion")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
proportionalLabel( -0.3,  0.5,   expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        

## Panel B: Y-linked inversion r = 0.002
        x        <-  yInvLinked$sm
        y        <-  yInvLinked$sf
        z        <-  round(yInvLinked$inv.eq, digits = 3)
        xcoords  <-  unique(x)
        ycoords  <-  unique(y)
        surface.matrix = matrix(z,nrow=length(xcoords),ncol=length(ycoords),byrow=FALSE)
par(new = "TRUE",plt = c(0.5,0.8,0.60,0.95),las = 1, cex.axis = 1, bty='o', xaxt='s', yaxt='s',xpd = NA)
filled.contour3(xcoords,ycoords,surface.matrix,xlab = "",ylab = "",xlim = c(min(xcoords),max(xcoords)),ylim = c(min(ycoords),max(ycoords)),zlim = c(min(surface.matrix),max(surface.matrix)), 
    col=COLS, axes=FALSE)
axis(2,labels = NA)
axis(1,labels = NA)
box()
proportionalLabel( 0.025,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel( 0.5,  1.2,   expression(paste(italic(r), " = ", 0.02)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel( 1.6,  0.5,   expression(paste(hat(italic(Y)))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

## Panel C: X-linked inversion Unlinked
        x        <-  xInvUnLinked$sm
        y        <-  xInvUnLinked$sf
        z        <-  round(xInvUnLinked$inv.eq, digits = 3)
        xcoords  <-  unique(x)
        ycoords  <-  unique(y)
        surface.matrix = matrix(z,nrow=length(xcoords),ncol=length(ycoords),byrow=FALSE)
par(new = "TRUE",plt = c(0.1,0.4,0.15,0.5),las = 1, cex.axis = 1, bty='o', xaxt='s', yaxt='s',xpd = NA)
filled.contour3(xcoords,ycoords,surface.matrix, levels = fibonacci.scale(20), col=COLS, xlab = "",ylab = "",xlim = c(min(xcoords),max(xcoords)),ylim = c(min(ycoords),max(ycoords)),zlim = c(min(surface.matrix),max(surface.matrix)), axes=FALSE)
axis(2,labels = TRUE)
axis(1,labels = TRUE)
box()
proportionalLabel( 0.025,  1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel( -0.45,  0.5,   expression(paste("X inversion")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
proportionalLabel( -0.3,  0.5,   expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
proportionalLabel( 0.5,  -0.25,  expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

#Bottom right plot:
        x        <-  xInvLinked$sm
        y        <-  xInvLinked$sf
        z        <-  round(xInvLinked$inv.eq, digits = 3)
        xcoords  <-  unique(x)
        ycoords  <-  unique(y)
        surface.matrix = matrix(z,nrow=length(xcoords),ncol=length(ycoords),byrow=FALSE)
par(new = "TRUE",plt = c(0.5,0.8,0.15,0.5),las = 1,cex.axis = 1, bty='o', xaxt='s', yaxt='s',xpd = NA)
filled.contour3(xcoords,ycoords,surface.matrix, levels = fibonacci.scale(20), col=COLS, xlab = "",ylab = "",xlim = c(min(xcoords),max(xcoords)),ylim = c(min(ycoords),max(ycoords)),zlim = c(min(surface.matrix),max(surface.matrix)), axes=FALSE)
axis(2,labels = NA)
axis(1,labels = TRUE)
box()
proportionalLabel( 0.025,  1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel( 0.5,  -0.25,  expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
proportionalLabel( 1.6,  0.5,   expression(paste(hat(italic(X)))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

#Add a legend:
par(new = "TRUE",plt = c(0.85,0.9,0.60,0.95),las = 1,cex.axis = 1)
filled.legend(xcoords,ycoords,surface.matrix,col = COLS,xlab = "",ylab = "",xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),zlim = c(0,1))
par(new = "TRUE",plt = c(0.85,0.9,0.15,0.5),las = 1,cex.axis = 1)
filled.legend(xcoords,ycoords,surface.matrix,levels = fibonacci.scale(20),col = COLS,xlab = "",ylab = "",xlim = c(min(xintercepts),max(xintercepts)),ylim = c(min(slopes),max(slopes)),zlim = c(0,1))


}








recombEffectSATwoLocus  <-  function() {

    ## import data.frames
    # Equal Selection
    EqualSel1    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-EqualSel0.01_sm0.01_hf0.5_hm0.5.csv")
    EqualSel2    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-EqualSel0.0125_sm0.0125_hf0.5_hm0.5.csv")
    EqualSel3    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-EqualSel0.0166667_sm0.0166667_hf0.5_hm0.5.csv")
    EqualSel4    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-EqualSel0.025_sm0.025_hf0.5_hm0.5.csv")
    EqualSel5    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-EqualSel0.05_sm0.05_hf0.5_hm0.5.csv")

    # Biased Selection
    MaleBiasedSel1    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-MaleBiasedSel0.01_sm0.010101_hf0.5_hm0.5.csv")
    MaleBiasedSel2    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-MaleBiasedSel0.0125_sm0.0126582_hf0.5_hm0.5.csv")
    MaleBiasedSel3    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-MaleBiasedSel0.0166667_sm0.0169492_hf0.5_hm0.5.csv")
    MaleBiasedSel4    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-MaleBiasedSel0.025_sm0.025641_hf0.5_hm0.5.csv")
    MaleBiasedSel5    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-MaleBiasedSel0.05_sm0.0526316_hf0.5_hm0.5.csv")

    FemaleBiasedSel1    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-BiasedSel0.00990099_sm0.01_hf0.5_hm0.5.csv")
    FemaleBiasedSel2    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-BiasedSel0.0123457_sm0.0125_hf0.5_hm0.5.csv")
    FemaleBiasedSel3    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-BiasedSel0.0163934_sm0.0166667_hf0.5_hm0.5.csv")
    FemaleBiasedSel4    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-BiasedSel0.0243902_sm0.025_hf0.5_hm0.5.csv")
    FemaleBiasedSel5    <-  read.csv(file = "./output/data/simResults/SA-TwoLocus-Add-Linked-sI-BiasedSel0.047619_sm0.05_hf0.5_hm0.5.csv")

    # set plot layout
    layout.mat <- matrix(c(1:2), nrow=1, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A: Neutral Inversions
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,5,5,2), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1/2), ylim = c(0,(1.05*0.013)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # sI
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=0.2), data=EqualSel1)
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=0.4), data=EqualSel2)
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=0.6), data=EqualSel3)
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=0.8), data=EqualSel4)
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=1.0), data=EqualSel5)
        # axes
        axis(1, labels=TRUE)
        axis(2, labels=TRUE, at=c(0,0.004,0.008,0.012))
#        axis(2, las=1, at=axTicks(2), labels=sciNotation(axTicks(2),1))
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.2,   expression(paste("Equal Selection")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(italic(s[f]), " = ", italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel(-0.25,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(r))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)



## Panel B: Beneficial Inverions
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1/2), ylim = c(0,(1.05*0.013)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # sI
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=0.2), data=MaleBiasedSel1)
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=0.4), data=MaleBiasedSel2)
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=0.6), data=MaleBiasedSel3)
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=0.8), data=MaleBiasedSel4)
        lines(sI ~ r, lwd=2, lty=1, col=transparentColor('#252525', opacity=1.0), data=MaleBiasedSel5)
        # axes
        axis(1, labels=TRUE)
        axis(2, labels=NA, at=c(0,0.004,0.008,0.012))
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.2,   expression(paste("Biased Selection")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.075, expression(paste(italic(s[f]), " = ", italic(s[m]), "/(1 - ", italic(s[m]), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)        
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(r))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(s[m]), " = ",1/20)),
                            expression(paste(italic(s[m]), " = ",1/40)),
                            expression(paste(italic(s[m]), " = ",1/60)),
                            expression(paste(italic(s[m]), " = ",1/80)),
                            expression(paste(italic(s[m]), " = ",1/100))),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor('#252525', opacity=1.0),
                            transparentColor('#252525', opacity=0.8),
                            transparentColor('#252525', opacity=0.6),
                            transparentColor('#252525', opacity=0.4),
                            transparentColor('#252525', opacity=0.2)),
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






expectedDistributionFig  <-  function(SLR_pos = 1/2) {

    # set constants
    x     <-  seq(0,1,by=0.005)
    Ny    <-  1*10^3
    Ud    <-  0.2
    sd    <-  0.02
    A     <-  1
    P     <-  0.05
    sb    <-  0.02
    UdSA    <-  0.1
    sdSA    <-  0.01
    sf    <-  0.05
    sm    <-  0.05
    qHatSA             <-  qHat_SAAdd(sf=sf, sm=sm)
    pFixSAunlinkedAdd  <-  pFix_SAUnlinkedAdd(x=x, sm=sm, qHat=qHatSA, A=A, Ud=Ud, sd=sd)
    lambda  <-  10


    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

### ROW 1: SLR_pos = 1/2
    SLR_pos  <-  SLR_pos[1]

    # Results for Random Breakpoint Model    
    neutralRBP  <-  pFixNeutralY(Ny=Ny)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvRBP(x=x)
    neutralRBP  <-  neutralRBP / sum(neutralRBP)
    beneficialRBP  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvRBP(x=x)
    beneficialRBP  <-  beneficialRBP / sum(beneficialRBP)
    sexAntagRBP  <-  pFix_SAUnlinkedAdd(x=x, sm=sm, qHat=qHatSA, A=A, Ud=UdSA, sd=sdSA)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvRBP(x=x)
    sexAntagRBP  <-  sexAntagRBP / sum(sexAntagRBP)

    # Results for Exponential Model    
    neutralEXP  <-  pFixNeutralY(Ny=Ny)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvEXP(x=x, lambda=10)
    neutralEXP  <-  neutralEXP / sum(neutralEXP)
    beneficialEXP  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvEXP(x=x, lambda=10)
    beneficialEXP  <-  beneficialEXP / sum(beneficialEXP)
    sexAntagEXP  <-  pFix_SAUnlinkedAdd(x=x, sm=sm, qHat=qHatSA, A=A, Ud=UdSA, sd=sdSA)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvEXP(x=x,lambda=10)
    sexAntagEXP  <-  sexAntagEXP / sum(sexAntagEXP)

## Panel A: Random Breakpoint 
    par(omi=c(0.1, 0.4, 0.1, 0.4), mar = c(2.5,3,3,1.5), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.05*(max(neutralRBP,beneficialRBP,sexAntagRBP))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(neutralRBP ~ x,    lwd=2, lty=1, col=transparentColor('#252525', opacity=1))
        lines(beneficialRBP ~ x, lwd=2, lty=2, col=transparentColor('#252525', opacity=1))
        lines(sexAntagRBP ~ x,   lwd=2, lty=3, col=transparentColor('#252525', opacity=1))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 1.18,  1.275,  expression(paste(SDR[loc], " = 1/2")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03,  1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.175,  expression(paste("Random Breakpoint")), cex=1.3, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.075, expression(paste(italic(f),"(", italic(x),") = 2(1 - ", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,    expression(paste("Density function, ", italic(g), "(", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        

## Panel B: Exponential Model
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.05*(max(neutralEXP,beneficialEXP,sexAntagEXP))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(neutralEXP ~ x,    lwd=2, lty=1, col=transparentColor('#252525', opacity=1))
        lines(beneficialEXP ~ x, lwd=2, lty=2, col=transparentColor('#252525', opacity=1))
        lines(sexAntagEXP ~ x,   lwd=2, lty=3, col=transparentColor('#252525', opacity=1))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(0,0.01,0.02,0.03))
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.175,   expression(paste("Exponential")), cex=1.3, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,   1.075, expression(paste(italic(f),"(", italic(x),") = ", italic(lambda)~italic(e)^-italic(lambda)~italic(x)/(1 - italic(e)^-italic(lambda)))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)


### ROW 2: SLR_pos = 1/10
    SLR_pos  <-  SLR_pos[2]

    # Results for Random Breakpoint Model    
    neutralRBP  <-  pFixNeutralY(Ny=Ny)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvRBP(x=x)
    neutralRBP  <-  neutralRBP / sum(neutralRBP)
    beneficialRBP  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvRBP(x=x)
    beneficialRBP  <-  beneficialRBP / sum(beneficialRBP)
    sexAntagRBP  <-  pFix_SAUnlinkedAdd(x=x, sm=sm, qHat=qHatSA, A=A, Ud=UdSA, sd=sdSA)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvRBP(x=x)
    sexAntagRBP  <-  sexAntagRBP / sum(sexAntagRBP)

    # Results for Exponential Model    
    neutralEXP  <-  pFixNeutralY(Ny=Ny)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvEXP(x=x, lambda=10)
    neutralEXP  <-  neutralEXP / sum(neutralEXP)
    beneficialEXP  <-  pFixBeneficialY(x=x, sI=sb, Ud=Ud, sd=sd)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvEXP(x=x, lambda=10)
    beneficialEXP  <-  beneficialEXP / sum(beneficialEXP)
    sexAntagEXP  <-  pFix_SAUnlinkedAdd(x=x, sm=sm, qHat=qHatSA, A=A, Ud=UdSA, sd=sdSA)*PrCatchSLR(SLR_pos=SLR_pos, x = x)*distNewInvEXP(x=x,lambda=10)
    sexAntagEXP  <-  sexAntagEXP / sum(sexAntagEXP)

## Panel C: Random Breakpoint 
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.05*(max(neutralRBP,beneficialRBP,sexAntagRBP))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(neutralRBP ~ x,    lwd=2, lty=1, col=transparentColor('#252525', opacity=1))
        lines(beneficialRBP ~ x, lwd=2, lty=2, col=transparentColor('#252525', opacity=1))
        lines(sexAntagRBP ~ x,   lwd=2, lty=3, col=transparentColor('#252525', opacity=1))
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 1.18,  1.2,   expression(paste(SDR[loc], " = 1/10")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.03,  1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Density function, ", italic(g), "(", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

## Panel D: Exponential Model
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1.05*(max(neutralEXP,beneficialEXP,sexAntagEXP))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # pInv
        lines(neutralEXP ~ x,    lwd=2, lty=1, col=transparentColor('#252525', opacity=1))
        lines(beneficialEXP ~ x, lwd=2, lty=2, col=transparentColor('#252525', opacity=1))
        lines(sexAntagEXP ~ x,   lwd=2, lty=3, col=transparentColor('#252525', opacity=1))
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.03,  1.075, expression(paste(bold(D))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,  expression(paste("Inversion size (", italic(x), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Neutral")),
                            expression(paste("Beneficial")),
                            expression(paste("Sex Antag."))),
               col     =  c(1,transparentColor('#252525', opacity=1)),
               lty     =  c(1,2,3),
               seg.len =  3,
               lwd     =  2,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
}
