## This file contains function code for making figures, it's
## intended to be sourced from the main file.


print.letter <- function(label="(a)",xy=c(0.1,0.925),...) {
    tmp <- par("usr")
    text.x <- tmp[1]+xy[1]*diff(tmp[1:2])   #x position, diff=difference
    text.y <- tmp[3]+xy[2]*diff(tmp[3:4])   #y position
    text(x=text.x, y=text.y, labels=label, ...)
}

make.file <- function(type=c("png","pdf", "none"), filename,
                      width, height, res){
    ## Pass it file type and dimensions in inches. It creates file.
    type <- match.arg(type)
    ## If no extension given, add one
    if(length(grep(type, filename))==0)
        filename <- paste0(filename,".",type)
    if(type=="png") png(filename, width=width,
       height=height, units="in", res=res)
    else if(type=="pdf"){pdf(filename, width=width, height=height)}
    else if(dev.cur()==1) dev.new(width=width, height=height)
}


figure1 <- function(type="png", res=500, width, height){
    on.exit(if(type!="none") dev.off())
    make.file(type=type, paste("Figures/Figure1", res,
              sep="_"), width=width, height=height, res=res)
    ## Plot the properties of the model for the two life
    ## histories
    par(mfcol=c(2,1), mar=c(2.8,2.8,.75,.75), tck=-0.02, oma=c(0,0,0,0),
        col.axis = gray(.3), mgp=c(1.25, .25, 0), cex.lab = 1.1,
        cex=.8, cex.axis = .8, col.axis=label.col)
    plot(0,0, type="n", main="", ylim=c(0, 1),
         xlim=c(0, 1), axes = FALSE,
         xlab="Spawning Biomass (% Equilibrium)",
         ylab="Recruits (% Equilibrium)")
    print.letter(label = "(a)", xy=c(.05, .95))
    for(SpeciesI in 1:3){
        xx <- Results.examples$lifehist.curves[[SpeciesI]]
        lwd.tmp <- 2
        tempx <- (sum(xx$SBPR_a)*xx$SR_alpha-1)/xx$SR_beta
        tempy <- tempx*xx$SR_alpha/(1+xx$SR_beta*tempx)
        lines(x=xx$SB/tempx, y=(xx$SR_alpha * xx$SB / (1 +
              xx$SR_beta*xx$SB))/tempy,
              lty=species.ltys[SpeciesI],
              col=species.cols[SpeciesI], lwd=lwd.tmp)
    }
    lines(x=c(0,1), y=c(0,1), col=gray(.5))
    axis(1, at=c(0,.5,1), col=border.col, mgp=par()$mgp)
    axis(2, at=c(0,.5,1),col=border.col, mgp=par()$mgp)
    legend("bottomright", legend=SpeciesSet_Example, col=species.cols,
           lty=species.ltys, lwd=lwd.tmp, bty="n", ncol = 1)
    box(col = border.col)
    ## Old plot of survival by age
    ## plot( x=xx$AgeSet, y=(xx$Surv_a), type="l", main="", pch=xx.pch,
    ##      col=xx.col, lty=xx.lty, lwd=lwd.tmp, ylim=c(0, .0004), xlim=c(0,10),
    ##      xlab="Age", ylab="Survival")
    ## lines( x=yy$AgeSet, y=(yy$Surv_a), type="l",  main="",
    ##      col=yy.col, lty=yy.lty, lwd=lwd.tmp, pch=yy.pch)
    par( col.axis=label.col)
    plot( 0,0, type='n', ylim=c(0, .2), xlim=c(0,100),
         ylab="Biomass per Recruit (Normalized)", xlab="Age", axes=FALSE)
    for(SpeciesI in 1:3){
        xx <- Results.examples$lifehist.curves[[SpeciesI]]
        lines( x=xx$AgeSet, y=xx$BPR_a/sum(xx$BPR_a),  col=species.cols[SpeciesI],
              lty=species.ltys[SpeciesI], lwd=lwd.tmp)
    }
    print.letter(label = "(b)", xy=c(.05, .95))
    axis(1, col=border.col, mgp=par()$mgp)
    axis(2, at=c(0, 1), col=border.col, mgp=par()$mgp)
    box(col = border.col)
}



## Figure 2,3 and 4 are the yield curves, different
figure2master <- function(MSY_type, SpeciesI, type,
                          res, width, height, name){
    on.exit(if(type!="none") dev.off())
    make.file(type=type, paste0("Figures/Figure",
              name, "_", res), width=width, height=height, res=res)
    ## crazy function to map the multipliesr to colors; for plotting only
    Multiplier.cols <- function(x) 0+
        (.7-0)*(x-min(Multiplier.values))/(max(Multiplier.values)-min(Multiplier.values))
    df <- droplevels(subset(Results.examples$yield.curves, Species==SpeciesSet_Example[SpeciesI] &
                            !Delta %in% c("Diff+", "Diff-")))
    df <- within(df, {Cmsy <- Cmsy/max(Catch); Catch <- Catch/max(Catch) })
    par(mar=.07*c(1,1,1,1), tck=-0.015, oma=c(2,2,1.25,.3),
        mgp=c(.5, .05,0), mfcol=c(3,4), cex=.8, cex.axis=.6,
        col.axis=label.col)
    xlim <- c(0,.4);
    ylim <- c(0, 1.1)
    lwd.tmp <- 1
    for(ParI in 1:length(ParOrder)){
        ## Make big plotting area
        plot(1, type="n", axes=FALSE, xlim=xlim, ylim=ylim,xlab="", ylab="")
        mtext(ParName[ParI], line=-1.1, cex=.6)
        ## Draw column headers
        text.temp <- rep(NA, len=length(ParOrder))
        text.temp[c(1,4,7,10)] <- c("Mortality", "Recruitment", "Size", "Maturity")
        mtext(side=3, outer=FALSE, line=.25, text=text.temp[ParI], cex=.8)
        ## Draw yield curves for different multipliers, and add MSY points
        for(mult in levels(df$Delta)){
            temp <- droplevels(subset(df, Delta==mult & F_Type=="Optimum" &
                                      Parameter==ParSet_Full[ParI] &
                                      Species==SpeciesSet_Example[SpeciesI]))
            col.temp <- gray(Multiplier.cols(temp$Multiplier[1]))
            ##      print(col.temp)
            ##  left off here, need to save more than just FMSY
            ## then call it here
            with(temp, lines(x=F, y=Catch, col=col.temp, lwd=1,
                             lty=ifelse(mult=="Default", 1,3)))
            with(temp, points(x=Fmsy, y=Cmsy, pch=16, col=col.temp))
        }
        youtside <- ifelse(ParI %in% c(1,2,3), TRUE, FALSE)
        xoutside <- ifelse(ParI %in% c(3,6,9,11), TRUE, FALSE)
        axis(2,  labels=youtside, col=border.col,
             mgp=par()$mgp, tck=ifelse(youtside, par()$tck, .015))
        axis(1,  labels=xoutside, col=border.col,
             mgp=par()$mgp, tck=ifelse(xoutside, par()$tck, .015))
        box(col=border.col)
        ## Draw axis labels
        mtext(side=1, outer=TRUE, line=.8, text="Fishing Effort (F)", cex=.7)
        mtext(side=2, outer=TRUE, line=.75, text="Normalized Catch", cex=.7)
    }
    legend("topright", legend= paste0("exp(", c(-.2, -.1, 0, .1, .2), ")"),
           col=gray(Multiplier.cols(as.numeric(Multiplier.values))),
           lty=c(3,3,1,3,3),  lwd=lwd.tmp, cex=.6, bty="n",
           pch=16)
}
figure2a <- function(MSY_type="Fmsy", type="png", res=500, width=width,
                    height=height){
    figure2master(MSY_type=MSY_type, SpeciesI=1, type=type, res=res, width=width,
                    height=height, name="2a")
}
figure2b <- function(MSY_type="Fmsy", type="none", res=500, width=width,
                    height=height){
    figure2master(MSY_type=MSY_type, SpeciesI=2, type=type, res=res, width=width,
                  height=height, name="2b")
}
figure2c <- function(MSY_type="Fmsy", type="none", res=500, width=width,
                    height=height){
    figure2master(MSY_type=MSY_type, SpeciesI=3, type=type, res=res, width=width,
                  height=height, name="2c")
}



## FIgure  is the  multiplier by MSY
figure3 <- function(type="png", res=500, width, height){
    on.exit(if(type!="none") dev.off())
    make.file(type=type, paste("Figures/Figure3", res,
              sep="_"), width=width, height=height, res=res)
    par(mar=.07*c(1,1,1,1), tck=-0.015, oma=c(2,2,1.25,.3),
        mgp=c(.5, .05,0), mfcol=c(3,4), cex=.8, cex.axis=.6,
        col.axis=label.col)
    xlim <- c(.78,1.27)
    F_I <- 1
    lwd.tmp <- 2
    ylim <- range(PlotResults[F_I,,,,'New','Cmsy'] /
                      PlotResults[F_I,,,,'Orig','Cmsy'], na.rm=TRUE )
    for(ParI in 1:length(ParOrder)){
        ## Make big plotting area
        plot(1, type="n", axes=FALSE, xlim=xlim, ylim=ylim)
        if(ParI==1)
            legend("right", legend=SpeciesSet_Example,
                   col=species.cols, cex=.6,
                   lty=species.ltys, lwd=lwd.tmp, bty="n", ncol = 1)
        mtext(ParName[ParI], line=-1.1, cex=.6)
        ## Draw column headers
        text.temp <- rep(NA, len=length(ParOrder))
        text.temp[c(1,4,7,10)] <- c("Mortality", "Recruitment", "Size", "Maturity")
        mtext(side=3, outer=FALSE, line=.25, text=text.temp[ParI], cex=.8)
        ## Draw elasticity curve
        for(SpeciesI in 1:length(SpeciesSet_Example)){
            ## plot(1, type="n", xlab="Multplier", ylab=, ylim=Ylim2, xlim=exp(c(-0.2,0.2))) # c(0.3,3)
            ## This was missing the Delta values so changed
            ## them on 2/6
            ##X = exp( c( -0.2, -0.1, 0, 0.1, 0.2 ) )
            X <- Multiplier.values
            Y <- PlotResults[F_I,SpeciesI,ParOrder[ParI],,'New','Cmsy'] /
                      PlotResults[F_I,SpeciesI,ParOrder[ParI],,'Orig','Cmsy']
            lines(x=X, y=Y, col=species.cols[SpeciesI],
                  lty=species.ltys[SpeciesI], lwd=lwd.tmp) # c(0.3,3)
            abline(v=1, h=1, lty="dashed", col=border.col, lwd=.5)
            youtside <- ifelse(ParI %in% c(1,2,3), TRUE, FALSE)
            xoutside <- ifelse(ParI %in% c(3,6,9,11), TRUE, FALSE)
            axis(2, at=c(0:5), labels=youtside,
                 col=border.col, mgp=par()$mgp, tck=ifelse(youtside, par()$tck, .015))
            axis(1, at=c(.8,1, 1.2), labels=xoutside,
                 col=border.col, mgp=par()$mgp, tck=ifelse(xoutside, par()$tck, .015))
            box(col=border.col)
        }
        ## Draw axis labels
        mtext(side=1, outer=TRUE, line=.75, text="Multiplier", cex=.7)
        mtext(side=2, outer=TRUE, line=.75, text=expression(MSY[Delta] / MSY[0]), cex=.7)
    }
}


## Figure 4 the ratios of SB and F MSY
figure4 <- function(type="png", res=500, width, height){
    on.exit(if(type!="none") dev.off())
    make.file(type=type, paste("Figures/Figure4", res,
              sep="_"), width=width, height=height, res=res)
    par(mar=.07*c(1,1,1,1), tck=-0.015, oma=c(2,2,1.25,.3),
        mgp=c(.5, .05,0), mfcol=c(3,4), cex=.8, cex.axis=.6,
        col.axis=label.col)
    lwd.tmp <- 1.5
    F_I=1
    xlim <- range(PlotResults[F_I,,,,'New','SBmsy'] /
                  PlotResults[F_I,,,,'Orig','SBmsy'], na.rm=TRUE)
    ylim <- range(PlotResults[F_I,,,,'New','Fmsy'] /
                  PlotResults[F_I,,,,'Orig','Fmsy'], na.rm=TRUE)
    for(ParI in 1:length(ParOrder)){
        ## Make big plotting area
        plot(1, type="n", axes=FALSE, xlim=xlim, ylim=ylim ,xlab="", ylab="") #
        if(ParI==1)
            legend("bottomright", legend=SpeciesSet_Example, lty=species.ltys,
                   col=species.cols, bty="n", cex=.6,
                   pch=species.pchs, lwd=lwd.tmp)
        mtext(ParName[ParI], line=-1.1, cex=.6)
        ## Draw column headers
        text.temp <- rep(NA, len=length(ParOrder))
        text.temp[c(1,4,7,10)] <- c("Mortality", "Recruitment", "Size", "Maturity")
        mtext(side=3, outer=FALSE, line=.25, text=text.temp[ParI], cex=.8)
        F_I=1
        ## Draw elasticity curve
        for(SpeciesI in 1:length(SpeciesSet_Example)){
            X1 <- PlotResults[F_I,SpeciesSet_Example[SpeciesI],ParOrder[ParI],,'New','SBmsy'] /
                PlotResults[F_I,SpeciesSet_Example[SpeciesI],ParOrder[ParI],,'Orig','SBmsy']
            Y1 <- PlotResults[F_I,SpeciesSet_Example[SpeciesI],ParOrder[ParI],,'New','Fmsy'] /
                PlotResults[F_I,SpeciesSet_Example[SpeciesI],ParOrder[ParI],,'Orig','Fmsy']
            X <- na.omit(cbind(X1,Y1))[,1];
            Y <- na.omit(cbind(X1,Y1))[,2]
            points(x=X[-length(X)], y=Y[-length(Y)], col=species.cols[SpeciesI],
                   pch=species.pchs[SpeciesI], cex=.65)
            shape::Arrows(x0=X[length(X)-1], x1=X[length(X)], y0=Y[length(Y)-1],
                          y1=Y[length(Y)], col=species.cols[SpeciesI],
                          cex=.1, arr.type="curved", arr.length=.2,
                          lty=1, segment=FALSE)
            segments(x0=X[-length(X)], x1=X[-1], y0=Y[-length(Y)], y1=Y[-1],
                     col=species.cols[SpeciesI],  lwd=lwd.tmp,
                     lty=species.ltys[SpeciesI])
        }
        abline(v=1, h=1, lty="dashed", col=border.col, lwd=.5)
        youtside <- ifelse(ParI %in% c(1,2,3), TRUE, FALSE)
        xoutside <- ifelse(ParI %in% c(3,6,9,11), TRUE, FALSE)
        axis(2, labels=youtside,
             col=border.col, mgp=par()$mgp, tck=ifelse(youtside, par()$tck, .015))
        axis(1, labels=xoutside,
             col=border.col, mgp=par()$mgp, tck=ifelse(xoutside, par()$tck, .015))
        box(col=border.col)
    }
    ## Draw axis labels
    mtext(side=1, outer=TRUE, line=1, text=expression(SB[MSY[Delta]] /
                                      SB[MSY[0]]), cex=.7)
    mtext(side=2, outer=TRUE, line=.5, text=expression(F[MSY[Delta]] /
                                        F[MSY[0]]), cex=.7)
}

# NEW FIGURE SUMMARIZING ELASTICITIES
figure_new <- function(type="png", res=500, width, height){
    on.exit(if(type!="none") dev.off())
    make.file(type=type, paste("Figures/FigureNew", res,
              sep="_"), width=width, height=height, res=res)
    ## Plot the properties of the model for the two life
    ## histories
    par(mar=.07*c(1,1,1,1), tck=-0.015, oma=c(2,2,1.25,.3),
        mgp=c(.5, .05,0), mfcol=c(3,4), cex=.8, cex.axis=.6,
        col.axis=label.col)
    lwd.tmp <- 1.5
    for(ParI in 1:length(ParOrder)){
        ## Make big plotting area
        plot(1, type="n", axes=FALSE, xlim=c(0.5,3.5), ylim=c(-12,8) ,xlab="", ylab="") #
        if(ParI==1)
            legend("bottomleft", legend=SpeciesSet_Example, lty=species.ltys,
                   col=species.cols, bty="n", cex=.6,
                   pch=species.pchs, lwd=lwd.tmp)
        mtext(ParName[ParI], line=-1.1, cex=.6)
        ## Draw column headers
        text.temp <- rep(NA, len=length(ParOrder))
        text.temp[c(1,4,7,10)] <- c("Mortality", "Recruitment", "Size", "Maturity")
        mtext(side=3, outer=FALSE, line=.25, text=text.temp[ParI], cex=.8)
        for(SpeciesI in 1:3){
            lines( x=1:3, y=Elas.table[c(0,3,6)+SpeciesI,ParI], col=species.cols[SpeciesI],
              lwd=lwd.tmp, lty=species.ltys[SpeciesI])
        }
        abline(h=0, lty="dashed", col=border.col, lwd=.5)
        youtside <- ifelse(ParI %in% c(1,2,3), TRUE, FALSE)
        xoutside <- ifelse(ParI %in% c(3,6,9,11), TRUE, FALSE)
        xlabels <- list(FALSE, c(expression(F[MSY]),expression(SB[MSY]),expression(MSY)))[[xoutside+1]]
        axis(2, labels=youtside,
             col=border.col, mgp=par()$mgp, tck=ifelse(youtside, par()$tck, .015))
        axis(1, at=1:3, labels=xlabels,
             col=border.col, mgp=par()$mgp, tck=ifelse(xoutside, par()$tck, .015))
        box(col=border.col)
    }
    ## Draw axis labels
    mtext(side=1, outer=TRUE, line=1, text="Management target")
    mtext(side=2, outer=TRUE, line=.5, text="Elasticities")
}


## I edited this function to remove plotting the z-scale
my.filled.contour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1,
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
    col = color.palette(length(levels) - 1), plot.title, plot.axes,
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
    axes = TRUE, frame.plot = axes, ...)
{
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
    ## mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    ## on.exit(par(par.orig))
    ## w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    ## layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    ## par(las = las)
    ## mar <- mar.orig
    ## mar[4L] <- mar[2L]
    ## mar[2L] <- 1
    ## par(mar = mar)
    ## plot.new()
    ## plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i",
    ##     yaxs = "i")
    ## rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    ## if (missing(key.axes)) {
    ##     if (axes)
    ##         axis(4)
    ## }
    ## else key.axes
    ## box()
    ## if (!missing(key.title))
    ##     key.title
    ## mar <- mar.orig
    ## mar[4L] <- 1
    ## par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    .filled.contour(x, y, z, levels, col)
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

## ## These are old plots from trying the contours
## ## Figure 6, the contour plots
## ## Figure 2 and 3 are the yield curves, different
## figure6master <- function(MSY_type, type, res, width, height){
##     on.exit(if(type!="none") dev.off())
##     make.file(type=type, paste0("Figures/Figure6", "_", res), width=width, height=height, res=res)
##     par(mar=.07*c(1,1,1,1), tck=-0.015, oma=c(2,2,1.25,.3),
##         mgp=c(.5, .05,0), mfcol=c(3,4), cex=.8, cex.axis=.6,
##         col.axis=label.col)
##     scaley <- 0*.01
##     scalex <- 0*.1
##     ylim <- c(0-scaley*(1-0), 1+scaley*(1-0))
##     xlim <- c(-1-scalex*(2+1), 2+scalex*(2+1))
##     zlim <- c(0,20)
##     range(sqrt(abs(Elas.grid.long$Elasticity)), na.rm=TRUE)
##     z.palette <- function (n, gamma = 2.2)
##         gray(seq.int(from = 1^gamma, to = 0^gamma, length.out = n)^(1/gamma))# gray.colors
##     for(ParI in 1:length(ParOrder)){
##         ## Make big plotting area
##         df.temp <- droplevels(subset(Elas.grid.long, Parameter==ParOrder[ParI] & Metric==MSY_type))
##         df.temp <- dcast(data=df.temp, formula=LMARR~M, value.var='Elasticity')
##         colnames(df.temp) <- row.names(df.temp) <- NULL
##         df.temp <- as.matrix(df.temp)[,-1]
##         df.temp[df.temp>20] <- 20
##         df.temp[df.temp< -20] <- -20
##         my.filled.contour(x=LMARR.seq, y=M.seq, z=(abs(df.temp)), axes=FALSE, zlim=zlim,
##                           color.palette=z.palette, nlevels=100, xaxs="r",
##                           xlim=xlim, ylim=ylim)
##         mtext(ParName[ParI], line=-1.1, cex=.6)
##         ## Draw column headers
##         text.temp <- rep(NA, len=length(ParOrder))
##         text.temp[c(1,4,7,10)] <- c("Mortality", "Recruitment", "Size", "Maturity")
##         mtext(side=3, outer=FALSE, line=.25, text=text.temp[ParI], cex=.8)
##         youtside <- ifelse(ParI %in% c(1,2,3), TRUE, FALSE)
##         xoutside <- ifelse(ParI %in% c(3,6,9,11), TRUE, FALSE)
##         axis(2, labels=youtside,  col=border.col, mgp=par()$mgp, tck=ifelse(youtside, par()$tck, .015))
##         axis(1, labels=xoutside, col=border.col, mgp=par()$mgp, tck=ifelse(xoutside, par()$tck, .015))

##         ## if(ParI %in% c(1,2,3)){
##         ##     axis(2, col=border.col,
##         ##          mgp=par()$mgp, tck=ifelse(youtside, par()$tck, .015))
##         ## } else if (ParI %in% c(3,6,9,11)) {
##         ##     axis(1,  col=border.col,
##         ##          mgp=par()$mgp, tck=ifelse(xoutside, par()$tck, .015))
##         ## }
##         box(col=border.col)
##         points(x=LMARR_Set_Example, y=M_Set_Example)
##     }
##     ## Draw axis labels
##     mtext(side=1, outer=TRUE, line=.8, text="LMARR", cex=.7)
##     mtext(side=2, outer=TRUE, line=.75, text="M", cex=.7)
## }
## figure6 <- function(type="none", res=500, width=width,
##                     height=height){
##     figure6master(MSY_type="Cmsy", type=type, res=res, width=width,
##                     height=height)
## }

