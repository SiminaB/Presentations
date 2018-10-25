##from http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# Multiple plot function - from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##function to replace zero values by minimum non-zero value
replace0 <-  function(x)
{
  xNonZeroMin <- min(x[x > 0])
  x[x==0] <- xNonZeroMin
  x
}

##function based on plotScan function in XCMS package
plotScanEdits <- function(object, scan, mzrange = numeric(),
                          ident = FALSE, 
                          mzlim = NULL, 
                          Intlim = NULL)
{
  if (scan<1 || scan>length(object@scanindex) ) {
    warning("scan out of range")
    return()
  }
  
  ## handle last spectrum
  if (scan == length(object@scanindex)) {
    followingScanIndex <- length(object@env$mz)
  } else {
    followingScanIndex <- object@scanindex[scan+1]
  }
  
  ## hendle empty spectra
  if (object@scanindex[scan] == length(object@env$mz) ||
      object@scanindex[scan] == followingScanIndex) {
    warning("empty scan")
    return()
  }
  
  idx <- (object@scanindex[scan]+1):min(followingScanIndex,
                                        length(object@env$mz), na.rm=TRUE)
  if (length(mzrange) >= 2) {
    mzrange <- range(mzrange)
    idx <- idx[object@env$mz[idx] >= mzrange[1] & object@env$mz[idx] <= mzrange[2]]
  }
  points <- cbind(object@env$mz[idx], object@env$intensity[idx])
  #   title = paste("Median m/z value: ", 
  #                 "Scan time: ", round(object@scantime[scan], 2), 
  #                 sep = "")
  if(!is.null(mzlim) & !is.null(Intlim))
  {
    plot(points, type="h", main = title, 
         ##xlab="m/z", ylab="Intensity",
         xlim=mzlim, ylim=Intlim, xlab="", ylab="")
    mtext("m/z", side=1, line=2)
    mtext("Intensity", side=2, line=2)    
  } else 
  {
    plot(points, type="h", main = title, 
         xlab="m/z", ylab="Intensity")
  }
  
  if (ident)
    return(identify(points, labels = round(points[,1], 1)))
  
  invisible(points)
}

##function to make correlation plots for CK vs creatine ratio
plotCKratioComp <- function(MasterFrame, isoform)
{ 
  ##get all correlations and p-values:
  corr=round(cor(MasterFrame[,isoform], MasterFrame$diffCreat),2)
  pval=signif(cor.test(MasterFrame[,isoform], MasterFrame$diffCreat)$p.value,2)
  
  corr1 = round(cor(MasterFrame[,isoform][as.character(MasterFrame$Category)=="4-7 years"], 
                    MasterFrame$diffCreat[as.character(MasterFrame$Category)=="4-7 years"]),2)
  pval1 = signif(cor.test(MasterFrame[,isoform][as.character(MasterFrame$Category)=="4-7 years"],
                          MasterFrame$diffCreat[as.character(MasterFrame$Category)=="4-7 years"])$p.value,2)
  corr2 = round(cor(MasterFrame[,isoform][as.character(MasterFrame$Category)==">7-11 years"], 
                    MasterFrame$diffCreat[as.character(MasterFrame$Category)==">7-11 years"]),2)
  pval2 = signif(cor.test(MasterFrame[,isoform][as.character(MasterFrame$Category)==">7-11 years"],
                          MasterFrame$diffCreat[as.character(MasterFrame$Category)==">7-11 years"])$p.value,2)
  corr3 = round(cor(MasterFrame[,isoform][as.character(MasterFrame$Category)==">11-18 years"], 
                    MasterFrame$diffCreat[as.character(MasterFrame$Category)==">11-18 years"]),2)
  pval3 = signif(cor.test(MasterFrame[,isoform][as.character(MasterFrame$Category)==">11-18 years"],
                          MasterFrame$diffCreat[as.character(MasterFrame$Category)==">11-18 years"])$p.value,2)
  corr4 = round(cor(MasterFrame[,isoform][as.character(MasterFrame$Category)==">18-29 years"],
                    MasterFrame$diffCreat[as.character(MasterFrame$Category)==">18-29 years"]),2)
  pval4 = signif(cor.test(MasterFrame[,isoform][as.character(MasterFrame$Category)==">18-29 years"],
                          MasterFrame$diffCreat[as.character(MasterFrame$Category)==">18-29 years"])$p.value,2)  
  
  gg1 <- ggplot(MasterFrame, aes_string(x="diffCreat", y=isoform, 
                                        color="Type.x", shape="Type.x")) +
    geom_point() +
    xlab("Creatine/creatinine ratio on the log2 scale") +
    ylab(paste(isoform, " on the log2 scale", sep="")) +
    scale_color_manual(name = "Class", 
                       breaks = c("Control", "DMD"),
                       labels = c("Control", "DMD"),
                       values = c(4,2)) +
    scale_shape_manual(name = "Class", 
                       breaks = c("Control", "DMD"),
                       labels = c("Control", "DMD"),
                       values = c(1,2)) +
    theme(plot.title = element_text(size = 15, hjust = 0.2, vjust=1.5),
          legend.title = element_text(size = 14),
          legend.text = element_text(size=14),
          axis.title = element_text(size=14)) 
  
  gg2 <- ggplot(MasterFrame, aes_string(x="diffCreat", y=isoform, 
                                        color="Category", shape="Category")) +
    geom_point() +
    xlab("Creatine/creatinine ratio on the log2 scale") +
    ylab(paste(isoform, " on the log2 scale", sep="")) +
    scale_color_manual(name = "Age category",
                       values=gg_color_hue(4)) +
    scale_shape_manual(name = "Age category", 
                       values = 1:4) +
    theme(plot.title = element_text(size = 15, hjust = 0.2, vjust=1.5),
          legend.title = element_text(size = 14),
          legend.text = element_text(size=14),
          axis.title = element_text(size=14)) 
  
  title <- bquote( paste( rho, " (overall) = " , .(corr),
                          " (", p, " = ", .(pval)  , ")"))
  title1 <- bquote( paste( rho, " (4-7 years) = " , .(corr1),
                           " (", p, " = ", .(pval1)  , ")"))
  title2 <- bquote( paste( rho, " (>7-11 years) = " , .(corr2),
                           " (", p, " = ", .(pval2)  , ")"))
  title3 <- bquote( paste( rho, " (>11-18 years) = " , .(corr3),
                           " (", p, " = ", .(pval3)  , ")"))
  title4 <- bquote( paste( rho, " (>18-29 years) = " , .(corr4),
                           " (", p, " = ", .(pval4)  , ")"))
  
  ##from http://stackoverflow.com/questions/8615530/place-title-of-multiplot-panel-with-ggplot2
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(6, 1, 
                                             heights = unit(c(rep(0.5,5), 8), "null"))))
  grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(title1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  grid.text(title2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
  grid.text(title3, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
  grid.text(title4, vp = viewport(layout.pos.row = 5, layout.pos.col = 1))
  
  ##print(gg1, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
  print(gg2, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
  
  
  
}