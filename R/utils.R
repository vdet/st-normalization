library(MASS)

### Helper function that converst a color vector with RGB values to RGB+alpha 75 values

transparentRGB <- function(col, alpha=75)
{
    tmp <- c(col2rgb(col), alpha, 255)
    names(tmp) <- c("red", "green", "blue", "alpha", "maxColorValue")
    do.call("rgb", as.list(tmp))
}
transparentRGB <- Vectorize(transparentRGB)


plotDensities <- function(d, col=1:length(as.list(d)),
                           xlim=range(unlist(sapply(as.list(d), function(z) z$x))),
                           ylim=range(unlist(sapply(as.list(d), function(z) z$y))), legend.pos="topright", ...) {
  d <- as.list(d)
  n <- names(d)
  d <- lapply(d, function(z) density(z, na.rm=TRUE))
  par(mar=c(6.1, 6.5, 4.1, 1.1))
  plot(d[[1]], type="n", xlim=xlim, ylim=ylim, cex.lab=2, cex.axis=2, cex.main=2, lwd=2, col=col[1], ...)
  sapply(1:length(d), function(x) lines(d[[x]], lwd=2, col=col[x]))
  if (length(d)>1)
    legend(x=legend.pos, legend=n, bty="n", text.col=col, cex=2)
}


plotScatter <- function(x, y, col=transparentRGB("black", 85), pch=19, log="", smooth=FALSE,
                         main="", xlab="", ylab="", fit="lm", cor=TRUE, cor.x=0.2, cor.y=0.9,
                         id.line=FALSE, id.x=0.4, id.y=0.5, xat=NULL, yat=NULL, file="", ...)
{
    if (file != "") {
        do.call(gsub(".*([a-z]+{3})$", "\\1", file), list(file))
        par(mar=c(6.1, 6.5, 4.1, 1.1))
    }

    if (smooth) {
      smoothScatter(x, y, axes=FALSE, frame=T, log=log,
                    main=main, xlab="", ylab="", cex.main=2.5, ...)
    } else {
      plot(x, y, col=col, pch=pch, cex=2, axes=FALSE, frame=T, log=log,
           main=main, xlab="", ylab="", cex.main=2.5, ...)
    }
    axis(1, at=xat, cex.axis=2.5, padj=0.25)
    axis(2, at=yat, cex.axis=2.5)
    mtext(side=1, text=xlab, line=4, cex=2.5)
    mtext(side=2, text=ylab, line=4, cex=2.5)

    if (fit != "") {
        if (log=="") {
            k <- y
            l <- x
        } else if (log=="x") {
            l <- log10(x)
            k <- y
        } else if (log=="y") {
            k <- log10(y)
            l <- x
        } else if (log=="xy" || log=="yx") {
            k = log10(y)
            l = log10(x)
        } else {
            cat("Error: incorrect 'log' argument\n")
            return()
        }
        if (fit == "lm") {
            abline(lm(k ~ l), col="red")
        } else {
            abline(lqs(k ~ l), col="red")
        } 
    }

    if (id.line) {
        abline(0, 1, col="green", lwd=2)
        pos.x <- min(x, na.rm=TRUE) + id.x*(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
        pos.y <- min(y, na.rm=TRUE) + id.y*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE)) 
        text(pos.x, pos.y, "x=y", col="green", cex=4)
    }
    
    z <- cor.test(x, y, method="spearman")
    rho <- signif(z$estimate, d=2)
    if ((p <- z$p.value)==0) {
        p <- paste("(p<2e-16)", sep="")
    } else {
        p <- paste("(p=", signif(p, d=1), ")", sep="")
    }
    if (cor) {
        pos.x <- min(x, na.rm=TRUE) + cor.x*(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
        pos.y <- min(y, na.rm=TRUE) + cor.y*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE)) 
        p.pos.y <- pos.y - 0.125*(max(y, na.rm=TRUE)-min(y, na.rm=TRUE))
        text(pos.x, pos.y, cex=4, labels=bquote(rho== .(signif(rho, d=2))))
        text(pos.x, p.pos.y, cex=3, labels=p)
    }
    
    if (file != "") {
        dev.off()
    }

    
    
    return(c(rho=rho, p=p))
}

### Helper function that generate a color ramp based on input values and input alpha
colorRampAlpha <- function(..., n, alpha) {
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}

### since original coords were measured on an reference image, this function converts coords pixels
### to account for image WxH of output image
### The "2*" values corrects for the fact that the inputPixels were calculated on a 2times smalle image than the original
convertCoordToPixel <- function(inputPixels, dimensionOriginal, currentImageDimension)
{
  newPixelCoord = 2 * currentImageDimension * inputPixels / dimensionOriginal 
  return(newPixelCoord)
}


legend.gradient <- function (pnts, cols = heat.colors(100), limits = c(0, 1), title = "Legend", ...) 
{
  pnts = try(as.matrix(pnts), silent = T)
  if (!is.matrix(pnts)) 
    stop("you must have a 4x2 matrix")
  if (dim(pnts)[1] != 4 || dim(pnts)[2] != 2) 
    stop("Matrix must have dimensions of 4 rows and 2 columms")
  if (length(cols) < 2) 
    stop("You must have 2 or more colors")
  yvals = seq(min(pnts[, 2]), max(pnts[, 2]), length = length(cols) + 
                1)
  for (i in 1:length(cols)) {
    polygon(x = pnts[, 1], y = c(yvals[i], yvals[i], yvals[i + 
                                                             1], yvals[i + 1]), col = cols[i], border = F)
  }
  text(max(pnts[, 1]), min(pnts[, 2]), labels = limits[1], 
       pos = 4, ...)
  text(max(pnts[, 1]), max(pnts[, 2]), labels = limits[2], 
       pos = 4, ...)
  text(min(pnts[, 1]) , max(pnts[, 2])*1.15, labels = title, adj = c(0, -1), ...)
}


### Function that will plot a vector with values (one for each tissue spot) interpreted as RGB values
### Uses information contained in object as obtained by loading st-slice.Rda
### @imageSize = vector of size two of the form c(WIDTH.HEIGHT) for the output size of the plot
### if NULL, max resolution will be used
plotColorOnSlide <- function(slice, colors, imageSize=NULL, title='', legend=NULL, file='', inputDir = input.dir)
{
  spotVector <- 1:length(slice$spotCoordinates)
  image <- image_read(paste0(inputDir, "/", slice$imagePath))
  imageInfo <- image_info(image)
  
  coordsN <- length(slice$spotCoordinates)
  nColors <- length(unique(colors))
  lineTypes <- rep(1,coordsN)
  lineTypes[is.na(colors)] <- 2
  
  w <- as.numeric(imageInfo$width)
  h <- as.numeric(imageInfo$height)
  oWidth <- slice$dimensionsOriginalImage[1]
  oHeight <- slice$dimensionsOriginalImage[2]
  
  legendXcoords <- c(w*0.03, w*0.1, w*0.1, w*0.03 )
  legendYcoords <- c(h*0.03, h*0.35, h*0.35, h*0.03 )
  
  drawImage <- image_draw(image)

  radius <- slice$pixelSizeFactor*w
  
  for(i in spotVector){
    draw.circle(convertCoordToPixel(slice$spotCoordinates[[i]][3], oWidth, w),
                convertCoordToPixel(slice$spotCoordinates[[i]][4], oHeight, h) ,
                radius = radius, border = 'black', col = colors[i], nv=20, lwd = 10, lty = lineTypes[i])
  }

  if (!is.null(legend)) {
    legendXcoords = c(w*0.03, w*0.1, w*0.1, w*0.03 )
    legendYcoords = c(h*0.03, h*0.35, h*0.35, h*0.03 )
    legend.gradient(cbind(legendXcoords, legendYcoords), 
                    cols = colorRampAlpha(legend$color, n=nColors, alpha=1), 
                    title = legend$title,
                    limits = signif(legend$limits, d=2),
                    cex = legend$cex)
  }
  
  text(w*0.5, h*0.1, labels = title, cex = 20)
  
  if(file != ""){
    if(is.null(imageSize)){
      image_write(drawImage, path = file, quality = 100, density = 600)
      dev.off()
    }else{
      finalImage = image_capture()
      dev.off()
      image_write(image_resize(finalImage, geometry = paste(imageSize[1],imageSize[2],sep="x")), path = file)
    }

  }else{
    finalImage = image_capture()
    dev.off()
    plot(finalImage)
  }
}

### Given a vector (scalar) with numeric values for each tissue spot, this function will generate a color ramp
### and plot in on the ST tissue slide using the plotColorOnSlide function
plotScalarOnSlide <- function(slice, scalar,
                              colorVectorForScale = c('lightblue', 'greenyellow', 'yellow', 'orange','red', 'purple'),
                              imageSize=NULL, title='', legendTitle='', file='')
{
  nColors <- ifelse(max(scalar, na.rm = T)<1000, 1000, max(scalar))
  reScaledData <- rescale(scalar, c(1,nColors))
  colors <- colorRampAlpha(colorVectorForScale, n=nColors, alpha=0.75)[reScaledData]

  legend <- list(colors = colorVectorForScale, title=legendTitle, limits=range(scalar, na.rm=T), cex=10)
  plotColorOnSlide(slice=slice, colors=colors, imageSize=imageSize, title=title, legend=legend, file=file)
}


regionBoxplot <- function(v, slice, yat=NULL, ylab="", main="", file="")
{
    if (file != "") {       
        do.call(gsub(".*([a-z]+{3})$", "\\1", file), list(file))
        par(mar=c(6.1, 6.5, 4.1, 1.1))
    }

    r <- rep(NA, length(v))
    names(r) <- names(v)
    for (z in names(slice$regions)) {
      r[slice$regions[[z]]] <- z
    }
    
    i <- is.na(r) | is.na(v)
    v <- v[!i]
    r <- r[!i]
    col <- transparentRGB(slice$regionColors[order(names(slice$regionColors))], 85)

    boxplot(v ~ r, outline=FALSE, xlab="", ylab="", axes=FALSE, frame=TRUE, lwd=1)
    stripchart(v ~ r, add=T, method="jitter", pch=19, vertical=T, col=col, cex=2)

    axis(1, at=1:5, labels=c('Cell.\nFibrosis', 'Epith.', 'Pure\nfibrosis', 'Imm.', 'Mix\nepi.+fib.'), cex.axis=1.5, padj=0.75)
    axis(2, at=yat, cex.axis=2.5)
    mtext(side=2, text=ylab, line=4, cex=2.5)

    if (file != "") {
        dev.off()
    }
}


plotGenes <- function(slice, genes, norm='dca', fileBase=genes[1], plotRegions=TRUE) {
  g <- intersect(genes, rownames(slice[[norm]]))
  if (length(g)==0) {
    cat("No gene found\n")
    return()
  } else if (length(g)==1) {
    x <- slice[[norm]][g,]
  } else {
      x <- apply(slice[[norm]][g,], 2, function(z) log10(sum(10**z)))
  }
      
 plotScalarOnSlide(slice, x, title=fileBase, legendTitle='Expression',
                   file=paste0(figures.dir, '/slice-', fileBase, '-', norm, '.jpg'))

 if (plotRegions)
    regionBoxplot(x, slice, ylab="Expression", main=fileBase,
                  file=paste0(figures.dir, "/boxplot-region-", fileBase, "-", norm, ".pdf"))
}


