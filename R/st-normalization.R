

# attemps to generate root Folder from working directory, supposing you downloaded the three main folders as given.
# will not work if WD =/= folder with this script. In that case, replace base.dir with root folder containing the three
# main folders (i.e. input, output and R)

# base.dir <- gsub("/R","",getwd())
# setwd("../")

base.dir <- '.'
code.dir <- paste0(base.dir, "/R")
input.dir <- paste0(base.dir, "/input")
output.dir <- paste0(base.dir, "/output")
figures.dir <- output.dir #paste0(output.dir, "/figures")


# Loads necessary packages
# This piece of code will try to load the packages names in the pkgVector object
# If necessary, it tries downloading and installing the packages using either install.packages or biocLite from
# bioconductor

pkgVector= c("magick", "plotrix", "MASS")

if(!all(lapply(pkgVector, require, character.only= T ))){
  notLoaded = which(lapply(pkgVector, require, character.only= T ) == FALSE)
  install.packages(pkgVector[notLoaded], repos = "http://cran.us.r-project.org");
  source("https://bioconductor.org/biocLite.R")
  notLoaded = which(lapply(pkgVector, require, character.only= T ) == FALSE)
  biocLite(pkgVector[notLoaded], ask = F) # will not prompt ask to update packages
  lapply(pkgVector, require, character.only= T )
} else{
  lapply(pkgVector, require, character.only= T )
}

source(paste0(code.dir, "/utils.R"))


##Load ST data object, list that contains all necessary data
load(paste0(input.dir, '/st-slice.Rda'))

##Plot regions
plotColorOnSlide(stSlice, stSlice$regionColorsVector, title='Regions', file=paste0(figures.dir, '/slice-regions.jpg'))


##Cell numbers
cellNumbers <- cbind(Total=apply(stSlice$segmentation[, c('cancer_cells', 'fibroblasts', 'other_cells')], 1, sum),
                     Epithelial=stSlice$segmentation[, 'cancer_cells'],
                     Fibroblasts=stSlice$segmentation[, 'fibroblasts'],
                     Others=stSlice$segmentation[, 'other_cells'])
summary(cellNumbers)
apply(cellNumbers, 2, sum)



##Total count
totalCounts <- 2**stSlice$rawCounts-1
totalCounts[totalCounts<1] <- 0
totalCounts <- colSums(totalCounts)
summary(totalCounts)

plotScalarOnSlide(stSlice, log2(totalCounts), title='Total Counts', legendTitle='Log2(counts)',
                  file=paste0(figures.dir, '/slice-log2-total-counts.jpg'))

plotScatter(cellNumbers[, 'Epithelial'], log2(totalCounts), xlab='# epithelial cells', ylab='# reads [log2]', fit='lqs',
            col=transparentRGB(stSlice$regionColorsVector, 85),
            cor.x=0.7, cor.y=0.2, file=paste0(figures.dir, '/scatter-number-epithelial-cells-vs-number-reads.pdf'))

plotScatter(cellNumbers[, 'Total'], log2(totalCounts), xlab='# cells', ylab='# reads [log2]', fit='lqs',
            col=transparentRGB(stSlice$regionColorsVector, 85),
            cor.x=0.7, cor.y=0.2, file=paste0(figures.dir, '/scatter-number-cells-vs-number-reads.pdf'))


fit <- lm(log2(totalCounts) ~ cellNumbers[, 'Epithelial'] + cellNumbers[, 'Fibroblasts'] + cellNumbers[, 'Others'])
summary(fit)
summary(fit)$coef


##Counts normalization
x <- (2**stSlice$rawCounts)-1
stSlice$adjCounts <- t(apply(x, 1, function(z) {
  z/totalCounts}))
stSlice$adjCounts <- log2(stSlice$adjCounts+1)
rownames(stSlice$adjCounts) <- rownames(stSlice$rawCounts)



##TG and VIM
plotGenes(stSlice, 'TG', norm='rawCounts')
plotGenes(stSlice, 'TG', norm='adjCounts')
plotGenes(stSlice, 'TG', norm='dca')

plotGenes(stSlice, 'VIM', norm='rawCounts')
plotGenes(stSlice, 'VIM', norm='adjCounts')
plotGenes(stSlice, 'VIM', norm='dca')


plotScatter(cellNumbers[,'Epithelial'], cellNumbers[,'Fibroblasts'],
            xlab='# epithelial cells', cor.x=0.75,
            ylab='# fibroblasts', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-number-epithelial-vs-number-fibroblasts.pdf'))

plotScatter(100*cellNumbers[,'Epithelial']/cellNumbers[, 'Total'], 100*cellNumbers[,'Fibroblasts']/cellNumbers[, 'Total'],
            xlab='% epithelial cells', cor.x=0.75,
            ylab='% fibroblasts', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-percent-epithelial-vs-percent-fibroblasts.pdf'))


plotScatter(stSlice$rawCounts['TG',], stSlice$rawCounts['VIM',],
            xlab='TG expression [log2]',
            ylab='VIM expression [log2]', fit='lqs',  col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-TG-vs-VIM-rawCounts.pdf'))
plotScatter(stSlice$adjCounts['TG',], stSlice$adjCounts['VIM',],
            xlab='TG expression [log2]', cor.x=0.75,
            ylab='VIM expression [log2]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-TG-vs-VIM-adjCounts.pdf'))
plotScatter(stSlice$dca['TG',], stSlice$dca['VIM',],
            xlab='TG expression [log2]', cor.x=0.75,
            ylab='VIM expression [log2]', fit='lqs',  col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-TG-vs-VIM-dca.pdf'))


##Correlations
g <- sample(1:nrow(stSlice$rawCounts), 500)
c.raw <- sapply(g, function(i)
  sapply(g, function(j)
        if (i==j) {NA} else {cor(stSlice$rawCounts[i,], stSlice$rawCounts[j,], method='spearman')}))

c.adj <- sapply(g, function(i)
  sapply(g, function(j)
        if (i==j) {NA} else {cor(stSlice$adjCounts[i,], stSlice$adjCounts[j,], method='spearman')}))

c.dca <- sapply(g, function(i)
  sapply(g, function(j)
        if (i==j) {NA} else {cor(stSlice$dca[i,], stSlice$dca[j,], method='spearman')}))

d <- list(c.raw, c.adj, c.dca)
names(d) <- c('Raw counts', 'Adj. counts', 'DCA')
pdf(paste0(figures.dir, '/density-gene-gene-correlation.pdf'))
par(mar=c(6.1, 6.5, 4.1, 1.1))
plotDensities(d, xlab='Gene-gene correlation', main='', legend.pos='topleft')
dev.off()


##Scatter plot for Supp. Fig. 1
q <- quantile(stSlice$adjCounts['TG',], 0.95)
stSlice$adjCounts['TG', stSlice$adjCounts['TG',]>q] <- q

plotGenes(stSlice, 'TG', norm='adjCounts')

plotScatter(stSlice$rawCounts['TG',], log2(totalCounts),
            xlab='TG raw counts [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Total counts [log2]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-TG-raw-counts-vs-total-counts.pdf'))
plotScatter(stSlice$adjCounts['TG',], log2(totalCounts),
            xlab='TG norm. counts [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Total counts [log2]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-TG-normalized-counts-vs-total-counts.pdf'))
plotScatter(stSlice$dca['TG',], log2(totalCounts),
            xlab='TG DCA [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Total counts [log2]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-TG-dca-vs-total-counts.pdf'))

plotScatter(stSlice$rawCounts['TG',], log10(cellNumbers[,'Epithelial']+1),
            xlab='TG raw counts [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Epith. cell per spot [log10]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-TG-raw-counts-vs-epith-cell-density.pdf'))
plotScatter(stSlice$adjCounts['TG',], log10(cellNumbers[,'Epithelial']+1),
            xlab='TG norm. counts [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Epith. cell per spot [log10]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-TG-normalized-counts-vs-epith-cell-density.pdf'))
plotScatter(stSlice$dca['TG',], log10(cellNumbers[,'Epithelial']+1),
            xlab='TG DCA [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Epith. cell per spot [log10]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-TG-dca-vs-epith-cell-density.pdf'))

q <- quantile(stSlice$adjCounts['VIM',], 0.95)
stSlice$adjCounts['VIM', stSlice$adjCounts['VIM',]>q] <- q

plotScatter(stSlice$rawCounts['VIM',], log2(totalCounts),
            xlab='VIM raw counts [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Total counts [log2]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-VIM-raw-counts-vs-total-counts.pdf'))
plotScatter(stSlice$adjCounts['VIM',], log2(totalCounts),
            xlab='VIM norm. counts [log2]', cor.x=0.70, cor.y=0.25,
            ylab='Total counts [log2]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-VIM-normalized-counts-vs-total-counts.pdf'))
plotScatter(stSlice$dca['VIM',], log2(totalCounts),
            xlab='VIM DCA [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Total counts [log2]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-VIM-dca-vs-total-counts.pdf'))

plotScatter(stSlice$rawCounts['VIM',], log10(cellNumbers[,'Epithelial']+1),
            xlab='VIM raw counts [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Epith. cell per spot [log10]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-VIM-raw-counts-vs-epith-cell-density.pdf'))
plotScatter(stSlice$adjCounts['VIM',], log10(cellNumbers[,'Epithelial']+1),
            xlab='VIM norm. counts [log2]', cor.x=0.70, cor.y=0.25,
            ylab='Epith. cell per spot [log10]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-VIM-normalized-counts-vs-epith-cell-density.pdf'))
plotScatter(stSlice$dca['VIM',], log10(cellNumbers[,'Epithelial']+1),
            xlab='VIM DCA [log2]', cor.x=0.75, cor.y=0.25,
            ylab='Epith. cell per spot [log10]', fit='lqs', col=transparentRGB(stSlice$regionColorsVector, 85),
            file=paste0(figures.dir, '/scatter-VIM-dca-vs-epith-cell-density.pdf'))


##Quit
sessionInfo()
q(save='no')
