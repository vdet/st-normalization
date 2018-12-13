

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
plotColorOnSlide(stSlice, stSlice$regionColorsVector, title='Regions', file=paste0(figures.dir, '/slice-regions.pdf'))


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
                  file=paste0(figures.dir, '/slice-log2-total-counts.pdf'))

plotScatter(cellNumbers[, 'Epithelial'], log2(totalCounts), xlab='# epithelial cells', ylab='# reads [log2]', fit='lqs',
            col=transparentRGB(stSlice$regionColorsVector, 85),
            cor.x=0.7, cor.y=0.2, file=paste0(figures.dir, '/scatter-number-epithelial-cells-vs-number-reads.pdf')) 


fit <- lm(log2(totalCounts) ~ cellNumbers[, 'Epithelial'] + cellNumbers[, 'Fibroblasts'] + cellNumbers[, 'Others'])
summary(fit)
summary(fit)$coef


##Counts adjustment
totalCountsLog2 <- log2(totalCounts)
stSlice$adjCounts <- t(apply(stSlice$rawCounts, 1, function(z) {
  lqs(z ~ totalCountsLog2)$residuals + mean(z)}))
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


##Quit
sessionInfo()
q(save='no')
