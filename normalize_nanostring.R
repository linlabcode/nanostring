library(affy)
library(graphics)
library(qualV)

setwd('~/Dropbox/MYCN_amp/Nanostring/')



#==================================================================
#==========================PARAMETERS==============================
#==================================================================

inputNanostringTable = read.delim('Runs/20150716_20150716-mycwt_RCC/20150717_SHEP_wildtype_myc_induction.txt')
name = '150717_shep_mycwt'



#inputNanostringTable = read.delim('Runs/20150716_20150716-mycwt_RCC/20150717_SHEP_wildtype_mycn_induction.txt')
#name = '150717_shep_mycnwt'

timePoints = c(0,4,8,24)

#==================================================================
#===========================FUNCTIONS==============================
#==================================================================




## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

loessNormalizeNano <- function(nanostringTable){
	
	erccRows = grep('ERCC',nanostringTable[,3] )

	#get a simple matrix
	nanoM = as.matrix(nanostringTable[,4:ncol(nanostringTable)])
	rownames(nanoM) = nanostringTable[,2]
	nanoM_norm = loess.normalize(nanoM,subset=erccRows,log.it=FALSE,family.loess='gaussian')
	return(list('norm'=nanoM_norm,'raw'=nanoM))
}

plotERCC <- function(nanoM_norm,name){
	
	
	negRows = grep('NEG',rownames(nanoM_norm))
	
	posRows = grep('POS',rownames(nanoM_norm))
	
	posConc = c(32,8,0.5,2,0.125,128)
	
	posOrder= order(posConc)
	
	posRows = posRows[posOrder]
	posConc = posConc[posOrder]
	
	erccRows = c(negRows,posRows)
	concVector = c(rep(0,length(negRows)),posConc)
	
	par(mfrow=c(4,3))
	for(i in erccRows){
		barplot(nanoM_norm[i,],main=paste(rownames(nanoM_norm)[i],name))
		}		
	
}

#==================================================================
#======================MAKING OUTPUTS==============================
#==================================================================

foo = loessNormalizeNano(inputNanostringTable)
nanoM_raw = foo$raw
nanoM_norm = foo$norm


pdf
plotERCC(nanoM_raw,'raw')
plotERCC(nanoM_norm,'norm')

className = 'rRNA'

plotByClass <- function(nanoM_norm,className){
	rRNArows = grep(className,rownames(nanoM_norm))
	nano_rRNA = nanoM_norm[rRNArows,]
	replicates = c('rep1','rep2','rep3')
	rRNA_norm_matrix = matrix(ncol = ncol(nanoM_norm),nrow=length(rRNArows))
	
	plot(timePoints,timePoints,cex=0,ylim=c(0,3),main=className)
	for(rep in replicates){
		
		repCols = grep(rep,colnames(nanoM_norm))
		nano_rRNA[,repCols] = nano_rRNA[,repCols]/nano_rRNA[,repCols[1]]

		
	}
	ticker = 1
	
	for(row in 1:nrow(nano_rRNA)){
		
		for(rep in replicates){
			repCols = grep(rep,colnames(nanoM_norm))
			lines(timePoints,nano_rRNA[row,repCols],type='b',col=ticker)
		}
		ticker = ticker+1	
			
	}
}



linePlotFileName = paste(name,'_linePlots.pdf',sep='')
pdf(file=linePlotFileName,width = 5,height = 5)
plotByClass(nanoM_norm,'rRNA')
plotByClass(nanoM_norm,'tRNA')
plotByClass(nanoM_norm,'mRNA')
dev.off()


