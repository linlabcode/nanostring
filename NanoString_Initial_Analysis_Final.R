#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = T)
print(args)
if (length(args) != 4) {
  stop('Need all 4 arguments (rcc path, group string, output folder, analysis name', call. = F)
} 

rcc_path = args[1]
group_string = args[2]
output_folder = args[3]
analysis_name = args[4]

#quit()

setwd('./')

library(affy)
library(graphics)
#library(qualV)

#==================================================================
#====================USER SUPPLIED PARAMETERS======================
#==================================================================

#rcc_path = '20190213_danpcp4hr12sample_RCC.txt' #command line supplied parameter

#group_string = "_01,_02,_03;_04,_05,_06;_07,_08,_09;_10,_11,_12" #to be replaced by a command line parameter where you
#give a unique substring for each sample. commas separate samples within a group, semicolon separates groups
#group_string = "_01;_02;_03;..." if you wanted it to be no groups
#need to make sure that each substring pulls out exactly one sample, if not throw error

#output folder
#output_folder = './'

#set the analysis name
#analysis_name = 'pan_cancer_4hr'

#n_ercc is always the same and they are always at bottom of table
n_ercc = 14
#==================================================================
#===============BRINGING IN DATA AND SANITY CHECKING===============
#==================================================================

#pulling in the table
rcc_table = read.delim(rcc_path, stringsAsFactors = F)

#interpreting the user supplied group string
group_list = unlist(strsplit(group_string,';'))

#defining the number of groups
n_groups = length(group_list)

#getting the sample names from the rcc table
sample_names = colnames(rcc_table)[4:ncol(rcc_table)]

#making a list of all sample names for later :)
all_samples_list = c()

#figuring out group column structure within the rcc table
col_list = list() #<- main list that will allow you to quickly find which columns correspond to which groups
#this is a list and not a matrix to allow raggedy group structure and design

#looping through each group of sample names
for(group in group_list){
  group_samples = unlist(strsplit(group,',')) # getting individual sample names
  group_cols = c() # to assign the columns for the whole group
  for(name in group_samples){ #looping through each name unique substring
    
    matching_samples = grep(name,sample_names) #getting columns of samples that match that
    if(length(matching_samples) > 1){ #making sure we don't match multiple samples
      print(samples_names[matching_samples])
      stop('YOU DONE MESSED UP A AARON')
      #quit() #properly throw a R error
    }else if(length(matching_samples)==0){ #making sure we match at least once
      stop('IS THERE A JAKWELLAN HERE')
      #quit()  
    }
    all_samples_list = c(all_samples_list,sample_names[matching_samples]) #all of your sample names
    group_cols = c(group_cols,matching_samples) #all of the columns for a group
  }
  col_list = c(col_list,list(group_cols)) #adding columns to the overall list
}

if(length(unique(all_samples_list)) != length(sample_names)){ #sanity check that # of samples is consistent
  stop('BALAKE')
  #throw an error
}else{
  print('ALL SAMPLES PROPERLY ASSIGNED')
}


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

#returns a vector of the concentrations for each ercc probe
plot_ercc <- function(expr_m,name){
  #first get the erccRows
  pos_rows = grep("POS_",rownames(expr_m))
  neg_rows = grep('NEG_',rownames(expr_m))
  ercc_rows = c(pos_rows,neg_rows)
  pos_conc = c(8,128,0.125,2,32,0.5)
  conc_vector = c(pos_conc,rep(0.01,length(neg_rows)))
  exp_order = order(conc_vector)
  
  #now let's do some cute plotting
  
  #get the expression dynamic range
  y_min = 0.1
  y_max = 1.2*max(expr_m)
  #set color palette
  palette = rainbow(ncol(expr_m),alpha=0.3)
  
  #first a blank plot
  plot(log10(conc_vector),log2(expr_m[ercc_rows,1]),
       cex=0,xlab='log10 attomoles/ul',ylab='log2 expression (counts)',
       main=paste(name,' spike-in expression',sep=""),xaxt='n',
       xlim = c(log10(0.01),log10(150)),ylim = c(log2(y_min),log2(y_max)))
  
  axis(1,log10(c(0.01,pos_conc)),labels = c('NA',as.character(pos_conc)))
  
  for(i in 1:ncol(expr_m)){
    #color = add.alpha(i,0.2)
    points(log10(conc_vector),log2(expr_m[ercc_rows,i]),pch=19,col =add.alpha(i,0.2),cex=1)
    lines(loess.smooth(log10(conc_vector[exp_order]),log2(expr_m[ercc_rows[exp_order],i])),lwd=2,col=i)	
  }				
  legend(log10(0.01),1*max(log2(expr_m[ercc_rows,])),colnames(expr_m),col=1:ncol(expr_m),lwd=2)
}


#panel function to do a scatter with a red diagonal line
panel.awesome <- function(x, y, col = par("col"), bg = NA, pch = par("pch"), 
                          cex = 1, col.smooth = "red", span = 2/3, iter = 3,n_ercc=14, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex,ylab='log2 expression (a.u.)',xlab='log2 expression (a.u.)')
  points(rev(x)[1:8], rev(y)[1:8], pch = 18, col = 'black',cex=1.5)
  points(rev(x)[9:14], rev(y)[9:14], pch = 18, col = 'blue',cex=1.5)
  
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    #lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = 'red', ...)
    abline(a=0,b=1,lwd=2,col='red')
}

#panel function to do correlation
#adapted from http://www.r-bloggers.com/five-ways-to-visualize-your-pairwise-comparisons/
panel.cor <- function(x,y,digits=2,prefix="",n_ercc=14,...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr=c(0,1,0,1))
  r <- abs(cor(x[1:(length(x)-n_ercc)],y[1:(length(y) - n_ercc)],method='spearman',use='complete'))
  txt <- round(r,4)
  txt <- paste(prefix,txt,sep="")
  cex <- 2
  test <- cor.test(x[1:(length(x)-n_ercc)],y[1:(length(y) - n_ercc)],method='spearman',use='complete')
  Signif <- symnum(test$p.value,corr=FALSE,na=FALSE,cutpoints = c(0,0.001,0.01,0.05,0.1,1),symbols = c("***","**","*",".","N.S"))
  text(0.5,0.5,txt,cex=cex*r)
  text(.8,.8,Signif,cex=cex,col=2)
  
}


#==================================================================
#====================MAKING PROCESSED TABLES=======================
#==================================================================

norm_table = loessNormalizeNano(rcc_table) #normalize based on ERCC
nanoM_raw = norm_table$raw
nanoM_norm = norm_table$norm

par(mar = c(1,1,1,1))


#clean up the raw and normalized tables
#output them into a meaningful table format
#rename columns based off of sample names

formatted_raw= matrix(nrow= nrow(nanoM_raw),ncol = length(all_samples_list))
formatted_norm= matrix(nrow= nrow(nanoM_raw),ncol = length(all_samples_list))

column_ticker = 1
column_names_list = c()
for(i in 1:n_groups){
  
  #first generate a group name
  group_name = paste('group_',as.character(i),sep='')
  #get the names of the samples in the group
  group_samples = unlist(strsplit(group_list[i],','))
  group_cols = unlist(col_list[i])
  if(length(group_samples) != length(group_cols)){
    print('GO SEE PRINCIPAL OSHAG HENNESEY')
  }
  group_size = length(group_samples)
  for(j in 1:group_size){
    formatted_sample_name = paste(group_name,'_',group_samples[j],sep='')
    column_names_list = c(column_names_list,formatted_sample_name)
    
    formatted_raw[,column_ticker] = as.numeric(nanoM_raw[,group_cols[j]])
    formatted_norm[,column_ticker] = round(as.numeric(nanoM_norm[,group_cols[j]]),2)
    column_ticker = column_ticker+1
  }
  
}

#set a floor of 0.1 for both datasets
formatted_norm[which(formatted_norm <= 0.1)] <- 0.1
formatted_raw[which(formatted_raw <= 0.1)] <- 0.1


colnames(formatted_raw) = column_names_list
rownames(formatted_raw) = rownames(nanoM_raw)
#formatted_raw = formatted_raw[order(row.names(formatted_raw)),]

colnames(formatted_norm) = column_names_list
rownames(formatted_norm) = rownames(nanoM_norm)
#formatted_norm = formatted_norm[order(row.names(formatted_norm)),]


#==================================================================
#===================WRITING PROCESSED TABLES=======================
#==================================================================


#make a sample key
sample_key_path = paste(output_folder,analysis_name,'_sample_key.txt',sep='')
sample_key_table = cbind(all_samples_list,column_names_list)
write.table(sample_key_table,sample_key_path,quote=FALSE,sep='\t',row.names=FALSE)

#writing out the raw and normalized expression data
raw_path = paste(output_folder,analysis_name,'_raw.txt',sep='')
norm_path = paste(output_folder,analysis_name,'_norm.txt',sep='')

output_raw_table = cbind(rownames(formatted_raw),formatted_raw)
colnames(output_raw_table) = c('GENE_NAME',colnames(formatted_raw))

output_norm_table = cbind(rownames(formatted_norm),formatted_norm)
colnames(output_norm_table) = c('GENE_NAME',colnames(formatted_norm))

write.table(output_raw_table,raw_path,quote=FALSE,sep='\t',row.names=FALSE)
write.table(output_norm_table,norm_path,quote=FALSE,sep='\t',row.names=FALSE)

#hard code which rows are endogenous and which rows are ercc



#==================================================================
#=======================QUALITY CONTROL============================
#==================================================================

#1. do replicate correlate?

#set the axis limits
axis_min = 0.1
axis_max = 1.2*max(formatted_norm)
replicate_scatter_path = paste(output_folder,analysis_name,'_replicate_scatter_norm.pdf',sep='')

pdf(file=replicate_scatter_path,width =9,height=10)

start_col = 1
for(i in 1:n_groups){
  group_samples = unlist(strsplit(group_list[i],','))
  group_size = length(group_samples)
  group_columns = start_col: (start_col + group_size - 1)
  start_col = start_col + group_size
  group_name = paste('group_',as.character(i),sep='')
  figure_title = paste(group_name,' pairwise scatter',sep='')
  pairs(log2(formatted_norm[,group_columns]),lower.panel=panel.awesome,upper.panel=panel.cor,
        pch=19,col=rgb(0.5,0.5,0.5,0.4),cex.labels =1.5,main=figure_title,
        xlim =c(log2(axis_min),log2(axis_max)),ylim = c(log2(axis_min),log2(axis_max)))
  
}

dev.off()


axis_min = 0.1
axis_max = 1.2*max(formatted_raw)
replicate_scatter_path = paste(output_folder,analysis_name,'_replicate_scatter_raw.pdf',sep='')

pdf(file=replicate_scatter_path,width =9,height=10)

start_col = 1
for(i in 1:n_groups){
  group_samples = unlist(strsplit(group_list[i],','))
  group_size = length(group_samples)
  group_columns = start_col: (start_col + group_size - 1)
  start_col = start_col + group_size
  group_name = paste('group_',as.character(i),sep='')
  figure_title = paste(group_name,' pairwise scatter',sep='')
  pairs(log2(formatted_raw[,group_columns]),lower.panel=panel.awesome,upper.panel=panel.cor,
        pch=19,col=rgb(0.5,0.5,0.5,0.4),cex.labels =1.5,main=figure_title,
        xlim =c(log2(axis_min),log2(axis_max)),ylim = c(log2(axis_min),log2(axis_max)))
  
}

dev.off()




#2. did the ercc normalization work?

ercc_path = paste(output_folder,analysis_name,'_ercc_line.pdf',sep='')
pdf(file=ercc_path,width = 10,height= 8)
plot_ercc(formatted_raw,'raw')
plot_ercc(formatted_norm,'norm')
dev.off()



#3. global differences in expression between samples?

#we want to do this w/ the norm data

#set expression cutoff using the pos control
pos_rows = grep('POS_',rownames(formatted_norm))
expr_cut = min(apply(formatted_norm[pos_rows,],1,mean))

#expressed genes are expressed in at least 1 sample above cutoff
endogenous_rows = 1:(nrow(formatted_norm) - n_ercc)

expr_rows = which(apply(formatted_norm[endogenous_rows,],1,max) > expr_cut)

formatted_norm = formatted_norm[expr_rows,]
formatted_raw = formatted_raw[expr_rows,]

formatted_norm = as.data.frame(formatted_norm)
formatted_raw = as.data.frame(formatted_raw)

for(i in 1:n_groups){
  group_name = paste('group_', as.character(i), sep = '')
  temp_df = formatted_norm[, grep(group_name, names(formatted_norm), value = T)]
  formatted_norm$avg = rowMeans(temp_df)
  colnames(formatted_norm)[length(formatted_norm)] = paste(group_name, '_avg', sep = '')
}

for(i in 1:n_groups){
  group_name = paste('group_', as.character(i), sep = '')
  temp_df = formatted_raw[, grep(group_name, names(formatted_raw), value = T)]
  formatted_raw$avg = rowMeans(temp_df)
  colnames(formatted_raw)[length(formatted_raw)] = paste(group_name, '_avg', sep = '')
}

formatted_raw = formatted_raw[order(row.names(formatted_raw)),]
formatted_norm = formatted_norm[order(row.names(formatted_norm)),]

formatted_raw_table = cbind(rownames(formatted_raw), formatted_raw)
colnames(formatted_raw_table) = c('Gene_ID', colnames(formatted_raw))

formatted_norm_table = cbind(rownames(formatted_norm), formatted_norm)
colnames(formatted_norm_table) = c('Gene_ID', colnames(formatted_norm))

write.table(formatted_raw_table, './Final_Formatted_Raw.csv', row.names = F, col.names = T, quote = F, sep = ',')
write.table(formatted_norm_table, './Final_Formatted_Norm.csv', row.names = F, col.names = T, quote = F, sep = ',')

#create a color vector based on groups
color_vector = c()

palette = rainbow(n_groups,alpha=0.3)

for(i in 1:n_groups){
  group_samples = unlist(strsplit(group_list[i],','))
  group_size = length(group_samples)
  color_vector = c(color_vector,rep(palette[i],group_size))
}


box_path = paste(output_folder,analysis_name,'_norm_boxplots.pdf',sep='')
pdf(file=box_path,width = 8,height =6)
par(mar=c(10,2,2,2))
boxplot(log2(formatted_norm[expr_rows,]),cex=0,col=color_vector,las=2,ylab='log2 counts')
abline(h=log2(expr_cut))
dev.off()

#4. did erccs capture expression dynamic range? (do we even detect ERCCs)




