rm(list=ls())
library(data.table)
library(stringr)

input=c("phylum","class","order","family","genus")

setwd(dir="/Users/James/Desktop/Dissertation/metadata/")
meta=read.table(file="china_rural_urban_metadata.csv",header=TRUE,sep = ",")
meta=meta[-c(81,82),]
meta=data.frame(meta$sample.ID,meta$ruralurban)
colnames(meta)=c("Sample_ID","ruralurban")
meta=meta[grepl(meta$Sample_ID,pattern = "A"),]

for (k in 1:length(input)){
  setwd(dir = "/Users/James/Desktop/Dissertation/16S_Vs_WGS/UrbanRuralChina-master/16SrRNA/inputData/RDP/")
  myT1=read.table(file=paste(input[k],'_taxaAsColumns.txt',sep=""),header=TRUE,sep = "\t")
  previous_myT1=myT1
  
  colnames(myT1)=gsub(pattern="\\.",replacement="_",colnames(myT1))
  
  myT1=myT1[,!grepl(pattern="Unclassified",x=colnames(myT1))]
  myT1=myT1[,!grepl(pattern="unclassified",x=colnames(myT1))]
  myT1=myT1[,!grepl(pattern="noname",x=colnames(myT1))]
  colnames(myT1)=gsub(pattern="[[]",x=colnames(myT1),replacement="")
  colnames(myT1)=gsub(pattern="[]]",x=colnames(myT1),replacement="")
  colnames(myT1)=gsub(pattern="sample",replacement = "Sample_ID",x=colnames(myT1))
  
  #Get rid of data that do not match WGS and reverse reads
  myT1=myT1[grep(pattern="B",myT1[,1],invert=TRUE),]
  myT1=myT1[grep(pattern="A_1",myT1[,1]),]
  
  myT1$Sample_ID=gsub(pattern="_1",replacement="",x=myT1$Sample_ID)
  myT1$True_order=1:nrow(myT1)
  temporary=merge(myT1,meta,by="Sample_ID")
  temporary=temporary[order(temporary$True_order),]
  RDP_ruralurban=temporary$ruralurban
  myT1$True_order=NULL
  
  boolean=TRUE
  boolean=append(boolean,colSums(myT1[,2:ncol(myT1)])>0)
  myT1=myT1[,boolean]
  
  setwd(dir="/Users/James/Desktop/Dissertation/16S_Vs_WGS/UrbanRuralChina-master/16SrRNA/inputData/qiime/")
  
  myT2=read.csv(file=paste('qiime_china_ruralurban_',input[k],'.tsv',sep=""),header=TRUE,sep = "\t",stringsAsFactors = FALSE)
  colnames(myT2)=gsub(pattern="\\.",replacement="_",colnames(myT2))  
  
  myT2$Sample_ID=gsub(pattern="_1",replacement="",x=myT2$Sample_ID)
  
  myT2=myT2[,!grepl(pattern="Unclassified",x=colnames(myT2))]
  myT2=myT2[,!grepl(pattern="unclassified",x=colnames(myT2))]
  myT2=myT2[,!grepl(pattern="noname",x=colnames(myT2))]
  colnames(myT2)=gsub(pattern="[[]",x=colnames(myT2),replacement="")
  colnames(myT2)=gsub(pattern="[]]",x=colnames(myT2),replacement="")
  colnames(myT2)=gsub(pattern="Sample ID",replacement = "Sample_ID",x=colnames(myT2))
  
  myT2$Sample_ID=gsub(pattern="_1",replacement="",x=myT2$Sample_ID)
  myT2$True_order=1:nrow(myT2)
  temporary=merge(myT2,meta,by="Sample_ID")
  temporary=temporary[order(temporary$True_order),]
  Qiime_ruralurban=temporary$ruralurban
  myT2$True_order=NULL
  
  boolean=TRUE
  boolean=append(boolean,colSums(myT2[,2:ncol(myT2)])>0)
  myT2=myT2[,boolean]
  
  Table1=myT1
  Table2=myT2
  
  Total=sum(rowSums(Table1[,2:ncol(Table1)]))
  RDP_Normalized=Table1
  for (i in 1:nrow(Table1)){
    RDP_Normalized[i,2:ncol(Table1)]=log10((Table1[i,2:ncol(Table1)]*1000000)/(rowSums(Table1[,2:ncol(Table1)])[i])+1)
  }
  
  Total=sum(rowSums(Table2[,2:ncol(Table2)]))
  Qiime_Normalized=Table2
  for (i in 1:nrow(Table2)){
    Qiime_Normalized[i,2:ncol(Table2)]=log10((Table2[i,2:ncol(Table2)]*1000000)/(rowSums(Table2[,2:ncol(Table2)])[i])+1)
  }
  
  #Normalized Tables
  
  Normalized_myT1=RDP_Normalized
  Normalized_myT2=Qiime_Normalized
  
  Normalized_Mean_myT1=log10(colMeans(10^Normalized_myT1[,2:ncol(Normalized_myT1)]-1)+1)
  Normalized_Mean_myT2=log10(colMeans(10^Normalized_myT2[,2:ncol(Normalized_myT2)]-1)+1)
  
  #Kruskal Wallis test for RDP
  pValuesruralurban_Normalized_myT1=vector()
  Name_Normalized_myT1 <- vector()
  
  for(i in 2:(ncol(Normalized_myT1))){
    ruralurban=factor(RDP_ruralurban)
    kruskal <- kruskal.test(formula=Normalized_myT1[,i]~ruralurban)
    pValuesruralurban_Normalized_myT1[i-1] <- kruskal$p.value
    Name_Normalized_myT1[i-1] <- names(Normalized_myT1)[i]
  }
  FDR_Corrected_pValuesruralurban_Normalized_myT1=p.adjust(pValuesruralurban_Normalized_myT1, method = "BH")
  
  #Kruskal Wallis test for Qiime
  pValuesruralurban_Normalized_myT2=vector()
  Name_Normalized_myT2 <- vector()
  
  for(i in 2:(ncol(Normalized_myT2))){
    ruralurban=factor(Qiime_ruralurban)
    kruskal <- kruskal.test(formula=Normalized_myT2[,i]~ruralurban)
    pValuesruralurban_Normalized_myT2[i-1] <- kruskal$p.value
    Name_Normalized_myT2[i-1] <- names(Normalized_myT2)[i]
  }
  FDR_Corrected_pValuesruralurban_Normalized_myT2=p.adjust(pValuesruralurban_Normalized_myT2, method = "BH")
  
  ordered_FDR_Normalized_myT1=FDR_Corrected_pValuesruralurban_Normalized_myT1[order(-Normalized_Mean_myT1)]
  ordered_FDR_Normalized_myT2=FDR_Corrected_pValuesruralurban_Normalized_myT2[order(-Normalized_Mean_myT2)]
  
  #Add two to ignore metadata
  ordered_Normalized_myT1=Normalized_myT1[,order(-Normalized_Mean_myT1)+1]
  ordered_Normalized_myT2=Normalized_myT2[,order(-Normalized_Mean_myT2)+1]
  
  ordered_Normalized_Mean_myT1=Normalized_Mean_myT1[order(-Normalized_Mean_myT1)]
  ordered_Normalized_Mean_myT2=Normalized_Mean_myT2[order(-Normalized_Mean_myT2)]
  
  ordered_myT1=myT1[,order(-Normalized_Mean_myT1)+1]
  ordered_myT2=myT2[,order(-Normalized_Mean_myT2)+1]

  ### RDP data is labelled as _1
  
  maxCors_1 <- vector(length=ncol(ordered_myT1))
  sums_1 <- apply(ordered_myT1[,1:ncol(ordered_myT1)], 2, sum)
  numCols_1 <- ncol(ordered_myT1)
  
  #Twice perfect negative correlation should be statisically unlikely or impossible thus the ideal default
  for( i in (1):numCols_1 )
  {
    maxCors_1[i] <- -2
  }
  
  #Ever decreasing lower average mean Taxa are compared against higher mean taxa starting at second most abundant
  for( i in (2):numCols_1)
  {
    stopVal <- min(i-1)
    
    for( j in (3):stopVal ){
      corVal <- cor(ordered_myT1[,i],ordered_myT1[,j],method="spearman" )
      
      if( ! is.na(corVal) ){	
        maxCors_1[i] <- max(  maxCors_1[i], corVal  )
      }
    }
  }
  maxCors_1[1] <- 1
  
  a1 <- 10
  # take the a number of top taxa and get the average depth after being ordered from largest to smallest
  a1_averageDepth <- ordered_Normalized_Mean_myT1[1:a1]
  aSum_1 <- 10^(a1_averageDepth)-1
  sumTop_1 <- sum(aSum_1-1)
  #Distribution of a number of top taxa asigned as probabiltiies 
  probs_1 = aSum_1/sumTop_1
  
  simMeans_1 <- vector()
  simCors_1 <- vector()
  
  for(i in 1:(numCols_1-a1))
  {
    backgroundErrorRate = runif(1)/1000
    
    if( backgroundErrorRate < 0) {
      backgroundErrorRate  = 0
    }
    
    # choose by relative abundance
    colIndex <- sample(1:a1,1,prob=probs_1)+1
    dataSimCol <- vector(length=nrow(ordered_Normalized_myT1))
    
    for( j in 1:nrow(ordered_Normalized_myT1)) {
      dataSimCol[j] <- rbinom(1,ceiling(10^ordered_Normalized_myT1[j,colIndex])-1,backgroundErrorRate)
    }
    
    simMeans_1[i] <- log10(mean(dataSimCol)+1) 
    
    if( is.nan(simMeans_1[i]) | simMeans_1[i] < 0 ) {
      simMeans_1[i] = 0
    }
    
    simCors_1[i] <- cor(dataSimCol,10^ordered_Normalized_myT1[,colIndex],method="spearman" )
  }
  
  #QIIME data is labelled as _2
  
  maxCors_2 <- vector(length=ncol(ordered_myT2))
  sums_2 <- apply(ordered_myT2[,1:ncol(ordered_myT2)], 2, sum)
  numCols_2 <- ncol(ordered_myT2)
  
  #Twice perfect negative correlation should be statisically unlikely or impossible thus the ideal default
  for( i in 1:numCols_2 )
  {
    maxCors_2[i] <- -2
  }
  
  #Ever decreasing lower average mean Taxa are compared against higher mean taxa starting at second most abundant
  for( i in (numCols_2:2))
  {
    stopVal <- min(i-1)
    
    for( j in stopVal:1){
      corVal <- cor(ordered_myT2[,i],ordered_myT2[,j],method="spearman" )
      
      if( ! is.na(corVal) ){	
        maxCors_2[i] <- max(  maxCors_2[i], corVal  )
      }
    }
  }
  maxCors_2[1] <- 1
  
  a2 <- 10
  # take the a number of top taxa and get the average depth after being ordered from largest to smallest
  a2_averageDepth <- ordered_Normalized_Mean_myT2[1:a2]
  aSum_2 <- 10^(a2_averageDepth)-1
  sumTop_2 <- sum(aSum_2-1)
  #Distribution of a number of top taxa asigned as probabiltiies 
  probs_2 = aSum_2/sumTop_2
  
  simMeans_2 <- vector()
  simCors_2 <- vector()
  
  for(i in (a2+1):(numCols_2))
  {
    backgroundErrorRate = runif(1)/1000
    
    if( backgroundErrorRate < 0) {
      backgroundErrorRate  = 0
    }
    
    # choose by relative abundance
    colIndex <- sample(1:a2,1,prob=probs_2)
    dataSimCol <- vector(length=nrow(ordered_Normalized_myT2))
    
    for( j in 1:nrow(ordered_Normalized_myT2)) {
      dataSimCol[j] <- rbinom(1,ceiling(10^ordered_Normalized_myT2[j,colIndex])-1,backgroundErrorRate)
    }
    
    simMeans_2[i] <- log10(mean(dataSimCol)+1) 
    
    if( is.nan(simMeans_2[i]) | simMeans_2[i] < 0 ) {
      simMeans_2[i] = 0
    }
    
    simCors_2[i] <- cor(dataSimCol,10^ordered_Normalized_myT2[,colIndex],method="spearman" )
  }
  
  maxCors_1[1] <- NA
  maxCors_2[1] <- NA
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/China/16S/")

  pdf(paste("China_Highest_Poisson_Simulated_Spearman_for_RDP_at_",input[k],".pdf",sep=""))
  plot(ordered_Normalized_Mean_myT1, maxCors_1,xlab = "Log Mean",
       ylab="Rho",main=paste("China RDP Spearman vs Simulated at ",input[k],sep=""),ylim=c(-1,1),xlim=c(0,6.5),
       col= ifelse(ordered_FDR_Normalized_myT1[-c(1)]>0.05,"black","orange"),pch=16)
  legend("bottom",c("Insignificant","Significant for RDP","Simulated Data"),
         pch = c(16, 16,16),cex=0.8,col=c("black","orange","green") )
  points(simMeans_1, simCors_1,col="green")
  dev.off()
  
  pdf(paste("China_Highest_Poisson_Simulated_Spearman_for_QIIME_at_",input[k],".pdf",sep=""))
  plot(ordered_Normalized_Mean_myT2, maxCors_2,xlab = "Log Mean",
       ylab="Rho",main=paste("China QIIME Spearman vs Simulated at ",input[k],sep=""),ylim=c(-1,1),xlim=c(0,6.5),
       col= ifelse(ordered_FDR_Normalized_myT2[-c(1)]>0.05,"black","pink"),pch=16)
  legend("bottom",c("Insignificant","Significant for QIIME","Simulated Data"),
         pch = c(16, 16,16),cex=0.8,col=c("black","pink","green") )
  points(simMeans_2, simCors_2,col="green")
  dev.off()
}
