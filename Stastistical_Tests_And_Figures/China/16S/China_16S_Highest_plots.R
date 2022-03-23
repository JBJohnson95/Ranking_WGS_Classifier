rm(list=ls())
library(data.table)
library(ecodist)
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
  maxCors_1[1] <- NA
  
  #Qiime data is labelled as _2
  
  maxCors_2 <- vector(length=ncol(ordered_myT2))
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
  maxCors_2[1] <- NA
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/China/16S/")

  pdf(paste("China_Normalized_RDP_Highest_Spearman_",input[k],".pdf",sep=""))
  plot(ordered_Normalized_Mean_myT1, maxCors_1,xlab = "Log Average value of Taxa",
  ylab="Rho",main=paste("China RDP Spearman test at ",input[k],sep=""),ylim=c(-1,1),
  col= ifelse(ordered_FDR_Normalized_myT1>0.05,"black","orange"),pch=16 )
  legend("bottom",inset=c(0,.01),c("Insignificant","Significant for RDP"), pch = c(16, 16),cex=0.8,
         col=c("black","orange") )
  dev.off()
  
  pdf(paste("China_Normalized_Qiime_Highest_Spearman_",input[k],".pdf",sep=""))
  plot(ordered_Normalized_Mean_myT2, maxCors_2,xlab = "Log Average value of Taxa",
       ylab="Rho",main=paste("China Qiime Spearman test at ",input[k],sep=""),ylim=c(-1,1),
       col= ifelse(ordered_FDR_Normalized_myT2>0.05,"black","dark green"),pch=16 )
  legend("bottom",inset=c(0,.01),c("Insignificant","Significant for Qiime"),pch = c(16, 16),cex=0.8,
         col=c("black","dark green") )
  dev.off()

}

