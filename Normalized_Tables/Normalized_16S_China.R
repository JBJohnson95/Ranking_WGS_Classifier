rm(list=ls())
library(ecodist)
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
  myT1=merge(meta,myT1,by="Sample_ID")
  myT1=myT1[order(myT1$True_order),]
  myT1$True_order=NULL
  
  boolean=c(TRUE,TRUE)
  boolean=append(boolean,colSums(myT1[,3:ncol(myT1)])>0)
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
  myT2=merge(meta,myT2,by="Sample_ID")
  myT2=myT2[order(myT2$True_order),]
  myT2$True_order=NULL
  
  boolean=c(TRUE,TRUE)
  boolean=append(boolean,colSums(myT2[,3:ncol(myT2)])>0)
  myT2=myT2[,boolean]
  
  Table1=myT1
  Table2=myT2
  
  Total=sum(rowSums(Table1[,3:ncol(Table1)]))
  RDP_Normalized=Table1
  for (i in 1:nrow(Table1)){
    RDP_Normalized[i,3:ncol(Table1)]=log10((Table1[i,3:ncol(Table1)]*1000000)/(rowSums(Table1[,3:ncol(Table1)])[i])+1)
  }
  
  Total=sum(rowSums(Table2[,3:ncol(Table2)]))
  Qiime_Normalized=Table2
  for (i in 1:nrow(Table2)){
    Qiime_Normalized[i,3:ncol(Table2)]=log10((Table2[i,3:ncol(Table2)]*1000000)/(rowSums(Table2[,3:ncol(Table2)])[i])+1)
  }
  
  #Normalized Tables
  
  Normalized_myT1=RDP_Normalized
  Normalized_myT2=Qiime_Normalized
  
  setwd('/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/16S_China/')
  
  fwrite(x=Normalized_myT1,file=paste("RDP_China_",input[k],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
  
  fwrite(x=Normalized_myT2,file=paste("Qiime_China_",input[k],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
  
}
