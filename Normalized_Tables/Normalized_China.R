rm(list=ls())
library(gdata)
library(data.table)
library(stringr)

input=c("phylum","class","order","family","genus")

setwd('/Users/James/Desktop/Dissertation/metadata/')
meta_initial=read.table('china_rural_urban_metadata.csv',header = TRUE,sep=",")
meta=data.frame(meta_initial$sample.ID,meta_initial$ruralurban)
colnames(meta)=c("Sample","ruralurban")

#Metaphlan2

for (j in 1:length(input)){
  setwd('/Users/James/Desktop/Dissertation/Farnaz_Kraken2_Unclassified_Reads/Tables/China/')
  keep(input,meta,j,sure =TRUE)
  
  myT=read.csv(file=paste("Metaphlan2_merged_",input[j],'.tsv',sep=""),sep="\t",header=FALSE,stringsAsFactors=FALSE)
  myT[1,which(myT[1,]=="ID")] <- "Sample"
  myT[1,]=gsub(pattern="_reads",replacement = "",myT[1,])
  
  #special transpose function
  transposed=transpose(myT)
  
  colnames(transposed)=transposed[1,]
  transposed=transposed[2:nrow(transposed),]
  
  transposed[,2:ncol(transposed)]=sapply(transposed[,2:ncol(transposed)],as.numeric)
  
  reduced=merge(meta,transposed,by="Sample")
  reduced=reduced[,!grepl(pattern="Unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="noname",x=colnames(reduced))]
  
  Total=sum(rowSums(reduced[,3:ncol(reduced)]))
  Normalized=reduced
  for (i in 1:nrow(reduced)){
    Normalized[i,3:ncol(reduced)]=log10((reduced[i,3:ncol(reduced)]*1000000)/(rowSums(reduced[,3:ncol(reduced)])[i])+1)
  }
  
  setwd('/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/China/')

  fwrite(x=Normalized,file=paste("Metaphlan2_China_",input[j],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
}

#Kraken2

for (j in 1:length(input)){
  setwd('/Users/James/Desktop/Dissertation/Farnaz_Kraken2_Unclassified_Reads/Tables/China/')
  keep(input,meta,j,sure =TRUE)
  myT=read.csv(file=paste("WGS_Forward_Kraken2_China_2021Mar14_taxaCount_",input[j],'.tsv',sep=""),sep="\t",header=TRUE)
  colnames(myT)[which(names(myT) == "sample.ID")] <- "Sample"
  reduced=merge(meta,myT,by="Sample")
  reduced=reduced[,!grepl(pattern="Unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="noname",x=colnames(reduced))]
  
  Total=sum(rowSums(reduced[,3:ncol(reduced)]))
  Normalized=reduced
  for (i in 1:nrow(reduced)){
    Normalized[i,3:ncol(reduced)]=log10((reduced[i,3:ncol(reduced)]*1000000)/(rowSums(reduced[,3:ncol(reduced)])[i])+1)
  }
  setwd('/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/China/')

  fwrite(x=Normalized,file=paste("Kraken2_China_",input[j],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
}