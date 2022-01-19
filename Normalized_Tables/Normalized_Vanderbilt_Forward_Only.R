rm(list=ls())
library(gdata)
library(data.table)
library(stringr)

input=c("phylum","class","order","family","genus")

setwd('/Users/James/Desktop/Dissertation/metadata/')
in_meta=read.table('Vanderbilt_metadata_all.csv',header = TRUE,sep=",")
metadata_columns=c("Sample_Names","sample_type","sex","BMI","collection_season")
colnames(in_meta)[which(names(in_meta) == "sample_id")] <- "Sample_Names"

meta=in_meta[metadata_columns]

#Metaphlan2
for (j in 1:length(input)){
  setwd('/Users/James/Desktop/Dissertation/Farnaz_Kraken2_Unclassified_Reads/Tables/Vanderbilt_Forward_Only/')
  keep(input,meta,j,sure =TRUE)
  myT=read.csv(file=paste("Metaphlan2_merged_",input[j],'.tsv',sep=""),sep="\t",header=FALSE,stringsAsFactors=FALSE)
  myT[1,which(myT[1,]=="ID")] <- "Sample_Names"
  myT[1,]=gsub(pattern="_reads",replacement = "",myT[1,])
  
  #special transpose function
  transposed=transpose(myT)

  colnames(transposed)=transposed[1,]
  transposed=transposed[2:nrow(transposed),]
  
  transposed[,2:ncol(transposed)]=sapply(transposed[,2:ncol(transposed)],as.numeric)
  
  reduced=merge(meta,transposed,by="Sample_Names")
  reduced=reduced[,!grepl(pattern="Unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="noname",x=colnames(reduced))]
  
  Total=sum(rowSums(reduced[,6:ncol(reduced)]))
  Normalized=reduced
  for (i in 1:nrow(reduced)){
    Normalized[i,6:ncol(reduced)]=log10((reduced[i,6:ncol(reduced)]*Total)/(rowSums(reduced[,6:ncol(reduced)])[i]*nrow(reduced))+1)
  }
  setwd('/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/Vanderbilt_Forward_Only/')
  
  fwrite(x=Normalized,file=paste("Metaphlan2_Vanderbilt_Forward_Only_",input[j],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
}

#Kraken2
for (j in 1:length(input)){
  setwd('/Users/James/Desktop/Dissertation/Farnaz_Kraken2_Unclassified_Reads/Tables/Vanderbilt_Forward_Only/')
  keep(input,meta,j,sure =TRUE)
  myT=read.csv(file=paste("WGS_Kraken2_Vanderbilt_Forward_2020Oct02_taxaCount_",input[j],'.tsv',sep=""),sep="\t",header=TRUE)
  colnames(myT)[which(names(myT) == "sample_id")] <- "Sample_Names"
  reduced=merge(meta,myT,by="Sample_Names")
  reduced=reduced[,!grepl(pattern="Unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="unclassified",x=colnames(reduced))]
  
  Total=sum(rowSums(reduced[,6:ncol(reduced)]))
  Normalized=reduced
  for (i in 1:nrow(reduced)){
    Normalized[i,6:ncol(reduced)]=log10((reduced[i,6:ncol(reduced)]*Total)/(rowSums(reduced[,6:ncol(reduced)])[i]*nrow(reduced))+1)
  }
  setwd('/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/Vanderbilt_Forward_Only/')
  
  fwrite(x=Normalized,file=paste("Kraken2_Vanderbilt_Forward_Only_",input[j],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
}
