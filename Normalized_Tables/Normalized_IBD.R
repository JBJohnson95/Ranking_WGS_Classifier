rm(list=ls())
library(gdata)
library(data.table)
library(stringr)

input=c("phylum","class","order","family","genus")

setwd('/Users/James/Desktop/Dissertation/metadata/')
in_meta=read.table('ibd_metadata.csv',header = TRUE,sep=",")
colnames(in_meta)[which(names(in_meta) == "Sample.ID")] <- "Sample_Names"
metadata_columns=c("Sample_Names","diagnosis","IBD")
meta=in_meta[metadata_columns]
meta=meta[!is.na(meta$IBD),]

#Metaphlan2
for (j in 1:length(input)){
  setwd('/Users/James/Desktop/Dissertation/Farnaz_Kraken2_Unclassified_Reads/Tables/IBD/')
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
  
  Total=sum(rowSums(reduced[,4:ncol(reduced)]))
  Normalized=reduced
  for (i in 1:nrow(reduced)){
    Normalized[i,4:ncol(reduced)]=log10((reduced[i,4:ncol(reduced)]*1000000)/(rowSums(reduced[,4:ncol(reduced)])[i])+1)
  }
  setwd('/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/IBD/')
  
  fwrite(x=Normalized,file=paste("Metaphlan2_IBD_",input[j],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
}

#Kraken2
for (j in 1:length(input)){
  setwd('/Users/James/Desktop/Dissertation/Farnaz_Kraken2_Unclassified_Reads/Tables/IBD/')
  keep(input,meta,j,sure =TRUE)
  myT=read.csv(file=paste("kraken_merged_abundance_table_",input[j],'.txt',sep=""),sep="\t",header=FALSE,stringsAsFactors=FALSE)
  myT[1,which(myT[1,]=="ID")] <- "Sample_Names"
  myT[1,]=gsub(pattern="_1_kneaddata",replacement = "",myT[1,])

  #special transpose function
  transposed=transpose(myT)

  colnames(transposed)=transposed[1,]
  transposed=transposed[2:nrow(transposed),]

  transposed[,2:ncol(transposed)]=sapply(transposed[,2:ncol(transposed)],as.numeric)

  reduced=merge(meta,transposed,by="Sample_Names")
  reduced=reduced[,!grepl(pattern="Unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="noname",x=colnames(reduced))]

  Total=sum(rowSums(reduced[,4:ncol(reduced)]))
  Normalized=reduced
  for (i in 1:nrow(reduced)){
    Normalized[i,4:ncol(reduced)]=log10((reduced[i,4:ncol(reduced)]*1000000)/(rowSums(reduced[,4:ncol(reduced)])[i])+1)
  }
  setwd('/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/IBD/')

  fwrite(x=Normalized,file=paste("Kraken2_IBD_",input[j],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
}

