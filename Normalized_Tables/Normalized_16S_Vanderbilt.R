rm(list=ls())
library(data.table)
library(stringr)

input=c("phylum","class","order","family","genus")

for (k in 1:length(input)){
  setwd(dir = "/Users/James/Desktop/Dissertation/Vanderbilt_16S data/")
  myT1=read.table(file=paste('rdpClassifications_',input[k],'Level.txt',sep=""),header=FALSE,sep = "\t")
  myT1=transpose(myT1)
  colnames(myT1)=myT1[1,]
  myT1=myT1[-1,]
  myT1=myT1[order(myT1[,1]),]
  for(i in 2:ncol(myT1)){
    myT1[,i]=as.numeric(unlist(myT1[,i]))
  }
  
  myT1=myT1[,!grepl(pattern="Unclassified",x=colnames(myT1))]
  myT1=myT1[,!grepl(pattern="unclassified",x=colnames(myT1))]
  myT1=myT1[,!grepl(pattern="noname",x=colnames(myT1))]
  colnames(myT1)=gsub(pattern="[[]",x=colnames(myT1),replacement="")
  colnames(myT1)=gsub(pattern="[]]",x=colnames(myT1),replacement="")
  colnames(myT1)=gsub(pattern="Sample ID",replacement = "Sample_ID",x=colnames(myT1))
  
  
  myT2=read.csv(file=paste('qiime_',input[k],'2Level.txt',sep=""),header=FALSE,skip=1,sep = "\t",stringsAsFactors = FALSE)
  header = read.csv(file=paste('qiime_',input[k],'2Level.txt',sep=""), header = F, nrows = 1,sep="\t",stringsAsFactors = FALSE)
  headers=sapply(header,as.character)
  new_header=append("Sample ID",headers)
  myT2=rbind(new_header,myT2)
  myT2=transpose(myT2)
  colnames(myT2)=myT2[1,]
  myT2=myT2[-1,]
  myT2=myT2[order(myT2[,1]),]
  for(i in 2:ncol(myT2)){
    myT2[,i]=as.numeric(unlist(myT2[,i]))
  }
  
  myT2=myT2[,!grepl(pattern="Unclassified",x=colnames(myT2))]
  myT2=myT2[,!grepl(pattern="unclassified",x=colnames(myT2))]
  myT2=myT2[,!grepl(pattern="noname",x=colnames(myT2))]
  colnames(myT2)=gsub(pattern="[[]",x=colnames(myT2),replacement="")
  colnames(myT2)=gsub(pattern="[]]",x=colnames(myT2),replacement="")
  colnames(myT2)=gsub(pattern="Sample ID",replacement = "Sample_ID",x=colnames(myT2))
  
  myT2=myT2[,grep(colnames(myT2),pattern="\\.1",invert=TRUE)]
  
  setwd(dir="/Users/James/Desktop/Dissertation/metadata/")
  meta=read.table(file="Vanderbilt_metadata_all.csv",header=TRUE,sep = ",")
  
  Correct_meta=meta[(meta$sample_id %in% myT2[,1]),]
  Sample_Type=Correct_meta$sample_type
  
  Sub_meta=Correct_meta[,c("sample_id","sample_type")]
  colnames(Sub_meta)=c("Sample_ID","sample_type")
  
  #Has correct sample rows
  
  Table1=myT1[ myT1[,1] %in% Correct_meta$sample_id,]
  Table2=myT2[ myT2[,1] %in% Correct_meta$sample_id,]
  
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
  
  myT1_Normalized=merge(Sub_meta,RDP_Normalized,by="Sample_ID")
  myT2_Normalized=merge(Sub_meta,Qiime_Normalized,by="Sample_ID")
  
  setwd('/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/16S_Vanderbilt/')
  
  fwrite(x=myT1_Normalized,file=paste("RDP_Vanderbilt_Forward_Only_",input[k],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
  fwrite(x=myT2_Normalized,file=paste("Qiime_Vanderbilt_Forward_Only_",input[k],"_Normalized.csv",sep=""),row.names = FALSE,sep=",")
  
}
  
