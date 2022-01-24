rm(list=ls())
library(data.table)
library(ecodist)
library(stringr)
library(gdata)

input=c("phylum","class","order","family","genus")

setwd('/Users/James/Desktop/Dissertation/metadata/')
in_meta=read.table('Vanderbilt_metadata_all.csv',header = TRUE,sep=",")
metadata_columns=c("Sample_Names","sample_type","sex","BMI","collection_season")
colnames(in_meta)[which(names(in_meta) == "sample_id")] <- "Sample_Names"

meta=in_meta[metadata_columns]

for (k in 1:length(input)){
  setwd('/Users/James/Desktop/Dissertation/Farnaz_Kraken2_Unclassified_Reads/Tables/Vanderbilt_Forward_Only/')
  
  #Kraken2 
  myT1=read.csv(file=paste("WGS_Kraken2_Vanderbilt_Forward_2020Oct02_taxaCount_",input[k],'.tsv',sep=""),sep="\t",header=TRUE)
  colnames(myT1)[which(names(myT1) == "sample_id")] <- "Sample_Names"
  reduced=merge(meta,myT1,by="Sample_Names")
  reduced=reduced[,!grepl(pattern="Unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="noname",x=colnames(reduced))]
  
  colnames(reduced)=str_replace_all(colnames(reduced),pattern="[[.]]",replacement="_")
  
  original_input_kraken2=myT1
  rm(myT1)
  rm(transposed)
  myT1=reduced
  rm(reduced)
  
  #Metaphlan2 
  myT2=read.csv(file=paste("Metaphlan2_merged_",input[k],'.tsv',sep=""),sep="\t",header=FALSE,stringsAsFactors=FALSE)
  myT2[1,which(myT2[1,]=="ID")] <- "Sample_Names"
  myT2[1,]=gsub(pattern="_reads",replacement = "",myT2[1,])
  
  #special transpose function
  transposed=transpose(myT2)
  
  colnames(transposed)=transposed[1,]
  transposed=transposed[2:nrow(transposed),]
  
  transposed[,2:ncol(transposed)]=sapply(transposed[,2:ncol(transposed)],as.numeric)
  
  reduced=merge(meta,transposed,by="Sample_Names")
  reduced=reduced[,!grepl(pattern="Unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="unclassified",x=colnames(reduced))]
  reduced=reduced[,!grepl(pattern="noname",x=colnames(reduced))]
  
  colnames(reduced)=str_replace_all(colnames(reduced),pattern="[[.]]",replacement="_")
  
  original_input_metaphlan2=myT2
  rm(myT2)
  myT2=reduced
  rm(reduced)
  
  Mean_myT1=colMeans(myT1[,6:ncol(myT1)])
  Mean_myT2=colMeans(myT2[,6:ncol(myT2)])
  
  ordered_Mean_myT1=Mean_myT1[order(-Mean_myT1)]
  ordered_Mean_myT2=Mean_myT2[order(-Mean_myT2)]
  
  ordered_myT1=myT1[order(-Mean_myT1)+5]
  ordered_myT2=myT2[order(-Mean_myT2)+5]
  
  ordered_myT1_colnames=names(ordered_Mean_myT1)
  ordered_myT2_colnames=names(ordered_Mean_myT2)
  
  #Kraken2
  highest=rep(0,length(ordered_Mean_myT1)-1)
  for (i in 2:length(ordered_Mean_myT1)){
    flag=i-1
    for (j in 1:flag){
      temp=cor.test(x=ordered_myT1[,i],y=ordered_myT1[,j],method="spearman",alternative="two.sided")$estimate
      if (abs(temp)>highest[i-1]){
        highest[i-1]=temp
      }
    }
  }
  
  #Metaphlan2
  highest_2=rep(0,length(ordered_Mean_myT2)-1)
  for (i in 2:length(ordered_Mean_myT2)){
    flag=i-1
    for (j in 1:flag){
      temp=cor.test(x=ordered_myT2[,i],y=ordered_myT2[,j],method="spearman",alternative="two.sided")$estimate
      if (abs(temp)>highest_2[i-1]){
        highest_2[i-1]=temp
      }
    }
  }
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/Vanderbilt_Forward_Only/")    
  
  Normalized_myT1=read.table(file=paste('Kraken2_Vanderbilt_Forward_Only_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  Normalized_myT2=read.table(file=paste('Metaphlan2_Vanderbilt_Forward_Only_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  
  colnames(Normalized_myT1)=str_replace_all(colnames(Normalized_myT1),pattern="[[.]]",replacement="_")
  colnames(Normalized_myT2)=str_replace_all(colnames(Normalized_myT2),pattern="[[.]]",replacement="_")
  
  Normalized_Mean_myT1=log10(colMeans(10^Normalized_myT1[,6:ncol(Normalized_myT1)]-1)+1)
  Normalized_Mean_myT2=log10(colMeans(10^Normalized_myT2[,6:ncol(Normalized_myT2)]-1)+1)
  
  #Kruskal Wallis test for Kraken2
  pValuesSampleType_Normalized_myT1=vector()
  Name_Normalized_myT1 <- vector()
  
  for(i in 6:(ncol(Normalized_myT1))){
    SampleType=factor(Normalized_myT1$sample_type)
    kruskal <- kruskal.test(formula=Normalized_myT1[,i]~SampleType)
    pValuesSampleType_Normalized_myT1[i-5] <- kruskal$p.value
    Name_Normalized_myT1[i-5] <- names(Normalized_myT1)[i]
  }
  FDR_Corrected_pValuesSampleType_Normalized_myT1=p.adjust(pValuesSampleType_Normalized_myT1, method = "BH")
  
  #Kruskal Wallis test for Metaphlan2
  pValuesSampleType_Normalized_myT2=vector()
  Name_Normalized_myT2 <- vector()
  
  for(i in 6:(ncol(Normalized_myT2))){
    SampleType=factor(Normalized_myT2$sample_type)
    kruskal <- kruskal.test(formula=Normalized_myT2[,i]~SampleType)
    pValuesSampleType_Normalized_myT2[i-5] <- kruskal$p.value
    Name_Normalized_myT2[i-5] <- names(Normalized_myT2)[i]
  }
  FDR_Corrected_pValuesSampleType_Normalized_myT2=p.adjust(pValuesSampleType_Normalized_myT2, method = "BH")
  
  ordered_FDR_Normalized_myT1=FDR_Corrected_pValuesSampleType_Normalized_myT1[order(-Mean_myT1)]
  ordered_FDR_Normalized_myT2=FDR_Corrected_pValuesSampleType_Normalized_myT2[order(-Mean_myT2)]
  
  #Add five to ignore metadata
  ordered_Normalized_myT1=Normalized_myT1[,order(-Mean_myT1)+5]
  ordered_Normalized_myT2=Normalized_myT2[,order(-Mean_myT2)+5]
  
  Original_Order_Normalized_Mean_myT1=Normalized_Mean_myT1[order(-Mean_myT1)]
  Original_Order_Normalized_Mean_myT2=Normalized_Mean_myT2[order(-Mean_myT2)]
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/Vanderbilt/")
  
  pdf(paste("Vanderbilt_Raw_Normalized_Kraken2_Spearman_",input[k],".pdf",sep=""))
  plot(Original_Order_Normalized_Mean_myT1[2:length(Original_Order_Normalized_Mean_myT1)],highest,xlab = "Log Mean",
       ylab="Rho",main=paste("Vanderbilt Kraken2 Spearman test at ",input[k]," level",sep=""),ylim=c(-1,1),
       col= ifelse(ordered_FDR_Normalized_myT1[-c(1)]>0.05,"black","blue"),pch=16)
  legend("bottom",inset=c(0,max(highest)*.1),c("Insignificant","Significant for Kraken2"), pch = c(16, 16),cex=0.8,
         col=c("black","blue") )
  dev.off()
  
  pdf(paste("Vanderbilt_Raw_Normalized_Metaphlan2_Spearman_",input[k],".pdf",sep=""))
  plot(Original_Order_Normalized_Mean_myT2[2:length(Original_Order_Normalized_Mean_myT2)],highest_2,xlab = "Log Mean",
       ylab="Rho",main=paste("Vanderbilt Metaphlan2 Spearman test at ",input[k]," level",sep=""),ylim=c(-1,1),
       col= ifelse(ordered_FDR_Normalized_myT2[-c(1)]>0.05,"black","red"),pch=16)
  legend("bottom",inset=c(0,max(highest_2)*.1),c("Insignificant","Significant for Metaphlan2"),
         pch = c(16, 16),cex=0.8,col=c("black","red") )
  dev.off()
  
}