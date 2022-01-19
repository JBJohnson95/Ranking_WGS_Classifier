rm(list=ls())
library(data.table)
library(ecodist)
library(stringr)
library(gdata)

input=c("phylum","class","order","family","genus")

setwd('/Users/James/Desktop/Dissertation/metadata/')
in_meta=read.table('ibd_metadata.csv',header = TRUE,sep=",")
colnames(in_meta)[which(names(in_meta) == "Sample.ID")] <- "Sample_Names"
metadata_columns=c("Sample_Names","diagnosis","IBD")
meta=in_meta[metadata_columns]
meta=meta[!is.na(meta$IBD),]

for (k in 1:length(input)){
  setwd('/Users/James/Desktop/Dissertation/Farnaz_Kraken2_Unclassified_Reads/Tables/IBD/')
  
  #Kraken2 
  myT1=read.csv(file=paste("kraken_merged_abundance_table_",input[k],'.txt',sep=""),sep="\t",header=FALSE,stringsAsFactors=FALSE)
  myT1[1,which(myT1[1,]=="ID")] <- "Sample_Names"
  myT1[1,]=gsub(pattern="_1_kneaddata",replacement = "",myT1[1,])
  
  #special transpose function
  transposed=transpose(myT1)
  
  colnames(transposed)=transposed[1,]
  transposed=transposed[2:nrow(transposed),]
  
  transposed[,2:ncol(transposed)]=sapply(transposed[,2:ncol(transposed)],as.numeric)
  
  reduced=merge(meta,transposed,by="Sample_Names")
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
  
  overlap_names=intersect(colnames(myT1[,4:ncol(myT1)]),colnames(myT2[,4:ncol(myT2)]))
  
  myT1_colnames=colnames(myT1[,4:ncol(myT1)])
  myT2_colnames=colnames(myT2[,4:ncol(myT2)])
  
  Mean_myT1=colMeans(myT1[,4:ncol(myT1)])
  Mean_myT2=colMeans(myT2[,4:ncol(myT2)])
  
  ordered_Mean_myT1=Mean_myT1[order(-Mean_myT1)]
  ordered_Mean_myT2=Mean_myT2[order(-Mean_myT2)]
  
  ordered_myT1=myT1[order(-Mean_myT1)+3]
  ordered_myT2=myT2[order(-Mean_myT2)+3]
  
  ordered_myT1_colnames=names(ordered_Mean_myT1)
  ordered_myT2_colnames=names(ordered_Mean_myT2)
  
  #Kraken2
  highest=rep(0,length(ordered_Mean_myT1)-1)
  for (i in 2:length(ordered_Mean_myT1)){
    if (ordered_Mean_myT1[i]>0.01){
      flag=i-1
      for (j in 1:flag){
        temp=cor.test(x=ordered_myT1[,i],y=ordered_myT1[,j],method="spearman",alternative="two.sided")$estimate
        if (abs(temp)>highest[i-1]){
          highest[i-1]=temp
        }
      }
    }
  }
  
  #Metaphlan2
  highest_2=rep(0,length(ordered_Mean_myT2)-1)
  for (i in 2:length(ordered_Mean_myT2)){
    if (ordered_Mean_myT2[i]>0.01){
      flag=i-1
      for (j in 1:flag){
        temp=cor.test(x=ordered_myT2[,i],y=ordered_myT2[,j],method="spearman",alternative="two.sided")$estimate
        if (abs(temp)>highest_2[i-1]){
          highest_2[i-1]=temp
        }
      }
    }
  }
  
  setwd('/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/IBD/')
  
  Normalized_myT1=read.table(file=paste('Kraken2_IBD_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  Normalized_myT2=read.table(file=paste('Metaphlan2_IBD_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  
  colnames(Normalized_myT1)=str_replace_all(colnames(Normalized_myT1),pattern="[[.]]",replacement="_")
  colnames(Normalized_myT2)=str_replace_all(colnames(Normalized_myT2),pattern="[[.]]",replacement="_")
  
  Normalized_Mean_myT1=log10(colMeans(10^Normalized_myT1[,4:ncol(Normalized_myT1)]-1)+1)
  Normalized_Mean_myT2=log10(colMeans(10^Normalized_myT2[,4:ncol(Normalized_myT2)]-1)+1)
  
  #Kruskal Wallis test for Kraken2
  pValuesDisease_Normalized_myT1=vector()
  Name_Normalized_myT1 <- vector()
  
  for(i in 4:(ncol(Normalized_myT1))){
    Disease=factor(Normalized_myT1$IBD)
    kruskal <- kruskal.test(formula=Normalized_myT1[,i]~Disease)
    pValuesDisease_Normalized_myT1[i-3] <- kruskal$p.value
    Name_Normalized_myT1[i-3] <- names(Normalized_myT1)[i]
  }
  FDR_Corrected_pValuesDisease_Normalized_myT1=p.adjust(pValuesDisease_Normalized_myT1, method = "BH")
  
  #Kruskal Wallis test for Metaphlan2
  pValuesDisease_Normalized_myT2=vector()
  Name_Normalized_myT2 <- vector()
  
  for(i in 4:(ncol(Normalized_myT2))){
    Disease=factor(Normalized_myT2$IBD)
    kruskal <- kruskal.test(formula=Normalized_myT2[,i]~Disease)
    pValuesDisease_Normalized_myT2[i-3] <- kruskal$p.value
    Name_Normalized_myT2[i-3] <- names(Normalized_myT2)[i]
  }
  FDR_Corrected_pValuesDisease_Normalized_myT2=p.adjust(pValuesDisease_Normalized_myT2, method = "BH")
  
  ordered_FDR_Normalized_myT1=FDR_Corrected_pValuesDisease_Normalized_myT1[order(-Mean_myT1)]
  ordered_FDR_Normalized_myT2=FDR_Corrected_pValuesDisease_Normalized_myT2[order(-Mean_myT2)]
  
  #Add three to ignore metadata
  ordered_Normalized_myT1=Normalized_myT1[,order(-Mean_myT1)+3]
  ordered_Normalized_myT2=Normalized_myT2[,order(-Mean_myT2)+3]
  
  Original_Order_Normalized_Mean_myT1=Normalized_Mean_myT1[order(-Mean_myT1)]
  Original_Order_Normalized_Mean_myT2=Normalized_Mean_myT2[order(-Mean_myT2)]
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/IBD/")
  
  pdf(paste("IBD_Raw_Normalized_Kraken2_Spearman_",input[k],".pdf",sep=""))
  plot(Original_Order_Normalized_Mean_myT1[2:length(Original_Order_Normalized_Mean_myT1)],highest,xlab = "Log Mean",
       ylab="Rho",main=paste("IBD Kraken2 Spearman test at ",input[k]," level",sep=""),ylim=c(-1,1),
       col= ifelse(ordered_FDR_Normalized_myT1[-c(1)]>0.05,"black","blue"),pch=16)
  dev.off()
  
  pdf(paste("IBD_Raw_Normalized_Metaphlan2_Spearman_",input[k],".pdf",sep=""))
  plot(Original_Order_Normalized_Mean_myT2[2:length(Original_Order_Normalized_Mean_myT2)],highest_2,xlab = "Log Mean",
       ylab="Rho",main=paste("IBD Metaphlan2 Spearman test at ",input[k]," level",sep=""),ylim=c(-1,1),
       col= ifelse(ordered_FDR_Normalized_myT2[-c(1)]>0.05,"black","blue"),pch=16)
  dev.off()
  
}
