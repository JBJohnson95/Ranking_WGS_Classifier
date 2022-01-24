rm(list=ls())
library(data.table)
library(ecodist)
library(vegan)
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
  
  #Define taxa that are present
  
  myT1_colnames=colnames(myT1[,2:ncol(myT1)])
  myT2_colnames=colnames(myT2[,2:ncol(myT2)])
  
  #Average Calculations based upon Tables
  
  Mean_myT1=log10(colMeans(10^RDP_Normalized[,2:ncol(RDP_Normalized)]-1)+1)
  Mean_myT2=log10(colMeans(10^Qiime_Normalized[,2:ncol(Qiime_Normalized)]-1)+1)
  
  ordered_Mean_myT1=Mean_myT1[order(-Mean_myT1)]
  ordered_Mean_myT2=Mean_myT2[order(-Mean_myT2)]
  
  ordered_myT1=RDP_Normalized[order(-Mean_myT1)+1]
  ordered_myT2=Qiime_Normalized[order(-Mean_myT2)+1]
  
  ordered_myT1_colnames=names(ordered_Mean_myT1)
  ordered_myT2_colnames=names(ordered_Mean_myT2)
  
  #Kruskal Wallis test for RDP
  pValuesSampleType_RDP=vector()
  Name_RDP <- vector()
  
  for(i in 2:(ncol(RDP_Normalized))){
    kruskal <- kruskal.test(formula=RDP_Normalized[,i]~Sample_Type)
    pValuesSampleType_RDP[i-1] <- kruskal$p.value
    Name_RDP[i-1] <- names(RDP_Normalized)[i]
  }
  FDR_Corrected_pValuesSampleType_RDP=p.adjust(pValuesSampleType_RDP, method = "BH")
  
  #Kruskal Wallis test for Qiime
  pValuesSampleType_Qiime=vector()
  Name_Qiime <- vector()
  
  for(i in 2:(ncol(Table2))){
    kruskal <- kruskal.test(formula=Table2[,i]~Sample_Type)
    pValuesSampleType_Qiime[i-1] <- kruskal$p.value
    Name_Qiime[i-1] <- names(Qiime_Normalized)[i]
  }
  FDR_Corrected_pValuesSampleType_Qiime=p.adjust(pValuesSampleType_Qiime, method = "BH")
  
  
  
  #WGS 
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/Vanderbilt_Forward_Only/")  
  
  myT3=read.table(file=paste('Kraken2_Vanderbilt_Forward_Only_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  myT4=read.table(file=paste('Metaphlan2_Vanderbilt_Forward_Only_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  
  colnames(myT3)=str_replace_all(colnames(myT3),pattern="[[.]]",replacement="_")
  colnames(myT4)=str_replace_all(colnames(myT4),pattern="[[.]]",replacement="_")
  
  overlap_names=intersect(colnames(myT3[,6:ncol(myT3)]),colnames(myT4[,6:ncol(myT4)]))
  
  myT3_colnames=colnames(myT3[,6:ncol(myT3)])
  myT4_colnames=colnames(myT4[,6:ncol(myT4)])
  
  Mean_myT3=log10(colMeans(10^myT3[,6:ncol(myT3)]-1)+1)
  Mean_myT4=log10(colMeans(10^myT4[,6:ncol(myT4)]-1)+1)
  
  #Kruskal Wallis test for Kraken2
  pValuesSampleType_myT3=vector()
  Name_myT3 <- vector()
  
  for(i in 6:(ncol(myT3))){
    SampleType=factor(myT3$sample_type)
    kruskal <- kruskal.test(formula=myT3[,i]~SampleType)
    pValuesSampleType_myT3[i-5] <- kruskal$p.value
    Name_myT3[i-5] <- names(myT3)[i]
  }
  FDR_Corrected_pValuesSampleType_myT3=p.adjust(pValuesSampleType_myT3, method = "BH")
  
  #Kruskal Wallis test for Metaphlan2
  pValuesSampleType_myT4=vector()
  Name_myT4 <- vector()
  
  for(i in 6:(ncol(myT4))){
    SampleType=factor(myT4$sample_type)
    kruskal <- kruskal.test(formula=myT4[,i]~SampleType)
    pValuesSampleType_myT4[i-5] <- kruskal$p.value
    Name_myT4[i-5] <- names(myT4)[i]
  }
  FDR_Corrected_pValuesSampleType_myT4=p.adjust(pValuesSampleType_myT4, method = "BH")
  
 
  #RDP is Table 1, Qiime is Table 2, Kraken is Table 3, Metaphlan is Table 4
  #This merges the different taxanomic databases together
  taxa_names=unique(c(myT1_colnames,myT2_colnames,myT3_colnames,myT4_colnames))
  RDP_means=rep(0,length(taxa_names))
  Qiime_means=rep(0,length(taxa_names))
  Kraken_means=rep(0,length(taxa_names))
  Metaphlan_means=rep(0,length(taxa_names))
  Kraken_Adjusted_Pvalue=rep(1,length(taxa_names))
  Metaphlan_Adjusted_Pvalue=rep(1,length(taxa_names))
  Qiime_Adjusted_Pvalue=rep(1,length(taxa_names))
  RDP_Adjusted_Pvalue=rep(1,length(taxa_names))
  
  for (x in 1:length(taxa_names)){
      if((taxa_names[x] %in% myT1_colnames)){
      RDP_means[x]=Mean_myT1[myT1_colnames==taxa_names[x]]
      RDP_Adjusted_Pvalue[x]=FDR_Corrected_pValuesSampleType_RDP[Name_RDP==taxa_names[x]]
      }
      if((taxa_names[x] %in% myT2_colnames)){
      Qiime_means[x]=Mean_myT2[myT2_colnames==taxa_names[x]]
      Qiime_Adjusted_Pvalue[x]=FDR_Corrected_pValuesSampleType_Qiime[Name_Qiime==taxa_names[x]]
      }
    if((taxa_names[x] %in% myT3_colnames)){
      Kraken_means[x]=Mean_myT3[myT3_colnames==taxa_names[x]]
      Kraken_Adjusted_Pvalue[x]=FDR_Corrected_pValuesSampleType_myT3[Name_myT3==taxa_names[x]]
    }
    if((taxa_names[x] %in% myT4_colnames)){
      Metaphlan_means[x]=Mean_myT4[myT4_colnames==taxa_names[x]]
      Metaphlan_Adjusted_Pvalue[x]=FDR_Corrected_pValuesSampleType_myT4[Name_myT4==taxa_names[x]]
    }
  }
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/Vanderbilt/16S/")
  
  R2=summary(lm(Kraken_means~RDP_means))$r.squared
  R2=format(round(R2,digits=3))
  
  pdf(paste("Kruskal_Wallis_Vanderbilt_Kraken2_vs_RDP_Average_Means_",input[k],"_16S.pdf",sep=""))
  plot(RDP_means,Kraken_means,xlim =c(0,8),ylim=c(0,8),xlab = "Log Normalized RDP",ylab="Log Normalized Kraken2",
       main = paste("Vanderbilt Kraken2 vs RDP Averge Value at ",input[k]," level\nR2 is ",R2,sep=""),
       col= ifelse(Kraken_Adjusted_Pvalue>0.05&RDP_Adjusted_Pvalue>0.05,"black",
                   ifelse(Kraken_Adjusted_Pvalue<=0.05&RDP_Adjusted_Pvalue>0.05,"blue",
                          ifelse(Kraken_Adjusted_Pvalue>0.05&RDP_Adjusted_Pvalue<=0.05,"darkorange","purple"))),pch=16)
  
  legend("topleft",inset=c(.1,0),c("Insignificant after correction for Both ",
                                   "FDR Significant only Kraken2","FDR Significant only RDP","FDR Significant for Both"),
         pch = c(16, 16),cex=0.8,
         col=c("black","blue","darkorange","purple")   )
  abline(0,1)
  
  labeled_index=order(-(RDP_means+Kraken_means))[1:5]
  text(Kraken_means[labeled_index]~RDP_means[labeled_index],labels=taxa_names[labeled_index],
       pos=4,cex=0.6,font=2)
  
  dev.off()
  
  
  
  R2=summary(lm(Metaphlan_means~RDP_means))$r.squared
  R2=format(round(R2,digits=3))
  
  pdf(paste("Kruskal_Wallis_Vanderbilt_Metaphlan2_vs_RDP_Average_Means_",input[k],"_16S.pdf",sep=""))
  plot(RDP_means,Metaphlan_means,xlim =c(0,8),ylim=c(0,8),xlab = "Log Normalized RDP",ylab="Log Normalized Metaphlan2",
       main = paste("Vanderbilt Metaphlan2 vs RDP Averge Value at ",input[k]," level\nR2 is ",sep=""),
       col= ifelse(Metaphlan_Adjusted_Pvalue>0.05&RDP_Adjusted_Pvalue>0.05,"black",
                   ifelse(Metaphlan_Adjusted_Pvalue<=0.05&RDP_Adjusted_Pvalue>0.05,"red",
                          ifelse(Metaphlan_Adjusted_Pvalue>0.05&RDP_Adjusted_Pvalue<=0.05,"darkorange","purple"))),pch=16)
  
  legend("topleft",inset=c(.1,0),c("Insignificant after correction for Both ",
                                   "FDR Significant only Metaphlan2","FDR Significant only RDP","FDR Significant for Both"),
         pch = c(16, 16),cex=0.8,
         col=c("black","red","darkorange","purple")   )
  abline(0,1)
  
  labeled_index=order(-(RDP_means+Metaphlan_means))[1:5]
  text(Metaphlan_means[labeled_index]~RDP_means[labeled_index],labels=taxa_names[labeled_index],
       pos=4,cex=0.6,font=2)
  
  dev.off()
  
  
  
  R2=summary(lm(Kraken_means~Qiime_means))$r.squared
  R2=format(round(R2,digits=3))
  
  pdf(paste("Kruskal_Wallis_Vanderbilt_Kraken2_vs_Qiime_Average_Means_",input[k],"_16S.pdf",sep=""))
  plot(Qiime_means,Kraken_means,xlim =c(0,8),ylim=c(0,8),xlab = "Log Normalized Qiime",ylab="Log Normalized Kraken2",
       main = paste("Vanderbilt Kraken2 vs Qiime Averge Value at ",input[k]," level\nR2 is ",R2,sep=""),
       col= ifelse(Qiime_Adjusted_Pvalue>0.05&Kraken_Adjusted_Pvalue>0.05,"black",
                   ifelse(Qiime_Adjusted_Pvalue<=0.05&Kraken_Adjusted_Pvalue>0.05,"darkgreen",
                          ifelse(Qiime_Adjusted_Pvalue>0.05&Kraken_Adjusted_Pvalue<=0.05,"blue","purple"))),pch=16)
  
  legend("topleft",inset=c(.1,0),c("Insignificant after correction for Both ",
                                   "FDR Significant only Qiime","FDR Significant only Kraken2","FDR Significant for Both"),
         pch = c(16, 16),cex=0.8,
         col=c("black","darkgreen","blue","purple")   )
  abline(0,1)
  
  labeled_index=order(-(Qiime_means+Kraken_means))[1:5]
  text(Kraken_means[labeled_index]~Qiime_means[labeled_index],labels=taxa_names[labeled_index],
       pos=4,cex=0.6,font=2)
  
  dev.off()
  
  
  
  R2=summary(lm(Metaphlan_means~Qiime_means))$r.squared
  R2=format(round(R2,digits=3))
  
  pdf(paste("Kruskal_Wallis_Vanderbilt_Metaphlan2_vs_Qiime_Average_Means_",input[k],"_16S.pdf",sep=""))
  plot(Qiime_means,Metaphlan_means,xlim =c(0,8),ylim=c(0,8),xlab = "Log Normalized Qiime",ylab="Log Normalized Metaphlan2",
       main = paste("Vanderbilt Metaphlan2 vs Qiime Averge Value at ",input[k]," level\nR2 is ",R2,sep=""),
       col= ifelse(Qiime_Adjusted_Pvalue>0.05&Metaphlan_Adjusted_Pvalue>0.05,"black",
                   ifelse(Qiime_Adjusted_Pvalue<=0.05&Metaphlan_Adjusted_Pvalue>0.05,"darkgreen",
                          ifelse(Qiime_Adjusted_Pvalue>0.05&Metaphlan_Adjusted_Pvalue<=0.05,"red","purple"))),pch=16)
  
  legend("topleft",inset=c(.1,0),c("Insignificant after correction for Both ",
                                   "FDR Significant only Qiime","FDR Significant only Metaphlan2","FDR Significant for Both"),
         pch = c(16, 16),cex=0.8,
         col=c("black","darkgreen","red","purple")   )
  abline(0,1)
  
  labeled_index=order(-(Qiime_means+Metaphlan_means))[1:5]
  text(Metaphlan_means[labeled_index]~Qiime_means[labeled_index],labels=taxa_names[labeled_index],
       pos=4,cex=0.6,font=2)
  
  dev.off()

  
}
