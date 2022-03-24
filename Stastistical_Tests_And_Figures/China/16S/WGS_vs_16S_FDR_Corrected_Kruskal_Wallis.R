rm(list=ls())
library(data.table)
library(ecodist)
library(vegan)
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
  
  #Define taxa that are present
  
  myT1_colnames=colnames(myT1[,2:ncol(myT1)])
  myT2_colnames=colnames(myT2[,2:ncol(myT2)])
  
  #Average Calculations
  
  Mean_myT1=log10(colMeans(10^RDP_Normalized[,2:ncol(RDP_Normalized)]-1)+1)
  Mean_myT2=log10(colMeans(10^Qiime_Normalized[,2:ncol(Qiime_Normalized)]-1)+1)
  
  #Kruskal Wallis test for RDP
  pValuesruralurban_RDP_Normalized=vector()
  Name_RDP_Normalized <- vector()
  
  for(i in 2:(ncol(RDP_Normalized))){
    ruralurban=factor(RDP_ruralurban)
    kruskal <- kruskal.test(formula=RDP_Normalized[,i]~ruralurban)
    pValuesruralurban_RDP_Normalized[i-1] <- kruskal$p.value
    Name_RDP_Normalized[i-1] <- names(RDP_Normalized)[i]
  }
  FDR_Corrected_pValuesruralurban_RDP_Normalized=p.adjust(pValuesruralurban_RDP_Normalized, method = "BH")
  
  #Kruskal Wallis test for QIIME
  pValuesruralurban_Qiime_Normalized=vector()
  Name_Qiime_Normalized <- vector()
  
  for(i in 2:(ncol(Qiime_Normalized))){
    ruralurban=factor(Qiime_ruralurban)
    kruskal <- kruskal.test(formula=Qiime_Normalized[,i]~ruralurban)
    pValuesruralurban_Qiime_Normalized[i-1] <- kruskal$p.value
    Name_Qiime_Normalized[i-1] <- names(Qiime_Normalized)[i]
  }
  FDR_Corrected_pValuesruralurban_Qiime_Normalized=p.adjust(pValuesruralurban_Qiime_Normalized, method = "BH")
  
  #WGS
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/China/")  
  
  myT3=read.table(file=paste('Kraken2_China_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  myT4=read.table(file=paste('Metaphlan2_China_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  
  colnames(myT3)=str_replace_all(colnames(myT3),pattern="[[.]]",replacement="_")
  colnames(myT4)=str_replace_all(colnames(myT4),pattern="[[.]]",replacement="_")
  
  overlap_names=intersect(colnames(myT3[,3:ncol(myT3)]),colnames(myT4[,3:ncol(myT4)]))
  
  myT3_colnames=colnames(myT3[,3:ncol(myT3)])
  myT4_colnames=colnames(myT4[,3:ncol(myT4)])
  
  Mean_myT3=log10(colMeans(10^myT3[,3:ncol(myT3)]-1)+1)
  Mean_myT4=log10(colMeans(10^myT4[,3:ncol(myT4)]-1)+1)
  
  #Kruskal Wallis test for Kraken2
  pValuesruralurban_myT3=vector()
  Name_myT3 <- vector()
  
  for(i in 3:(ncol(myT3))){
    ruralurban=factor(myT3$ruralurban)
    kruskal <- kruskal.test(formula=myT3[,i]~ruralurban)
    pValuesruralurban_myT3[i-2] <- kruskal$p.value
    Name_myT3[i-2] <- names(myT3)[i]
  }
  FDR_Corrected_pValuesruralurban_myT3=p.adjust(pValuesruralurban_myT3, method = "BH")
  
  #Kruskal Wallis test for Metaphlan2
  pValuesruralurban_myT4=vector()
  Name_myT4 <- vector()
  
  for(i in 3:(ncol(myT4))){
    ruralurban=factor(myT4$ruralurban)
    kruskal <- kruskal.test(formula=myT4[,i]~ruralurban)
    pValuesruralurban_myT4[i-2] <- kruskal$p.value
    Name_myT4[i-2] <- names(myT4)[i]
  }
  FDR_Corrected_pValuesruralurban_myT4=p.adjust(pValuesruralurban_myT4, method = "BH")
  
  
  #
  taxa_names=unique(c(myT1_colnames,myT2_colnames,myT3_colnames,myT4_colnames))
  
  Kraken_means=rep(0,length(taxa_names))
  Metaphlan_means=rep(0,length(taxa_names))
  Metaphlan_Adjusted_Pvalue=rep(1,length(taxa_names))
  Kraken_Adjusted_Pvalue=rep(1,length(taxa_names))

  RDP_means=rep(0,length(taxa_names))
  Qiime_means=rep(0,length(taxa_names))
  Qiime_Adjusted_Pvalue=rep(1,length(taxa_names))
  RDP_Adjusted_Pvalue=rep(1,length(taxa_names))
  
  for (x in 1:length(taxa_names)){
    if((taxa_names[x] %in% myT1_colnames)){
      RDP_means[x]=Mean_myT1[myT1_colnames==taxa_names[x]]
      RDP_Adjusted_Pvalue[x]=FDR_Corrected_pValuesruralurban_RDP_Normalized[Name_RDP_Normalized==taxa_names[x]]
    }
    if((taxa_names[x] %in% myT2_colnames)){
      Qiime_means[x]=Mean_myT2[myT2_colnames==taxa_names[x]]
      Qiime_Adjusted_Pvalue[x]=FDR_Corrected_pValuesruralurban_Qiime_Normalized[Name_Qiime_Normalized==taxa_names[x]]
    }
    if((taxa_names[x] %in% myT3_colnames)){
        Kraken_means[x]=Mean_myT3[myT3_colnames==taxa_names[x]]
        Kraken_Adjusted_Pvalue[x]=FDR_Corrected_pValuesruralurban_myT3[Name_myT3==taxa_names[x]]
      }
    if((taxa_names[x] %in% myT4_colnames)){
        Metaphlan_means[x]=Mean_myT4[myT4_colnames==taxa_names[x]]
        Metaphlan_Adjusted_Pvalue[x]=FDR_Corrected_pValuesruralurban_myT4[Name_myT4==taxa_names[x]]
      }
    }
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/China/16S/")
  
  R2=summary(lm(Kraken_means~RDP_means))$r.squared
  R2=format(round(R2,digits=3))
  
  pdf(paste("Kruskal_Wallis_China_Kraken2_vs_RDP_Average_Means_",input[k],"_16S.pdf",sep=""))
  plot(RDP_means,Kraken_means,xlim =c(0,8),ylim=c(0,8),xlab = "Log Normalized RDP",ylab="Log Normalized Kraken2",
       main = paste("China Kraken2 vs RDP Averge Value at ",input[k]," level\nR2 is ",R2,sep=""),
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
  
  pdf(paste("Kruskal_Wallis_China_Metaphlan2_vs_RDP_Average_Means_",input[k],"_16S.pdf",sep=""))
  plot(RDP_means,Metaphlan_means,xlim =c(0,8),ylim=c(0,8),xlab = "Log Normalized RDP",ylab="Log Normalized Metaphlan2",
       main = paste("China Metaphlan2 vs RDP Averge Value at ",input[k]," level\nR2 is ",R2,sep=""),
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
  
  pdf(paste("Kruskal_Wallis_China_Kraken2_vs_QIIME_Average_Means_",input[k],"_16S.pdf",sep=""))
  plot(Qiime_means,Kraken_means,xlim =c(0,8),ylim=c(0,8),xlab = "Log Normalized QIIME",ylab="Log Normalized Kraken2",
       main = paste("China Kraken2 vs QIIME Averge Value at ",input[k]," level\nR2 is ",R2,sep=""),
       col= ifelse(Qiime_Adjusted_Pvalue>0.05&Kraken_Adjusted_Pvalue>0.05,"black",
                   ifelse(Qiime_Adjusted_Pvalue<=0.05&Kraken_Adjusted_Pvalue>0.05,"pink",
                          ifelse(Qiime_Adjusted_Pvalue>0.05&Kraken_Adjusted_Pvalue<=0.05,"blue","purple"))),pch=16)
  
  legend("topleft",inset=c(.1,0),c("Insignificant after correction for Both ",
                                   "FDR Significant only QIIME","FDR Significant only Kraken2","FDR Significant for Both"),
         pch = c(16, 16),cex=0.8,
         col=c("black","pink","blue","purple")   )
  abline(0,1)
  
  labeled_index=order(-(Qiime_means+Kraken_means))[1:5]
  text(Kraken_means[labeled_index]~Qiime_means[labeled_index],labels=taxa_names[labeled_index],
       pos=4,cex=0.6,font=2)
  
  dev.off()
  
  
  
  R2=summary(lm(Metaphlan_means~Qiime_means))$r.squared
  R2=format(round(R2,digits=3))
  
  pdf(paste("Kruskal_Wallis_China_Metaphlan2_vs_QIIME_Average_Means_",input[k],"_16S.pdf",sep=""))
  plot(Qiime_means,Metaphlan_means,xlim =c(0,8),ylim=c(0,8),xlab = "Log Normalized QIIME",ylab="Log Normalized Metaphlan2",
       main = paste("China Metaphlan2 vs QIIME Averge Value at ",input[k]," level\nR2 is ",R2,sep=""),
       col= ifelse(Qiime_Adjusted_Pvalue>0.05&Metaphlan_Adjusted_Pvalue>0.05,"black",
                   ifelse(Qiime_Adjusted_Pvalue<=0.05&Metaphlan_Adjusted_Pvalue>0.05,"pink",
                          ifelse(Qiime_Adjusted_Pvalue>0.05&Metaphlan_Adjusted_Pvalue<=0.05,"red","purple"))),pch=16)
  
  legend("topleft",inset=c(.1,0),c("Insignificant after correction for Both ",
                                   "FDR Significant only QIIME","FDR Significant only Metaphlan2","FDR Significant for Both"),
         pch = c(16, 16),cex=0.8,
         col=c("black","pink","red","purple")   )
  abline(0,1)
  
  labeled_index=order(-(Qiime_means+Metaphlan_means))[1:5]
  text(Metaphlan_means[labeled_index]~Qiime_means[labeled_index],labels=taxa_names[labeled_index],
       pos=4,cex=0.6,font=2)
  
  dev.off()
}