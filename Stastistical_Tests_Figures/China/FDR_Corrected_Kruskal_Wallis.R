rm(list=ls())
library(data.table)
library(ecodist)
library(vegan)
library(stringr)

input=c("phylum","class","order","family","genus")

for (k in 1:length(input)){
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/China/")  
  
  myT1=read.table(file=paste('Kraken2_China_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  myT2=read.table(file=paste('Metaphlan2_China_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
  
  colnames(myT1)=str_replace_all(colnames(myT1),pattern="[[.]]",replacement="_")
  colnames(myT2)=str_replace_all(colnames(myT2),pattern="[[.]]",replacement="_")
  
  overlap_names=intersect(colnames(myT1[,3:ncol(myT1)]),colnames(myT2[,3:ncol(myT2)]))
  
  myT1_colnames=colnames(myT1[,3:ncol(myT1)])
  myT2_colnames=colnames(myT2[,3:ncol(myT2)])
  
  Mean_myT1=log10(colMeans(10^myT1[,3:ncol(myT1)]-1)+1)
  Mean_myT2=log10(colMeans(10^myT2[,3:ncol(myT2)]-1)+1)
  
  #Kruskal Wallis test for Kraken2
  pValuesruralurban_myT1=vector()
  Name_myT1 <- vector()
  
  for(i in 3:(ncol(myT1))){
    ruralurban=factor(myT1$ruralurban)
    kruskal <- kruskal.test(formula=myT1[,i]~ruralurban)
    pValuesruralurban_myT1[i-2] <- kruskal$p.value
    Name_myT1[i-2] <- names(myT1)[i]
  }
  FDR_Corrected_pValuesruralurban_myT1=p.adjust(pValuesruralurban_myT1, method = "BH")
  
  #Kruskal Wallis test for Metaphlan2
  pValuesruralurban_myT2=vector()
  Name_myT2 <- vector()
  
  for(i in 3:(ncol(myT2))){
    ruralurban=factor(myT2$ruralurban)
    kruskal <- kruskal.test(formula=myT2[,i]~ruralurban)
    pValuesruralurban_myT2[i-2] <- kruskal$p.value
    Name_myT2[i-2] <- names(myT2)[i]
  }
  FDR_Corrected_pValuesruralurban_myT2=p.adjust(pValuesruralurban_myT2, method = "BH")
  
  #
  taxa_names=unique(c(myT1_colnames,myT2_colnames))
  Kraken_means=rep(0,length(taxa_names))
  Metaphlan_means=rep(0,length(taxa_names))
  Metaphlan2_Adjusted_Pvalue=rep(1,length(taxa_names))
  Kraken2_Adjusted_Pvalue=rep(1,length(taxa_names))
  
  for (x in 1:length(taxa_names)){
    if ((taxa_names[x] %in% myT1_colnames)&&(taxa_names[x] %in% myT2_colnames)){
      Kraken_means[x]=Mean_myT1[myT1_colnames==taxa_names[x]]
      Metaphlan_means[x]=Mean_myT2[myT2_colnames==taxa_names[x]]
      Kraken2_Adjusted_Pvalue[x]=FDR_Corrected_pValuesruralurban_myT1[Name_myT1==taxa_names[x]]
      Metaphlan2_Adjusted_Pvalue[x]=FDR_Corrected_pValuesruralurban_myT2[Name_myT2==taxa_names[x]]
    }
    else if((taxa_names[x] %in% myT1_colnames)){
      Kraken_means[x]=Mean_myT1[myT1_colnames==taxa_names[x]]
      Kraken2_Adjusted_Pvalue[x]=FDR_Corrected_pValuesruralurban_myT1[Name_myT1==taxa_names[x]]
    }
    else if((taxa_names[x] %in% myT2_colnames)){
      Metaphlan_means[x]=Mean_myT2[myT2_colnames==taxa_names[x]]
      Metaphlan2_Adjusted_Pvalue[x]=FDR_Corrected_pValuesruralurban_myT2[Name_myT2==taxa_names[x]]
    }
  }
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/China/")
  
  R2=summary(lm(Metaphlan_means~Kraken_means))$r.squared
  R2=format(round(R2,digits=3))
  
  pdf(paste("Kruskal_Wallis_China_Metaphlan_vs_Kraken_Average_Means_",input[k],".pdf",sep=""))
  plot(Kraken_means,Metaphlan_means,xlim =c(0,8),ylim=c(0,8),xlab = "Log Normalized Kraken2",ylab="Log Normalized Metaphlan2",
       main = paste("China Metaphlan2 vs Kraken2 Averge Value at ",input[k]," level\nR2 is ",R2,sep=""),
       col= ifelse(Metaphlan2_Adjusted_Pvalue>0.05&Kraken2_Adjusted_Pvalue>0.05,"black",
                   ifelse(Metaphlan2_Adjusted_Pvalue<=0.05&Kraken2_Adjusted_Pvalue>0.05,"red",
                          ifelse(Metaphlan2_Adjusted_Pvalue>0.05&Kraken2_Adjusted_Pvalue<=0.05,"blue","purple"))),pch=16)
  
  legend("topleft",inset=c(.1,0),c("Insignificant after correction for Both ",
                                   "FDR Significant only Metaphlan2","FDR Significant only Kraken2","FDR Significant for Both"),
         pch = c(16, 16),cex=0.8,
         col=c("black","red","blue","purple")   )
  abline(0,1)
  
  labeled_index=order(-(Metaphlan_means+Kraken_means))[1:5]
  text(Metaphlan_means[labeled_index]~Kraken_means[labeled_index],labels=taxa_names[labeled_index],
       pos=4,cex=0.6,font=2)
  
 dev.off()
  
}