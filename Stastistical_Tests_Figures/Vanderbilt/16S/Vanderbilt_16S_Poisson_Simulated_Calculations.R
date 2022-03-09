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
  
  boolean=TRUE
  boolean=append(boolean,colSums(myT1[,2:ncol(myT1)])>0)
  myT1=myT1[,boolean]
  
  
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
  
  boolean=TRUE
  boolean=append(boolean,colSums(myT2[,2:ncol(myT2)])>0)
  myT2=myT2[,boolean]
  
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
  
  Mean_myT1=colMeans(myT1[,2:ncol(myT1)])
  Mean_myT2=colMeans(myT2[,2:ncol(myT2)])
  
  ordered_Mean_myT1=Mean_myT1[order(-Mean_myT1)]
  ordered_Mean_myT2=Mean_myT2[order(-Mean_myT2)]
  
  ordered_myT1=myT1[order(-Mean_myT1)+1]
  ordered_myT2=myT2[order(-Mean_myT2)+1]
  
  ordered_myT1_colnames=names(ordered_Mean_myT1)
  ordered_myT2_colnames=names(ordered_Mean_myT2)
  
  Normalized_myT1=RDP_Normalized
  Normalized_myT2=Qiime_Normalized
  
  Normalized_Mean_myT1=log10(colMeans(10^Normalized_myT1[,2:ncol(Normalized_myT1)]-1)+1)
  Normalized_Mean_myT2=log10(colMeans(10^Normalized_myT2[,2:ncol(Normalized_myT2)]-1)+1)
  
  Original_Order_Normalized_Mean_myT1=Normalized_Mean_myT1[order(-Mean_myT1)]
  Original_Order_Normalized_Mean_myT2=Normalized_Mean_myT2[order(-Mean_myT2)]
  
  #RDP
  highest=rep(0,length(ordered_Mean_myT1)-1)
  for (i in 2:length(ordered_Mean_myT1)){
    if (ordered_Mean_myT1[i]>0.01){
      for (j in 2:length(ordered_Mean_myT1)){
        if (i!=j && ordered_Mean_myT1[j]>0.01){
      temp=cor.test(x=ordered_myT1[,i],y=ordered_myT1[,j],method="spearman",alternative="two.sided")$estimate
      if (abs(temp)>highest[i-1]){
        highest[i-1]=temp
                  } } }
      }   }

  #Qiime
  highest_2=rep(0,length(ordered_Mean_myT2)-1)
  for (i in 2:length(ordered_Mean_myT2)){
    if (ordered_Mean_myT2[i]>0.01){
      for (j in 2:length(ordered_Mean_myT2)){
        if (i!=j && ordered_Mean_myT2[j]>0.01){
      temp=cor.test(x=ordered_myT2[,i],y=ordered_myT2[,j],method="spearman",alternative="two.sided")$estimate
      if (abs(temp)>highest_2[i-1]){
        highest_2[i-1]=temp
                    } } } 
    }   }
  
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
  
  ordered_FDR_Normalized_myT1=FDR_Corrected_pValuesSampleType_RDP[order(-Mean_myT1)]
  ordered_FDR_Normalized_myT2=FDR_Corrected_pValuesSampleType_Qiime[order(-Mean_myT2)]
  
  #added 1 to ignore metadata
  
  ordered_Normalized_myT1=Normalized_myT1[,order(-Mean_myT1)+1]
  ordered_Normalized_myT2=Normalized_myT2[,order(-Mean_myT2)+1]
  
  #Grab the top ten Taxa for both datasets
  a<-10
  probs_1<-vector(length=a)
  SumTop_1<-sum((10^Original_Order_Normalized_Mean_myT1[1:a])-1)
  
  probs_2<-vector(length=a)
  SumTop_2<-sum((10^Original_Order_Normalized_Mean_myT2[1:a])-1)
  
  for( i in 1:a)  {
    #Use distribution of top ten taxa as probabilities of phantom taxa for RDP
    a_Sum_1<-(10^Original_Order_Normalized_Mean_myT1[i])-1
    probs_1[i]=a_Sum_1/SumTop_1
    #Use distribution of top ten taxa as probabilities of phantom taxa for Qiime
    a_Sum_2<-(10^Original_Order_Normalized_Mean_myT2[i])-1
    probs_2[i]=a_Sum_2/SumTop_2
  }
  
  #RDP
  
  simMeans_RDP<-vector()
  simCors_RDP<-vector()
  
  for( i in 1:(ncol(ordered_Normalized_myT1)-a)){
    #random Error Rate that is less than .01
    backgroundErrorRate = runif(1)/ 1000
    
    if( backgroundErrorRate < 0) {
      backgroundErrorRate  = 0
    }
    
    # choose by relative abundance for first top ten
    colIndex <- sample(1:a,1,prob=probs_1)
    
    dataSimCol<-vector(length=nrow(ordered_Normalized_myT1))
    
    for( j in 1:nrow(ordered_Normalized_myT1)) {
      dataSimCol[j]<-rpois(1,(ceiling((10^ordered_Normalized_myT1[j,colIndex]))-1)*backgroundErrorRate )
    }
    
    simMeans_RDP[i]<-log10(mean(dataSimCol)+1) 
    
    if( is.nan(simMeans_RDP[i])|simMeans_RDP[i]<0 ) {
      simMeans_RDP[i]=0
    }
    simCors_RDP[i]<-cor( dataSimCol,ordered_Normalized_myT1[,colIndex],method="spearman" )
  }
  
  #Qiime
  
  simMeans_Qiime<-vector()
  simCors_Qiime<-vector()
  
  for( i in 1:(ncol(ordered_Normalized_myT2)-a)){
    #random Error Rate that is less than .001
    backgroundErrorRate = runif(1)/ 1000
    
    if( backgroundErrorRate < 0) {
      backgroundErrorRate  = 0
    }
    
    # choose by relative abundance for first top ten
    colIndex <- sample(1:a,1,prob=probs_2)
    
    dataSimCol<-vector(length=nrow(ordered_Normalized_myT2))
    
    for( j in 1:nrow(ordered_Normalized_myT2)) {
      dataSimCol[j]<-rpois(1,(ceiling((10^ordered_Normalized_myT2[j,colIndex]))-1)*backgroundErrorRate )
    }
    
    simMeans_Qiime[i]<-log10(mean(dataSimCol)+1) 
    
    if( is.nan(simMeans_Qiime[i])|simMeans_Qiime[i]<0 ) {
      simMeans_Qiime[i]=0
    }
    simCors_Qiime[i]<-cor( dataSimCol,ordered_Normalized_myT2[,colIndex],method="spearman" )
  }
  
  setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/Vanderbilt/16S/")
  
  pdf(paste("Vanderbilt_16S_Poisson_Simulated_Spearman_for_RDP_at_",input[k],".pdf",sep=""))
  plot(Original_Order_Normalized_Mean_myT1[2:length(Original_Order_Normalized_Mean_myT1)],highest,xlab = "Log Mean",
       ylab="Rho",main=paste("Vanderbilt RDP Spearman vs Simulated at ",input[k],sep=""),ylim=c(-1,1),xlim=c(0,6.5),
       col= ifelse(ordered_FDR_Normalized_myT1[-c(1)]>0.05,"black","orange"),pch=16)
  legend("bottom",inset=c(0,.01),c("Insignificant","Significant for RDP","Simulated Data"), 
         pch = c(16, 16,16),cex=0.8,col=c("black","orange","green") )
  points(simMeans_RDP, simCors_RDP,col="green")
  dev.off()
  
  
  pdf(paste("Vanderbilt_16S_Poisson_Simulated_Spearman_for_Qiime_at_",input[k],".pdf",sep=""))
  plot(Original_Order_Normalized_Mean_myT2[2:length(Original_Order_Normalized_Mean_myT2)],highest_2,xlab = "Log Mean",
       ylab="Rho",main=paste("Vanderbilt Qiime Spearman vs Simulated at ",input[k],sep=""),ylim=c(-1,1),xlim=c(0,6.5),
       col= ifelse(ordered_FDR_Normalized_myT2[-c(1)]>0.05,"black","purple"),pch=16)
  legend("bottom",inset=c(0,.01),c("Insignificant","Significant for Qiime","Simulated Data"),
         pch = c(16, 16,16),cex=0.8,col=c("black","purple","green") )
  points(simMeans_Qiime, simCors_Qiime,col="green")
  dev.off()
  
}  
  