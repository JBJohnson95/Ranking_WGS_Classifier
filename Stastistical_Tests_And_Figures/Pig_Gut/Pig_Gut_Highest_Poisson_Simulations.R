rm(list=ls())
library(data.table)
library(ecodist)
library(stringr)

input=c("phylum","class","order","family","genus")

setwd('/Users/James/Desktop/Dissertation/metadata/')
in_meta=read.table('PIG_GUT_Metadata_NCBI_paired.tsv',header = TRUE,sep="\t")
metadata_columns=c("Sample_Names","Alias","Geography")
colnames(in_meta)[which(names(in_meta) == "Sample")] <- "Sample_Names"
meta=in_meta[metadata_columns]

meta_col=3

for (k in 1:length(input)){

setwd(dir = "/Users/James/Desktop/Dissertation/Farnaz_Kraken2_Unclassified_Reads/Tables/Pig_Gut/")  

#Kraken2 
myT1=read.csv(file=paste("WGS_Kraken2_Pig_Gut_2021Apr14_taxaCount_",input[k],'.tsv',sep=""),sep="\t",header=TRUE)
colnames(myT1)[which(names(myT1) == "Sample")] <- "Sample_Names"
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

overlap_names=intersect(colnames(myT1[,(meta_col+1):ncol(myT1)]),colnames(myT2[,(meta_col+1):ncol(myT2)]))

setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/Normalized_Tables/Pig_Gut/")  

Normalized_myT1=read.table(file=paste('Kraken2_Pig_Gut_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")
Normalized_myT2=read.table(file=paste('Metaphlan2_Pig_Gut_',input[k],'_Normalized.csv',sep=""),header=TRUE,sep = ",")

colnames(Normalized_myT1)=str_replace_all(colnames(Normalized_myT1),pattern="[[.]]",replacement="_")
colnames(Normalized_myT2)=str_replace_all(colnames(Normalized_myT2),pattern="[[.]]",replacement="_")

#Convert log normalized means to regular space to average and then back to log base 10 space

Normalized_Mean_myT1=log10(colMeans(10^Normalized_myT1[,(meta_col+1):ncol(Normalized_myT1)]-1)+1)
Normalized_Mean_myT2=log10(colMeans(10^Normalized_myT2[,(meta_col+1):ncol(Normalized_myT2)]-1)+1)

#Kruskal Wallis test for Kraken2
pValuesGeography_Normalized_myT1=vector()
Name_Normalized_myT1 <- vector()

for(i in (meta_col+1):(ncol(Normalized_myT1))){
  Geography=factor(Normalized_myT1$Geography)
  kruskal <- kruskal.test(formula=Normalized_myT1[,i]~Geography)
  pValuesGeography_Normalized_myT1[i-meta_col] <- kruskal$p.value
  Name_Normalized_myT1[i-meta_col] <- names(Normalized_myT1)[i]
}
FDR_Corrected_pValuesGeography_Normalized_myT1=p.adjust(pValuesGeography_Normalized_myT1, method = "BH")

#Kruskal Wallis test for Metaphlan2
pValuesGeography_Normalized_myT2=vector()
Name_Normalized_myT2 <- vector()

for(i in (meta_col+1):(ncol(Normalized_myT2))){
  Geography=factor(Normalized_myT2$Geography)
  kruskal <- kruskal.test(formula=Normalized_myT2[,i]~Geography)
  pValuesGeography_Normalized_myT2[i-meta_col] <- kruskal$p.value
  Name_Normalized_myT2[i-meta_col] <- names(Normalized_myT2)[i]
}

FDR_Corrected_pValuesGeography_Normalized_myT2=p.adjust(pValuesGeography_Normalized_myT2, method = "BH")

ordered_FDR_Normalized_myT1=FDR_Corrected_pValuesGeography_Normalized_myT1[order(-Normalized_Mean_myT1)]
ordered_FDR_Normalized_myT2=FDR_Corrected_pValuesGeography_Normalized_myT2[order(-Normalized_Mean_myT2)]

#Add two to ignore metadata
ordered_Normalized_myT1=Normalized_myT1[,order(-Normalized_Mean_myT1)+meta_col]
ordered_Normalized_myT2=Normalized_myT2[,order(-Normalized_Mean_myT2)+meta_col]

ordered_Normalized_Mean_myT1=Normalized_Mean_myT1[order(-Normalized_Mean_myT1)]
ordered_Normalized_Mean_myT2=Normalized_Mean_myT2[order(-Normalized_Mean_myT2)]

ordered_myT1=myT1[,order(-Normalized_Mean_myT1)+meta_col]
ordered_myT2=myT2[,order(-Normalized_Mean_myT2)+meta_col]

### Kraken2 data is labelled as _1

maxCors_1 <- vector(length=ncol(ordered_myT1))
sums_1 <- apply(ordered_myT1[,meta_col:ncol(ordered_myT1)], 2, sum)
numCols_1 <- ncol(ordered_myT1)

#Twice perfect negative correlation should be statisically unlikely or impossible thus the ideal default
for( i in 1:numCols_1)
{
  maxCors_1[i] <- -2
}

#Ever decreasing lower average mean Taxa are compared against higher mean taxa starting at second most abundant
for( i in (numCols_1:2))
{
  stopVal <- min(i-1)
  
  for( j in stopVal:1){
    corVal <- cor(ordered_myT1[,i],ordered_myT1[,j],method="spearman" )
    
    if( ! is.na(corVal) ){	
      maxCors_1[i] <- max(  maxCors_1[i], corVal  )
    }
  }
}
maxCors_1[1] <- 1

a1 <- 10
# take the a number of top taxa and get the average depth after being ordered from largest to smallest
a1_averageDepth <- ordered_Normalized_Mean_myT1[1:a1]
aSum_1 <- 10^(a1_averageDepth)-1
sumTop_1 <- sum(aSum_1-1)
#Distribution of a number of top taxa asigned as probabiltiies 
probs_1 = aSum_1/sumTop_1

simMeans_1 <- vector()
simCors_1 <- vector()

for(i in (a1+1):(numCols_1))
{
  backgroundErrorRate = runif(1)/1000
  
  if( backgroundErrorRate < 0) {
    backgroundErrorRate  = 0
  }
  
  # choose by relative abundance
  colIndex <- sample(1:a1,1,prob=probs_1)
  dataSimCol <- vector(length=nrow(ordered_Normalized_myT1))
  
  for( j in 1:nrow(ordered_Normalized_myT1)) {
    dataSimCol[j] <- rbinom(1,ceiling(10^ordered_Normalized_myT1[j,colIndex])-1,backgroundErrorRate)
  }
  
  simMeans_1[i] <- log10(mean(dataSimCol)+1) 
  
  if( is.nan(simMeans_1[i]) | simMeans_1[i] < 0 ) {
    simMeans_1[i] = 0
  }
  
  simCors_1[i] <- cor(dataSimCol,10^ordered_Normalized_myT1[,colIndex],method="spearman" )
}


#Metaphlan2 data is labelled as _2

maxCors_2 <- vector(length=ncol(ordered_myT2))
sums_2 <- apply(ordered_myT2[,meta_col:ncol(ordered_myT2)], 2, sum)
numCols_2 <- ncol(ordered_myT2)

#Twice perfect negative correlation should be statisically unlikely or impossible thus the ideal default
for( i in 1:numCols_2 )
{
  maxCors_2[i] <- -2
}

#Ever decreasing lower average mean Taxa are compared against higher mean taxa starting at second most abundant
for( i in (numCols_2:2))
{
  stopVal <- min(i-1)
  
  for( j in stopVal:1){
    corVal <- cor(ordered_myT2[,i],ordered_myT2[,j],method="spearman" )
    
    if( ! is.na(corVal) ){	
      maxCors_2[i] <- max(  maxCors_2[i], corVal  )
    }
  }
}
maxCors_2[1] <- 1

a2 <- 10
# take the a number of top taxa and get the average depth after being ordered from largest to smallest
a2_averageDepth <- ordered_Normalized_Mean_myT2[1:a2]
aSum_2 <- 10^(a2_averageDepth)-1
sumTop_2 <- sum(aSum_2-1)
#Distribution of a number of top taxa asigned as probabiltiies 
probs_2 = aSum_2/sumTop_2

simMeans_2 <- vector()
simCors_2 <- vector()

for(i in (a2+1):(numCols_2))
{
  backgroundErrorRate = runif(1)/1000
  
  if( backgroundErrorRate < 0) {
    backgroundErrorRate  = 0
  }
  
  # choose by relative abundance
  colIndex <- sample(1:a2,1,prob=probs_2)
  dataSimCol <- vector(length=nrow(ordered_Normalized_myT2))
  
  for( j in 1:nrow(ordered_Normalized_myT2)) {
    dataSimCol[j] <- rbinom(1,ceiling(10^ordered_Normalized_myT2[j,colIndex])-1,backgroundErrorRate)
  }
  
  simMeans_2[i] <- log10(mean(dataSimCol)+1) 
  
  if( is.nan(simMeans_2[i]) | simMeans_2[i] < 0 ) {
    simMeans_2[i] = 0
  }
  
  simCors_2[i] <- cor(dataSimCol,10^ordered_Normalized_myT2[,colIndex],method="spearman" )
}

maxCors_1[1] <- NA
maxCors_2[1] <- NA

setwd(dir = "/Users/James/Desktop/Dissertation/Final_Figures/Pig_Gut/")

pdf(paste("Pig_Gut_Highest_Poisson_Simulated_Spearman_for_Kraken2_at_",input[k],".pdf",sep=""))
plot(ordered_Normalized_Mean_myT1, maxCors_1,xlab = "Log Mean",
     ylab="Rho",main=paste("Pig Gut Kraken2 Spearman vs Simulated at ",input[k],sep=""),ylim=c(-1,1),xlim=c(0,6.5),
     col= ifelse(ordered_FDR_Normalized_myT1[-c(1)]>0.05,"black","blue"),pch=16)
legend("bottom",c("Insignificant","Significant for Kraken2","Simulated Data"),
       pch = c(16, 16,16),cex=0.8,col=c("black","blue","green") )
points(simMeans_1, simCors_1,col="green")
dev.off()

pdf(paste("Pig_Gut_Highest_Poisson_Simulated_Spearman_for_Metaphlan2_at_",input[k],".pdf",sep=""))
plot(ordered_Normalized_Mean_myT2, maxCors_2,xlab = "Log Mean",
     ylab="Rho",main=paste("Pig Gut Metaphlan2 Spearman vs Simulated at ",input[k],sep=""),ylim=c(-1,1),xlim=c(0,6.5),
     col= ifelse(ordered_FDR_Normalized_myT2[-c(1)]>0.05,"black","red"),pch=16)
legend("bottom",c("Insignificant","Significant for Metaphlan2","Simulated Data"),
       pch = c(16, 16,16),cex=0.8,col=c("black","red","green") )
points(simMeans_2, simCors_2,col="green")
dev.off()

}
