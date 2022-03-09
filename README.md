# Transparency and Reproducibility

This repository contains:

All the R code and input files to complete the stastical tests employed in the preprint detailed below. 
Please Note that all directory paths will need to be to updated based upon your local machine along with the locations of the input files and how you want to store the output files.

# **Manuscript**

The manuscript for this project is currently available through              bioarhive placeholder

Title: Systematic classification error profoundly impacts inference in high-depth WGS datasets

Authors:
James Johnson, Shan Sun, Anthony A. Fodor

doi:          placeholder

This article is a preprint and has not been peer-reviewed.

# **Abstract**

There is little consensus in the literature as to which approach for classification of Whole Genome Shotgun (WGS) sequences is best. In this paper, we examine two of the most popular algorithms, Kraken2 and Metaphlan2 utilizing four publicly available datasets. As expected from previous literature, we found that Kraken2 reports more overall taxa while Metaphlan2 reports fewer taxa while classifying fewer overall reads. To our surprise, however, Kraken 2 reported not only more taxa but many more taxa that were significantly associated with metadata. This implies that either Kraken2 is more sensitive to taxa that are biologically relevant and are simply missed by Metaphlan2, or that Kraken2’s classification errors are generated in such a way to impact inference. To discriminate between these two possibilities, we compared Spearman correlations coefficients of each taxa against all other taxa from the same dataset. We found that Kraken2, but not Metaphlan2, showed a consistent pattern of classifying low abundance taxa that generated high correlation coefficients with other taxa. Neither Metaphlan2, nor 16S sequences that were available for two of our four datasets, showed this pattern.  Simple simulations based on a variable Poisson error rate sampled from the uniform distribution with an average error rate of 0.0005 showed strikingly strong concordance with the observed correlation patterns from Kraken2.  Our results strongly suggest that Kraken2 consistently misclassifies high abundance taxa into the same erroneous low abundance taxa creating “phantom” taxa have a similar pattern of inference as the high abundance source. Because of the large sequencing depths of modern WGS cohorts, these “phantom” taxa will appear statistically significant in statistical models even with a low overall rate of classification error from Kraken.  Our simulations suggest that this can occur with average error rates as low as 1 in 2,000 reads.   These data suggest a novel metric for evaluating classifier accuracy and suggest that the pattern of classification errors should be considered in addition to overall classification error rate since consistent classification errors have a more profound impact on inference compared to classification errors that do not always result in assignment to the same erroneous taxa.  This work highlights fundamental questions on how classifiers function and interact with large sequencing depth and statistical models that still need to be resolved for WGS, especially if correlation coefficients between taxa are to be used to build covariance networks. Our work also suggests that despite its limitations, 16S sequencing may still be useful as neither of the two most popular 16S classifiers showed these patterns of inflated correlation coefficients between taxa.
