# **Warning**

After observing these finding you might become BFs with the Huttenhower Lab and become a biobaker

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

There are many different approaches to classifying WSG sequences in the literature but little consensus as to which approach is best. In this paper, we benchmark WGS classifiers by utilizing four publicly available datasets with significant associations with metadata of varying effect sizes. As expected from previous literature, in comparing the two most popular classifiers Kraken2 and Metaphlan2, we found that Kraken2 reports more overall taxa while Metaphlan2 reports fewer taxa based on classifying fewer reads. To our surprise, however, Kraken 2 reported not only more taxa but many more taxa that were significantly associated with metadata. This suggests that either Kraken 2 is more sensitive to taxa that are biologically relevant and are simply missed by Metaphlan2, or that Kraken 2’s classification errors are generated in such a way to impact inference. To discriminate between these two possibilities, we compared Spearman correlations coefficients of each taxa with more abundant taxa from the same dataset. We found that Kraken 2, but not Metaphlan 2, showed a consistent pattern of classifying low abundance taxa that generated high correlation coefficients with higher abundance taxa. Neither Metaphlan2, nor 16S sequences that were available for two of our four datasets, showed this pattern. Our results strongly suggest that Kraken 2 consistently misclassifies high abundance taxa into the same erroneous low abundance taxa creating “phantom” taxa have a similar pattern of inference as the high abundance source. Because of the large sequencing depths of modern WGS cohorts, these “phantom” taxa will appear statistically significant in statistical models even with a low overall rate of classification error from Kraken. These data suggest a novel metric for evaluating classifier accuracy and highlights fundamental questions on how classifiers work and interact with large sequencing depth and statistical models that still need to be resolved for WGS, especially if correlation coefficients between taxa are to be used to build covariance networks. Our work also suggests that despite its limitations, 16S sequencing may still be useful as neither of the two most popular 16S classifiers showed these patterns of inflated correlation coefficients between taxa.
