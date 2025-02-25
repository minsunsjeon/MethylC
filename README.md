# MethylC

R package for meanshift clusering of CpG sites

Author: Minsun S. Jeon, Casey Pei

Description
-------

MethylC is a bioinformatics package designed for identifying and analyzing CpG regions in the genome. We employ the meanshift clustering algorithm to identify CpG-rich regions in an unsupervised manner. 
Unlike traditional methods that impose fixed thresholds, meanshift clustering leverages the natural density distribution of CpG sites across the genome, allowing for adaptive identification of CpG clusters. 
This approach ensures that regions with a high density of CpGs are detected without relying on arbitrary cutoffs, making it more flexible and capable of capturing non-canonical CpG clusters. 
By identifying these regions based on their intrinsic clustering patterns, our method provides a more objective and statistically robust alternative for defining CpG-enriched regions.

**Key Features:**
- Efficient detection of CpG regions in genomic data
- Customizable thresholds and parameters for cluster identification
- Support for genome-wide analysis
- Exportable results for downstream epigenetic studies

**Dependencies:** data.table, dplyr, parallel, Rcpp, RcppArmadillo

Installation
--------
`devtools::install_github("minsunsjeon/MethylC")`


Available functions:
---------
|Code| Function |
|:-|:-|
|applyMeanshift|Applies the meanshift algorithm to the methylation dataset|
|formatIslands|Takes the result of the meanshift algorithm and make it into a format that includes the list of islands and single sites|
