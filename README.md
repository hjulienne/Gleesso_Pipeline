# Gleesso Pipeline



We implemented the Gleesso pipeline to adapt the SPIEC-EASI pipeline to shotgun MGS abundances. From the Matrix of MGS abundance the pipeline goes through the following step:

* Select MGS that are present in at least 5% of samples and have a mean FPKM abundance of 10−7 across samples.
* Infer the graph with the SPIEC EASI pipeline.
* Infer the community structure with the walktrap algorithm from the igraph R package http: //igraph.org/r/.
* Output the graph structure to a Gephi https://gephi.org/ compatible format.
* Compute the community abundances by summing the abundance of MGS in the community.


# Installation

The dependency packages phyloseq and SPIEC-EASI (https://github.com/zdk123/SpiecEasi) need to be manually installed before installing Gleesso. The Gleesso package specifically depends on the 0.1.4 version of SPIEC-EASI.


The following code lines install both the dependencies and the gleesso package.

```r
install.packages('devtools')
library(devtools)

source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")

devtools::install_github("zdk123/SpiecEasi@v0.1.4")

library(SpiecEasi)
install_github("hjulienne/Gleesso_Pipeline/Gleesso")
```

# Usage example

Consult the Gleesso walkthrough notebook :

https://github.com/hjulienne/Gleesso_Pipeline/blob/master/Gleesso_walkthrough/Gleesso_walkthrough.ipynb
