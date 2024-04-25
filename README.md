# CASCC: A co-expression assisted single-cell RNA-seq data clustering method


CASCC, a clustering method designed to improve biological accuracy using gene co-expression features identified using an unsupervised adaptive attractor algorithm. 

The users can refer to our published manuscript, available at [this link](https://doi.org/10.1093/bioinformatics/btae283), for a detailed description of the method. 

<!-- overview -->

# Installation
```R
require(devtools)
remotes::install_version("Seurat", "4.3.0") # dependency
install_github("weiyi-bitw/cafr") # dependency
install_github("LingyiC/CASCC")
CASCC::checkDependencies() 
```
# Tutorial
Vignette: https://rpubs.com/Lingyi/CASCC

# Benchmarking
https://github.com/LingyiC/CASCC_benchmark

<!-- # Citation 
If you use CASCC in your research, please consider citing:
 -->
