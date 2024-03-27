# CASCC

CASCC, a clustering method designed to improve biological accuracy using gene co-expression features identified using an unsupervised adaptive attractor algorithm. 

The users can refer to our manuscript for a detailed description of the modeling and applications. 

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
