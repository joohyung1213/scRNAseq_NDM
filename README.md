# scRNA

Constructing Network Degree Matrix of Cell-specific Network constructed by single-cell RNA sequencing data

Original method from D. Hao, et al. (2019) *Cell-specific network constructed by single-cell RNA sequencing data." Nucleic Acids Research*., **47**, e62-e62

Construction of network degree matrix with reduced main memory usage.

I**nput File** : CSV file of Gene Expression Matrix
each row : Gene
each column : Cell

**Output File** : Network Degree Matrix
each row : Gene
each column : Cell

**System Inputs**
File Path : input file path. (i.e. D://Data//data.csv)
Buffer Size : integer Gigabytes.
boxsize : Size of neighborhood. (GreyBox)
alpha : Significant level.

**System Requirements**
NVIDIA CUDA available GPU
Installation of CUDA runtime API and CUDA thrust API
Installation of OpenMP API

**Compile and Execute**
Windows
Linux : $ nvcc [NDM.cu](http://ndm.cu) -std=c++14 -Xcompiler=-fopenmp
