# scRNA Cell Specific Network Degree Matrix (with GPU, OpenMP)

Constructing Network Degree Matrix (NDM) of Cell-specific Network (CSN) constructed by single-cell RNA sequencing data

Original method from D. Hao, et al. (2019) *Cell-specific network constructed by single-cell RNA sequencing data." Nucleic Acids Research*., **47**, e62-e62

Construction of network degree matrix with reduced main memory usage. Leverages GPU and multi-threading by OpenMP for speedup.

## Input File

CSV file of **Gene Expression Matrix** without column name and row name (only number values)

Each row : Gene

Each column : Cell

## Output File

CSV file of **Network Degree Matrix**

Output file has same size with Input file

Each row : Gene

Each column : Cell

## Input Parameters

**File Path** : Input File Path (i.e. D://Data//data.csv)

**boxsize** : Size of neighborhood. (GreyBox)

**alpha** : Significant Level

## System Requirements

NVIDIA CUDA available **GPU**

Installation of **CUDA runtime API** and **CUDA thrust API**

Installation of **OpenMP API**

## Compile and Execute

**Windows** : Execute by Visual Studio CUDA Project  
**Linux** : $ nvcc NDM.cu -std=c++14 -Xcompiler=-fopenmp -O2

-O2 is compile optimization option of g++ in linux. ([https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html)) ([https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#ptxas-options-opt-level](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#ptxas-options-opt-level))


