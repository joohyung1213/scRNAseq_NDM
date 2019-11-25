# scRNA Cell Specific Network Degree Matrix (with GPU, OpenMP)

Constructing Network Degree Matrix (NDM) of Cell-specific Network (CSN) constructed by single-cell RNA sequencing data

Original method from D. Hao, et al. (2019) *Cell-specific network constructed by single-cell RNA sequencing data." Nucleic Acids Research*., **47**, e62-e62

Construction of network degree matrix with reduced main memory usage. Leverages GPU and multi-threading by OpenMP for speedup.

Original NDM constructing consumes memory size of Input File and Output File.

## NDM_v1.cu

**NDM_v1.cu** reads Input File by buffer size and construct NDM of read data. Then reads next buffer of input file and constructs NDM again.
It requires memory size of only buffer size and size of output file (same as input file size).

## NDM_v2.cu

**NDM_v2.cu** reads Input File by buffer size same as NDM_v1.cu. In addition, writes the result values in the temporary result files for reduced usage of main memory.
It requires memory size of only buffer size for reading input file.

NDM_v1.cu is recommended if system has enough main memory to store output file.

NDM_v1.cu is recommended if system has less main memory to store output file.

## Input File

CSV file of Gene Expression Matrix without column name and row name (only number values)

Each row : Gene

Each column : Cell

## Output File

CSV file of Network Degree Matrix
Output file has same size with Input file

Each row : Gene

Each column : Cell

## Input Parameters

**File Path** : Input File Path (i.e. D://Data//data.csv)

**Buffer Size** : Integer value in Gigabytes.

**boxsize** : Size of neighborhood. (GreyBox)

**alpha** : Significant Level

## System Requirements

NVIDIA CUDA available **GPU**
Installation of **CUDA runtime API** and **CUDA thrust API**
Installation of **OpenMP API**

## Compile and Execute

**Windows** : Execute by Visual Studio CUDA Project
**Linux** : $ nvcc NDM_v1.cu -std=c++14 -Xcompiler=-fopenmp -O2

-O2 is compile optimization option of g++ in linux. ([https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html)) ([https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#ptxas-options-opt-level](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#ptxas-options-opt-level))

## System Requirements

NVIDIA CUDA available GPU
Installation of CUDA runtime API and CUDA thrust API
Installation of OpenMP API

## Performance

Performance varies by number of CPU core, main memory size, GPU, and I/O speed of secondary memory (SSD, HDD).