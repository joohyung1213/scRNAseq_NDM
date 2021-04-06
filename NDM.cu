#include <iostream>
#include <string>

#include "preproc.h"
#include "processNDM.cuh"

#include <omp.h>
#include <cuda_runtime.h>

int main() {
	std::cout << "=== CSN NDM Implementation (cuBLAS) ===\n";
	std::cout << "=== DATA READ (1) THREAD, Partition Sort, Reduced Memory Sort ===\n";
	std::cout << "OpenMP Thread: " << omp_get_max_threads() << "\n";
	cudaSetDevice(0);

	std::string input_file_path;
	std::string file_directory;
	double boxsize_ratio = 0.1;
	double sig_level = 0.01; // alpha

	input_file_path = FileOpen();
	file_directory = InitializeFunction(input_file_path, boxsize_ratio, sig_level);
	DataRead(input_file_path);
	ProcessNDM();
	RecordResult(file_directory);

	return 0;
}
