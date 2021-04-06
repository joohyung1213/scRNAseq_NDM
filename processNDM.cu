#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/generate.h>

#include <omp.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"

using namespace std;
using namespace thrust;

/*
FILEPATH_EXAMPLE "D://CSN//D3.csv"
boxsize_EXAMPLE 0.1
sig_level_EXAMPLE 0.01
icdf 2.3236
*/

int gene_num, cell_num, boxsize;
float icdf;
double ro_temp1, ro_temp2;

vector <vector <float>> data_vector;

vector <vector<int>> greybox;
vector <int> greybox_line;

vector <vector<int>> NDM;


int my_mod_start1 = 0;
int my_mod_start2 = 0;
int my_mod1() {
	return (my_mod_start1++) / cell_num;
}
int my_mod2() {
	return (my_mod_start2++) / boxsize;
}

void ProcessNDM() {
	cout << "[Process] Start\n";

	thrust::host_vector <int> data_vector_sort_idx;
	thrust::host_vector <float> data_vector_sort;

	thrust::device_vector <int> data_vector_sort_idx_dev;
	thrust::device_vector <float> data_vector_sort_dev;

	thrust::device_vector<float> d_result2;
	thrust::host_vector<int> h_segments;
	thrust::device_vector<int> d_segments;

	thrust::device_vector<int> greybox_dev;

	thrust::host_vector<int> h_segments1;
	thrust::device_vector<int> d_segments2;

	thrust::host_vector<int> cooColIndex_host;
	thrust::host_vector<int> cooRowIndex_host;

	int cell_iter1, cell_iter2, gene_iter1, gene_iter2;

	int sort_partition = 1;
	int gene_per_partition = (int)(gene_num / sort_partition);
	int size_per_partition = (int)(gene_per_partition * cell_num);

	cudaError_t cudaStat;
	cublasStatus_t stat;
	cublasHandle_t handle;
	stat = cublasCreate(&handle);
	if (stat != CUBLAS_STATUS_SUCCESS) {
		printf("CUBLAS initialization failed\n");
	}

	const float alpha = 1.0;
	const float beta = 0.0;
	double ro;

	for (cell_iter1 = 0; cell_iter1 < cell_num; cell_iter1++) {
		cout << "CEll : " << cell_iter1 << endl;
		/*
		*
		*	[ Part1 ]
		*	In (GENE G1) calculate distance of (CELL C1 (Criteria) ) and (CELL C2)
		*	Get H Number of index of (CELL C2) which is closest with (CELL C1)
		*	GrayBox size(GENENUM, H) is (index of CELL) which is closest with (CELL C1) in each GENE
		*
		*	Each partition is Chunk of Columns
		*	Chunk of Columns(GENE) are copied to Device, then the Device sorts each Column(GENE)
		*	data_vector_sort, data_vector_sort_dev : Chunk of Columns
		*	data_vector_sort_idx, data_vector_sort_idx_dev : Indices of Cells in each Column
		*/

		for (int part_iter = 0; part_iter < sort_partition; part_iter++) {
			data_vector_sort.resize(size_per_partition);
			data_vector_sort_idx_dev.resize(size_per_partition);

			for (gene_iter1 = 0; gene_iter1 < gene_per_partition; gene_iter1++) {
				// Set the idx vector(device) in sequence
				thrust::sequence(data_vector_sort_idx_dev.begin() + (gene_iter1 * cell_num), data_vector_sort_idx_dev.begin() + ((gene_iter1 + 1) * cell_num));

				// Set the data vector(host) with distance
				int gene_idx = part_iter * gene_per_partition + gene_iter1;
				float cri_data = data_vector[gene_idx][cell_iter1];
#pragma omp parallel for
				for (cell_iter2 = 0; cell_iter2 < cell_num; cell_iter2++) {
					data_vector_sort[gene_iter1 * cell_num + cell_iter2] = abs(cri_data - data_vector[gene_idx][cell_iter2]);
				}
			}

			// Copy data vector(host) to data vector(device)
			data_vector_sort_dev = data_vector_sort;


			// Vectorized Batch sort. Sorts Chunk of Columns in column by column
			// Increasing amount of work per kernal call
			thrust::stable_sort_by_key(data_vector_sort_dev.begin(), data_vector_sort_dev.end(), data_vector_sort_idx_dev.begin());
			data_vector_sort_dev.clear();
			data_vector_sort_idx = data_vector_sort_idx_dev;
			data_vector_sort_idx_dev.clear();

			d_result2 = data_vector_sort;
			data_vector_sort.clear();

			my_mod_start1 = 0;
			h_segments.resize(size_per_partition);
			thrust::generate(h_segments.begin(), h_segments.end(), my_mod1);
			d_segments = h_segments;
			h_segments.clear();

			thrust::stable_sort_by_key(d_result2.begin(), d_result2.end(), d_segments.begin());
			d_result2.clear();

			data_vector_sort_idx_dev = data_vector_sort_idx;
			data_vector_sort_idx.clear();
			thrust::stable_sort_by_key(d_segments.begin(), d_segments.end(), data_vector_sort_idx_dev.begin());
			d_segments.clear();

			for (gene_iter1 = 0; gene_iter1 < gene_per_partition; gene_iter1++) {
				greybox_dev.insert(greybox_dev.end(), data_vector_sort_idx_dev.begin() + (gene_iter1 * cell_num), data_vector_sort_idx_dev.begin() + (gene_iter1 * cell_num) + boxsize);
			}
			data_vector_sort_idx_dev.clear();
		}

		my_mod_start2 = 0;
		h_segments1.resize(boxsize * gene_num);
		thrust::generate(h_segments1.begin(), h_segments1.end(), my_mod2);
		d_segments2 = h_segments1;
		h_segments1.clear();

		thrust::stable_sort_by_key(greybox_dev.begin(), greybox_dev.end(), d_segments2.begin());
		thrust::stable_sort_by_key(d_segments2.begin(), d_segments2.end(), greybox_dev.begin());

		cooColIndex_host = greybox_dev;
		cooRowIndex_host = d_segments2;
		greybox_dev.clear();
		d_segments2.clear();

		/*
		*
		*	[ Part2 ]
		*	Count number of cells in intersection of two GreyBox (nxyk)
		*	Matrix Multiplication of Greybox and Transpose of Greybox 
		*
		*/

		// Ready Data
		// greybox includes boolean values, but cuBlas doesn't supports boolean type.
		float* greybox;
		float* result;
		float* dev_greybox;
		float* dev_result;

		// Host
		greybox = (float*)malloc(gene_num * cell_num * sizeof(float));
		result = (float*)malloc(gene_num * gene_num * sizeof(float));
		std::fill_n(greybox, gene_num * cell_num, 0);
		std::fill_n(result, gene_num * gene_num, 0);

		for (gene_iter1 = 0; gene_iter1 < gene_num * boxsize; gene_iter1++) {
			greybox[cooColIndex_host[gene_iter1] * gene_num + cooRowIndex_host[gene_iter1]] = 1;
		}


		// Device
		cudaStat = cudaMalloc((void**)&dev_greybox, gene_num * cell_num * sizeof(float));
		if (cudaStat != cudaSuccess) {
			printf("dev_greybox : device memory allocation failed");
		}
		cudaStat = cudaMalloc((void**)&dev_result, gene_num * gene_num * sizeof(float));
		if (cudaStat != cudaSuccess) {
			printf("dev_result : device memory allocation failed");
		}

		stat = cublasSetMatrix(gene_num, cell_num, sizeof(float), greybox, gene_num, dev_greybox, gene_num);
		if (stat != CUBLAS_STATUS_SUCCESS) {
			printf("dev_greybox : data download failed");
		}
		stat = cublasSetMatrix(gene_num, gene_num, sizeof(float), result, gene_num, dev_result, gene_num);
		if (stat != CUBLAS_STATUS_SUCCESS) {
			printf("dev_result : data download failed");
		}

		cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, gene_num, gene_num, cell_num, &alpha, dev_greybox, gene_num, dev_greybox, gene_num, &beta, dev_result, gene_num);
		cudaDeviceSynchronize();


		stat = cublasGetMatrix(gene_num, gene_num, sizeof(float), dev_result, gene_num, result, gene_num);
		if (stat != CUBLAS_STATUS_SUCCESS) {
			printf("result : data upload failed");
		}

		for (gene_iter1 = 0; gene_iter1 < gene_num; gene_iter1++) {
			for (gene_iter2 = gene_iter1 + 1; gene_iter2 < gene_num; gene_iter2++) {
				ro = result[gene_iter1 * gene_num + gene_iter2] * ro_temp1 - ro_temp2;
				if (ro > icdf) {
					NDM[gene_iter1][cell_iter1]++;
					NDM[gene_iter2][cell_iter1]++;
				}
			}
		}

		cudaFree(dev_greybox);
		cudaFree(dev_result);

		free(greybox);
		free(result);

	}
	cublasDestroy(handle);
}