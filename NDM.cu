#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <thread>
#include <condition_variable>

#include <vector>
#include <string>
#include <unordered_map>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>

#include <omp.h>

/*
FILEPATH_EXAMPLE "D://CSVData//data.csv"
BOXSIZE_EXAMPLE 0.1
alpha_EXAMPLE 0.01
icdf 2.3236
*/

#define GigaByte_in_Byte 1073741824

using namespace std;
using namespace thrust;

int buffer_size_in_GB = 0;

int gene_num = 0;
int cell_num = 0;
int boxsize = 0;
float icdf = 0;
double ro_temp1 = 0;
double ro_temp2 = 0;

int criterion_gene_idx = 0;
int process_gene_idx = 0;
int buffer_size_gene_num = 0;

bool reading = true;

vector <vector <float>> data_vector;
vector <float> data_vector_line;

vector <vector<int>> greybox;
vector <int> greybox_line;

vector <vector<int>> NDM;

unordered_map <int, int> greybox_idx_map;

double RationalApproximation(double t);
double NormalCDFInverse(double p);
string SplitFileDirectory(string filepath, char sep);
void ConcatResultFiles();
string InitializeFunction(string filepath, double boxsize_ratio, double alpha);
void DataRead(string filepath, condition_variable* cv, mutex* m);
void ProcessNDM(string filepath, condition_variable* cv, mutex* m);



double RationalApproximation(double t) {
	/*
	Made by John D. Cook, PhD, President
	https://www.johndcook.com/blog/cpp_phi_inverse/
	*/

	// Abramowitz and Stegun formula 26.2.23.
	// The absolute value of the error should be less than 4.5 e-4.
	double c[] = { 2.515517, 0.802853, 0.010328 };
	double d[] = { 1.432788, 0.189269, 0.001308 };
	return t - ((c[2] * t + c[1])*t + c[0]) /
		(((d[2] * t + d[1])*t + d[0])*t + 1.0);
}
double NormalCDFInverse(double p) {
	/*
	Made by John D. Cook, PhD, President
	https://www.johndcook.com/blog/cpp_phi_inverse/
	*/

	if (p <= 0.0 || p >= 1.0) {
		std::stringstream os;
		os << "Invalid input argument (" << p
			<< "); must be larger than 0 but less than 1.";
		throw std::invalid_argument(os.str());
	}

	// See article above for explanation of this section.
	if (p < 0.5) {
		// F^-1(p) = - G^-1(p)
		return -RationalApproximation(sqrt(-2.0*log(p)));
	}
	else {
		// F^-1(p) = G^-1(1-p)
		return RationalApproximation(sqrt(-2.0*log(1 - p)));
	}
}
string SplitFileDirectory(string filepath, char sep) {
	vector<string> paths;
	stringstream stream(filepath);
	string path_chunk;
	string directory;

	while (getline(stream, path_chunk, sep)) {
		paths.push_back(path_chunk);
	}
	for (int i = 0; i < paths.size() - 1; i++) {
		directory = directory + paths[i] + "/";
	}
	return directory;
}
void ConcatResultFiles(string file_directory) {
	int gene_iter;
	string result_file_path = file_directory + "NDM_Result.csv";
	ofstream result_file;
	result_file.open(result_file_path);

	if (result_file.fail())
		cout << strerror(errno) << endl;
	else if (result_file.is_open()) {
		for (gene_iter = 0; gene_iter < gene_num; gene_iter++) {
			string one_line;
			string temp_result_file_path = file_directory + "Result_Temp_" + to_string(gene_iter) + ".csv";
			ifstream temp_result_file;

			temp_result_file.open(temp_result_file_path);

			if (temp_result_file.fail())
				cout << strerror(errno) << endl;
			else if (temp_result_file.is_open()) {
				getline(temp_result_file, one_line);
				result_file << one_line << endl;
				temp_result_file.close();

				// remove concatenated temp result file
				vector<char> filepath_v(temp_result_file_path.begin(), temp_result_file_path.end());
				filepath_v.push_back('\0');
				char* charArrayPtr = &filepath_v[0];
				remove(charArrayPtr);
			}

		}

	}
	result_file.close();
}
string InitializeFunction(string filepath, double boxsize_ratio, double alpha) {
	int i = 0;
	string file_directory;
	ifstream file(filepath);

	if (file.fail()) cout << "Fail to Read File" << endl;
	else if (file.is_open()) {
		cout << "Initializing..." << endl;

		/*
		*
		*	[ Scan Data File Once to Get (Row, Column) Size ]
		*	[ Initialize Global Variables ]
		*
		*/
		// Count Column Number of input file
		string one_line;
		getline(file, one_line);
		stringstream line_stream(one_line);
		string one_data;
		while (getline(line_stream, one_data, ',')) {
			stringstream convertor(one_data);
			cell_num++;
		}
		gene_num++;

		// Count Row Number of input file
		while (getline(file, one_line)) {
			gene_num++;
		}

		// Initialize Global Variables
		boxsize = round(boxsize_ratio * cell_num) + 1;
		icdf = -NormalCDFInverse(alpha);

		ro_temp1 = cell_num * sqrt((double)(cell_num - 1)) / ((boxsize - 1) * (cell_num - (boxsize - 1)));
		ro_temp2 = sqrt((double)(cell_num - 1)) * (boxsize - 1) / (cell_num - (boxsize - 1));

		buffer_size_gene_num = ceil(GigaByte_in_Byte / (cell_num * sizeof(double)) * buffer_size_in_GB);
		buffer_size_gene_num = buffer_size_gene_num < gene_num ? buffer_size_gene_num : gene_num;

		vector <int> NDM_line(cell_num, 0);
		for (i = 0; i < gene_num; i++)NDM.push_back(NDM_line);

	}

	// Get the directory of input file to make result file
	file_directory = SplitFileDirectory(filepath, '/');

	cout << "Creating Result Files..." << endl;
	// Initialize temporary result files
#pragma omp parallel for schedule(guided)
	for (i = 0; i < gene_num; i++) {
		string file_path = file_directory + "Result_Temp_" + to_string(i) + ".csv";
		ofstream file;
		
		file.open(file_path, ios::out);
		if (file.fail()) cout << "Output Create File Fail (" << strerror(errno) << ")" << endl;
		else if (file.is_open()) {
			file << 0;
			for (int j = 1; j < cell_num; j++) {
				file << "," << 0;
			}
			file.close();
		}
	}
	
	cout << "Number of Gene : " << gene_num << endl;
	cout << "Number of Cell : " << cell_num << endl;
	cout << "GreyBox Size: " << boxsize << endl;
	cout << "icdf (alpha) : " << icdf << endl;
	cout << "Reading Buffer Size : " << buffer_size_gene_num << " rows(genes)\n" << endl;
	
	return file_directory;

}

void DataRead(string filepath, condition_variable* cv, mutex* m) {
	cout << "[DataRead] Thread Start" << endl;
	int i = 0;
	int gene_iter = 0;
	int buffer_counter = 0;

	ifstream file(filepath);
	string one_line;
	double one_data;

	unique_lock<mutex> lock(*m);



	if (file.fail()) cout << "Fail to Read File" << endl;
	else if (file.is_open()) {
		cout << "[DataRead] File Readed\n";

		// Wait for Process Thread Starts
		std::this_thread::sleep_for(3s);

		cout << "[DataRead] File Read Start\n";

		// Start Reading Data
		for (gene_iter = 0; gene_iter < gene_num; gene_iter++) {
			// Initialize File Read Pointer
			file.clear();
			file.istream::seekg(0, ios::beg);

			process_gene_idx = 0;
			buffer_counter = 0;
			data_vector_line.clear();
			data_vector.clear();

			/*
			*
			*	[ Set Criterion GENE and Processed GENE ]
			*
			*
			*/

			// Reach to Criterion GENE
			for (i = 0; i < criterion_gene_idx; i++) {
				getline(file, one_line);
			}

			process_gene_idx = criterion_gene_idx + 1;

			cout << "Criterion GENE :" << criterion_gene_idx << endl;

			// Read Criterion GENE
			getline(file, one_line);
			stringstream one_line_stream(one_line);
			string one_data_string;

			while (getline(one_line_stream, one_data_string, ',')) {
				stringstream one_data_stream(one_data_string);
				one_data_stream >> one_data;
				data_vector_line.push_back(one_data);
			}

			data_vector.push_back(data_vector_line);
			data_vector_line.clear();


			/*
			*
			*	[ Read Data which will processed much as buffer size ]
			*
			*
			*/

			// Read until Buffer is FUll or END of DATA
			while (getline(file, one_line)) {
				stringstream one_line_stream(one_line);
				string one_data_string;

				while (getline(one_line_stream, one_data_string, ',')) {
					stringstream one_data_stream(one_data_string);
					one_data_stream >> one_data;
					data_vector_line.push_back(one_data);
				}

				data_vector.push_back(data_vector_line);
				data_vector_line.clear();

				buffer_counter++;

				// if Buffer if FULL
				if (buffer_counter >= buffer_size_gene_num) {
					cv->notify_one();
					cv->wait(lock);

					process_gene_idx += buffer_counter;
					buffer_counter = 0;
					data_vector.erase(data_vector.begin() + 1, data_vector.end());
				}

			}
			// if Reading File is done
			if (data_vector.size() > 1) {
				cv->notify_one();
				cv->wait(lock);
			}
			criterion_gene_idx++;

		}

		data_vector.clear();
		reading = false;
		cv->notify_one();

	}

}

void ProcessNDM(string filepath, condition_variable* cv, mutex* m) {
	cout << "[ProcessNDM] Thread Start\n";
	int cell_iter1, cell_iter2, gene_iter1, gene_iter2, boxsize_iter, nxyk, criterion_gene_NDM_count;
	double ro;

	cout << "[ProcessNDM] CSN READY\n";
	unique_lock<mutex> lock(*m);

	cout << "[ProcessNDM] CSN Start\n";

	while (reading) {

		int process_Column_Size = data_vector.size();

		greybox_line.resize(boxsize);

		thrust::host_vector <int> data_vector_sort_idx;
		thrust::host_vector <double> data_vector_sort(cell_num * process_Column_Size);

		thrust::device_vector <int> data_vector_sort_idx_dev(cell_num * process_Column_Size);
		thrust::device_vector <double> data_vector_sort_dev;

		unordered_map <int, int>::iterator it;
		
		for (cell_iter1 = 0; cell_iter1 < cell_num; cell_iter1++) {
			/*
			*
			*	[ Part1 ]
			*	In (GENE gene_iter1) calculate distance of (CELL cell_iter1 (Criteria) ) and (CELL cell_iter2)
			*	Get indices of (CELL gell_iter2)s which is closest with (CELL cell_iter1).
			*
			*/
			data_vector_sort_idx.clear();
			data_vector_sort_dev.clear();

			// Set the idx vector(device) in sequence
			for (gene_iter1 = 0; gene_iter1 < process_Column_Size; gene_iter1++) {
				thrust::sequence(data_vector_sort_idx_dev.begin() + (gene_iter1 * cell_num), data_vector_sort_idx_dev.begin() + ((gene_iter1 + 1) * cell_num));

				// Set the data vector(host) with distance
				for (cell_iter2 = 0; cell_iter2 < cell_num; cell_iter2++) {
					data_vector_sort[gene_iter1 * cell_num + cell_iter2] = abs(data_vector[gene_iter1][cell_iter1] - data_vector[gene_iter1][cell_iter2]);
				}
			}

			// Copy data vector(host) to data vector(device)
			data_vector_sort_dev = data_vector_sort;

			// Sort data vector(device) and idx vector(device)
			for (gene_iter1 = 0; gene_iter1 < process_Column_Size; gene_iter1++) {
				thrust::sort_by_key(data_vector_sort_dev.begin() + (gene_iter1 * cell_num), data_vector_sort_dev.begin() + ((gene_iter1 + 1) * cell_num), data_vector_sort_idx_dev.begin() + (gene_iter1 * cell_num));
			}

			// Copy result idx vector(device) to idx vector(host)
			data_vector_sort_idx = data_vector_sort_idx_dev;

			// Get H Number of idx of closest C2
			for (gene_iter1 = 0; gene_iter1 < process_Column_Size; gene_iter1++) {
				for (boxsize_iter = 0; boxsize_iter < boxsize; boxsize_iter++) {
					greybox_line[boxsize_iter] = data_vector_sort_idx[(gene_iter1 * cell_num) + boxsize_iter];
				}
				greybox.push_back(greybox_line);
			}


			/*
			*
			*	[ Part2 ]
			*	Count number of cells in intersection of two GreyBox (nxyk)
			*
			*/
			greybox_idx_map.clear();
			criterion_gene_NDM_count = 0;

			for (boxsize_iter = 0; boxsize_iter < boxsize; boxsize_iter++) {
				greybox_idx_map.insert(unordered_map<int, int>::value_type(greybox[0][boxsize_iter], 1));
			}
			
			// get 'ro' value and update result NDM
#pragma omp parallel for reduction(+:criterion_gene_NDM_count) private(cell_iter2, nxyk, ro) schedule(guided)
			for (gene_iter2 = 1; gene_iter2 < process_Column_Size; gene_iter2++) {
				nxyk = 0;
				ro = 0;
				unordered_map<int, int>::iterator it;

				string one_line;
				int one_data;
				string file_path = filepath + "Result_Temp_" + to_string(process_gene_idx + gene_iter2 - 1) + ".csv";
				fstream file;
				vector<int> temp_result_v;
				temp_result_v.clear();

				for (boxsize_iter = 0; boxsize_iter < boxsize; boxsize_iter++) {
					it = greybox_idx_map.find(greybox[gene_iter2][boxsize_iter]);
					if (it != greybox_idx_map.end())
						nxyk++;
				}

				ro = nxyk * ro_temp1 - ro_temp2;
				if (ro > icdf) {
					// update NDM value in result files
					criterion_gene_NDM_count++;

					file.open(file_path, fstream::in|fstream::out);
					
					if (file.fail())
						cout << strerror(errno) << endl;
					else if (file.is_open()) {
						getline(file, one_line);
						stringstream one_line_stream(one_line);
						string one_data_string;

						while (getline(one_line_stream, one_data_string, ',')) {
							stringstream one_data_stream(one_data_string);
							one_data_stream >> one_data;
							temp_result_v.push_back(one_data);
						}
						++temp_result_v[cell_iter1];

						file.clear();
						file.fstream::seekg(0, ios::beg);

						file << temp_result_v[0];
						for (cell_iter2 = 1; cell_iter2 < temp_result_v.size(); cell_iter2++) {
							file << "," << temp_result_v[cell_iter2];
						}
						file.close();

					}

				}

			}

			// update NDM value on criterion gene in result files
			string one_line;
			int one_data;
			string file_path = filepath + "Result_Temp_" + to_string(criterion_gene_idx) + ".csv";
			fstream file;
			vector<int> temp_result_v;
			temp_result_v.clear();
			file.open(file_path, fstream::in | fstream::out);
			if (file.fail())
				cout << "Fail to Read File " << criterion_gene_idx << endl;
			else if (file.is_open()) {
				getline(file, one_line);
				stringstream one_line_stream(one_line);
				string one_data_string;

				while (getline(one_line_stream, one_data_string, ',')) {
					stringstream one_data_stream(one_data_string);
					one_data_stream >> one_data;
					temp_result_v.push_back(one_data);
				}
				temp_result_v[cell_iter1] += criterion_gene_NDM_count;

				file.clear();
				file.fstream::seekg(0, ios::beg);

				file << temp_result_v[0];
				for (cell_iter2 = 1; cell_iter2 < temp_result_v.size(); cell_iter2++) {
					file << "," << temp_result_v[cell_iter2];
				}

				file.close();
			}

			greybox.clear();
		}

		cv->notify_one();
		cv->wait(lock);
	}

}

int main() {
	cout << "=== CSN NDM Implementation ===" << endl;
	cout << "OpenMP Thread: " << omp_get_max_threads() << endl;

	bool file_open_failed = true;
	string input_file_path;
	string file_directory;
	double boxsize_ratio;
	double alpha;

	condition_variable cv;
	mutex m;

	cudaSetDevice(0);

	while (file_open_failed) {
		cout << "File Path : ";
		cin >> input_file_path;
		ifstream file(input_file_path);
		if (file.fail()) cout << "Failed Opening File. Try Again." << endl;
		else file_open_failed = false;
	}
	cout << "Buffer Size (GB) : ";
	cin >> buffer_size_in_GB;
	cout << "boxsize (0. ~ 1) : ";
	cin >> boxsize_ratio;
	cout << "alpha : ";
	cin >> alpha;

	/*
	 *
	 *	[ Initialization ]
	 *
	 */
	file_directory = InitializeFunction(input_file_path, boxsize_ratio, alpha);

	/*
	 *
	 *	[ Start Processing ]
	 *
	 */
	thread thrd_DataRead(DataRead, input_file_path, &cv, &m);
	thread thrd_ProcessNDM(ProcessNDM, file_directory, &cv, &m);
	thrd_ProcessNDM.join();
	thrd_DataRead.join();
	
	/*
	*
	*	[ Concatenate Result Files]
	*
	*/
	ConcatResultFiles(file_directory);

	return 0;
}