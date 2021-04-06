#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include "utils.h"

using namespace std;

extern int gene_num, cell_num, boxsize;
extern float icdf;
extern double ro_temp1, ro_temp2;
extern vector <vector <float>> data_vector;
vector <float> data_vector_line;
extern vector <vector<int>> NDM;

string FileOpen() {
	string input_file_path;
	bool file_open_failed = true;
	while (file_open_failed) {
		cout << "FILE Path : ";
		cin >> input_file_path;
		ifstream file(input_file_path);
		if (file.fail()) cout << "Failed Opening File. Try Again." << endl;
		else file_open_failed = false;
	}
	return input_file_path;
}

string SplitFileDirectory(string filepath) {
	vector<string> tokens;
	stringstream stream(filepath);
	string path_chunk;
	string directory;

	// Directory form with / slash
	if (filepath.find('/') != string::npos) {
		while (getline(stream, path_chunk, '/'))
			tokens.push_back(path_chunk);
		tokens.pop_back();

		for (auto token : tokens)
			directory = directory + token + '/';
	}
	// Directory form with \ backslash
	else if (filepath.find("\\") != string::npos) {
		while (getline(stream, path_chunk, '\\'))
			tokens.push_back(path_chunk);
		tokens.pop_back();

		for (auto token : tokens)
			directory = directory + token + '\\';
	}

	// directory without filename
	return directory;
}

string InitializeFunction(string filepath, double boxsize_ratio, double sig_level) {
	int i = 0;
	string file_directory;
	ifstream file(filepath);

	if (file.fail()) cout << "Fail to Read File" << endl;
	else if (file.is_open()) {
		cout << "File Readed" << endl;
		cout << "Initializing..." << endl;

		/*
		*
		*	[ Scan Data File Once to Get (Row, Column) Size ]
		*	[ Initialize Global Variables ]
		*
		*/

		// Count Column Number
		string one_line;
		getline(file, one_line);
		stringstream line_stream(one_line);
		string one_data;
		while (getline(line_stream, one_data, ',')) {
			stringstream convertor(one_data);
			cell_num++;
		}
		gene_num++;

		// Count Row Number
		while (getline(file, one_line)) {
			gene_num++;
		}

		// Initialize Global Variables
		boxsize = round(boxsize_ratio * cell_num) + 1;
		icdf = -NormalCDFInverse(sig_level);

		ro_temp1 = cell_num * sqrt((double)(cell_num - 1)) / ((boxsize - 1) * (cell_num - (boxsize - 1)));
		ro_temp2 = sqrt((double)(cell_num - 1)) * (boxsize - 1) / (cell_num - (boxsize - 1));

		vector <int> NDM_line(cell_num, 0);
		for (i = 0; i < gene_num; i++) NDM.push_back(NDM_line);

	}

	file_directory = SplitFileDirectory(filepath);

	cout << "Number of Gene : " << gene_num << endl;
	cout << "Number of Cell : " << cell_num << endl;
	cout << "GreyBox Size: " << boxsize << endl;
	cout << "icdf (Significant Level) : " << icdf << endl;
	// cout << "Reading Buffer Size : " << buffer_size_gene_num << " rows(genes)\n" << endl;

	return file_directory;

}

void DataRead(string filepath) {
	cout << "[DataRead] Start" << endl;

	ifstream file(filepath);
	string one_line;
	double one_data;

	if (file.fail()) cout << "[DataRead] Fail to Read File" << endl;
	else if (file.is_open()) {
		cout << "[DataRead] File Readed\n";

		// Start Reading Data
		data_vector_line.clear();
		data_vector.clear();

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
		}
		cout << "[DataRead] Read Done" << endl;

	}

}

void RecordResult(string file_directory) {
	ofstream result_file;
	string result_file_path = file_directory + "NDM_Result.csv";
	result_file.open(result_file_path);

	for (int gene_iter = 0; gene_iter < gene_num; gene_iter++) {
		for (int cell_iter = 0; cell_iter < cell_num; cell_iter++) {
			result_file << NDM[gene_iter][cell_iter] << ",";
		}
		result_file << "\n";
	}
}