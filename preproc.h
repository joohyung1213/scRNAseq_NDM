#ifndef PREPROC_H
#define PREPROC_H

#include <string>

std::string FileOpen();
std::string SplitFileDirectory(std::string filepath);
std::string InitializeFunction(std::string filepath, double boxsize_ratio, double sig_level);
void DataRead(std::string filepath);
void RecordResult(std::string file_directory);

#endif