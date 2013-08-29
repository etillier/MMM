// classDistanceMatrixFile.h
// Version 2010.01.29
// (c) 2010, Author: Robert L. Charlebois

#ifndef CLASS_DISTANCE_MATRIX_FILE_H
#define CLASS_DISTANCE_MATRIX_FILE_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

template<typename T> class DistanceMatrixFile {
protected:
	std::string fileName_;
	std::vector<std::string> names_;
	std::vector<std::vector<T> > distances_;
public:
	explicit DistanceMatrixFile(const std::string& matrixFileName) :
		fileName_(matrixFileName)
	{
		if (!matrixFileName.empty()) {
			std::ifstream matrixFile(matrixFileName.c_str());
			if (matrixFile.is_open()) {
				try {
					std::string label;
					T dist;
					std::size_t numEntries;
					matrixFile >> numEntries;
					distances_.resize(numEntries);
					names_.reserve(numEntries);
					for (std::size_t j, i = 0; matrixFile && i < numEntries && matrixFile; ++i) {
						matrixFile >> label;
						names_.push_back(label);
						std::vector<T>& row = distances_[i];
						row.reserve(numEntries);
						for (j = 0; matrixFile && j < numEntries; ++j) {
							matrixFile >> dist;
							row.push_back(dist);
						}
						if (row.size() != numEntries) throw 1;
					}
					if (names_.size() != numEntries) throw 2;
				} catch (...) {
					names_.clear();
					distances_.clear();
					std::cout << "Error reading " << matrixFileName << std::endl;
				}
			} else {
				std::cout << "Error: Cannot open " << matrixFileName << std::endl;
			}
		} else {
			std::cout << "Error: No matrix file name specified" << std::endl;
		}
	}
	
	virtual ~DistanceMatrixFile() { }
	
	bool hasData() const { return !distances_.empty(); }
	std::size_t size() const { return distances_.size(); }
	const std::string& fileName() const { return fileName_; }
	const std::string& name(int i) const { return names_.at(i); }
	T distance(int i, int j) const { return distances_.at(i).at(j); }
};

#endif
