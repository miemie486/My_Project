#ifndef CSVREADER
#define CSVREADER

#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

/*
 * A class to read data from a csv file.
 */
class CSVReader
{
	std::string fileName;
	std::string delimeter;

 public:
 CSVReader(std::string filename, std::string delm = ",") :
  fileName(filename), delimeter(delm)
    { }

	// Function to fetch data from a CSV File
	// std::vector<std::vector<std::string> > getData();
  bool getData(std::vector<std::vector<std::string>> &);
};

#endif
