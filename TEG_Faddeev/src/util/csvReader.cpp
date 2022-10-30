#include "csvReader.h"
#include <iostream>
/*
 * Parses through csv file line by line and returns the data
 * in vector of vector of strings.
 */
bool CSVReader::getData(std::vector<std::vector<std::string>>& dataList)
{
	std::ifstream file(fileName);
  // std::vector<std::vector<std::string> > dataList;
  std::string line = "";
  bool successful;

	if (file.fail()) 
  {
    //std::cout << "Can\'t open " << fileName << '\n';
    successful = false;
  } 
  else 
  {

    // Iterate through each line and split the content using delimeter
    while (getline(file, line))
      {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        std::stringstream lineStream (line);
        std::string word = "";
        std::vector<std::string> vec;
        while(std::getline(lineStream, word, ','))
          vec.push_back(word);
        dataList.push_back(vec);
      }
    successful = true;
  }
	// Close the File
	file.close();
  return successful;
}
