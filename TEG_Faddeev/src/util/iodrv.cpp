// Bingwei Long  09/03/2020

#include "iodrv.h"

// using namespace std;

void iodrv::parse_line_to_reals(string &str_nums, vector<REAL> &num_array, char delim) {

  istringstream ss(str_nums);
  string token;
  while(getline(ss, token, delim)) {
    num_array.push_back(stod(token));
  }
}
