#ifndef IODRV
#define IODRV

#include <sstream>
#include <fstream>
#include <string.h>
#include <vector>

#include "def.h"

using namespace std;

namespace iodrv {
  void parse_line_to_reals(string &str_nums, vector<REAL> &num_array, char delim);
  void discard_comments(char cmmnt_char, ifstream &input_file);
}

#endif
