/* Read inputs in input.csv, which includes parameters required by solving the Faddeev Equation */
#include "csvReader.h"

class Fstr
{
  public:
  std::vector<double> maxL;
  std::vector<double> maxOpT;
  std::vector<double> opJ;
  std::vector<double> parity;
  std::vector<double> channels;
  std::vector<double> mambda;
  std::vector<double> NP;
  std::vector<std::string> potential;
  std::vector<double> NX;
  std::vector<double> CMRatio;
  std::vector<std::string> RFMTHD;
  std::vector<double> BE01;
  std::vector<double> error;
  std::vector<double> phi;
  std::vector<double> ratioNQP; // NQ/NP in Faddeev eqn

  //Judge that whether those data reading successful or not

  bool maxL_j=false;
  bool maxOpT_j=false;
  bool opJ_j=false;
  bool parity_j=false;
  bool channels_j=false;
  bool mambda_j=false;
  bool NP_j=false;
  bool potential_j=false;
  bool NX_j=false;
  bool CMRatio_j=false;
  bool RFMTHD_j=false;
  bool BE01_j=false;
  bool error_j=false;
  bool phi_j=false;
  bool ratioNQP_j = false;

  // Read input file can collect required data
  bool cookInput(std::string path);
  void PrintParameter();
};
