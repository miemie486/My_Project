/* Find corrections to the binding energy of the triton. */

// #include <complex>
#include "matrixK.h"
#include "mkl.h"
#include "kMatGen.h"
#include <fstream>
#include <iostream>
#include "timer.h"
#include "channels.h"
#include "fstr.h"
#include <string.h>

using namespace std;

/* Usage: sublH3 x filename
where x is normally a small number < 1.0, filename is the name of the input
file. If filename not provided, sublH3_input.csv will be used */
int main(int argc, char* argv[]){

  // Read input from filename or sublH3_input.csv
  Fstr inputTable; // Input parameters
  //Set default parameters
  inputTable.maxL={2};
  inputTable.maxOpT={0.5};
  inputTable.opJ={0.5};
  inputTable.parity={0};
  inputTable.channels={0,1,12,2,13,5};
  inputTable.mambda={400};
  inputTable.NP={3,8,11};
  inputTable.potential={"chengdu_MMWLY_400"};
  inputTable.NX={7};
  inputTable.CMRatio={0.6};
  inputTable.RFMTHD={"Scant"};
  inputTable.BE01={-15.0,-3.0};
  inputTable.error={1.0E-2};

  string fname;
  if ( (argc < 3) || (argv[2][0] == '\0') )
  {
    fname = "sublH3_input.csv";
  }
  else
  {
    fname = argv[2];
  }
  if (! inputTable.cookInput(fname)) return 0;

  // Generalize channel table
  Chs chsTable {};

  chsTable.writeMaxL(inputTable.maxL[0]); // Several parameters are needed for generating channels table
  chsTable.writeMaxOpT(inputTable.maxOpT[0]);
  chsTable.writeOpJ(inputTable.opJ[0]);
  chsTable.writeParity(inputTable.parity[0]);
  chsTable.genChs();
  chsTable.showChannels(); // Show channel table for user to select

  // Initialize the selected channels
  int nChannels = inputTable.channels.size();
  int labelChs[nChannels];
  for(int i = 0; i < nChannels; i++)
    labelChs[i] = inputTable.channels[i];
  REAL **channels = new REAL*[nChannels];
  for(int i = 0; i < nChannels; i++){
    channels[i] = new REAL[_N_QMN_CHN_];
  }
  chsTable.outChannels(channels, nChannels, labelChs);

  // Print out the channels the user have seleted
  printf("You have selected channels:\n");
  printf("%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n", "alpha", "l", "s", "j", "lambda", "I", "T",  "Coupling", "opT");
  for(int i = 0; i < nChannels; i++){
    printf("%-10d", labelChs[i]);
    for(int j = 0; j < _N_QMN_CHN_; j++)
      printf("%-10.2f", channels[i][j]);
    printf("\n");
  }

  double x, potpara[10] {0.0};
  int uptoQn = 1, numPara = 2;
  if (argc >= 2) {
    x = strtod(argv[1], NULL);
  } else {
    x = 0.0;
  }
  potpara[0] = uptoQn;
  potpara[1] = x;

  // Initialize k matrix, note that all configurations can be changed later.
  int NP = 15; // numer of p mesh
  REAL mambda = inputTable.mambda[0]; // Cutoff
  REAL _CTanMesh = mambda * inputTable.CMRatio[0];

  cout << "/////////////////////////////////" << '\n';
  cout << "Input file = " << fname << '\n';
  cout << "Potential = " << inputTable.potential[0] << '\n';
  cout << "Mambda = " << mambda << '\n';
  cout << "CTanMesh/Mambda = " << inputTable.CMRatio[0] << '\n';
  cout << "x at cmd line = " << x << '\n';
  cout << "/////////////////////////////////" << '\n';
  fflush(stdout);

  kMatrix<REAL> kMatR {NP, nChannels, _CTanMesh, channels, inputTable.opJ[0]}; // The final one is opJ
  kMatR.setFacContour(1.0);
  kMatR.setPotName(inputTable.potential[0]); // potential
  kMatR.loadTwoBodyPara(numPara, potpara);
  kMatR.setCTanMesh(_CTanMesh);
  // Define the mesh points' range
  REAL dPara = inputTable.NP[0], para0 = inputTable.NP[1], para1 = inputTable.NP[2];
  size_t lenPara = (para1 - para0) / dPara + 1;
  REAL yArr[lenPara], xArr[lenPara]; // record the number of mesh points in xArr, and the correspond binding energy in yArr
  int idx = 0;

  // Timer clock (lenPara + 1); // Test the speed

  // clock.record("");
  // Find the binding energy
  printf("%-20s%-20s\n", "Mesh Points", "Binding Energy (MeV)");
  for(NP = para0; NP < para1; NP += dPara){
    kMatR.setNP(NP);
    kMatR.setNX(6);
    kMatR.smartUpdate();
    yArr[idx] = kMatR.findBE(-5, -10, 0.01, "Secant"); // You can change the method "Secant" to "Newton", which is slower but more stable
    xArr[idx] = NP;
    printf("%-20.1f", xArr[idx]);
    printf("%-20.4e\n", yArr[idx]);
    // clock.print();
    // clock.record("");
    fflush(stdout);
    idx++;
  }
  // clock.record("");
  // clock.print();

  for(int i = 0; i < nChannels; i++){
    delete[] channels[i];
  }
  delete[] channels;
}

