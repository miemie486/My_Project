/* Reade a csv-formatted file, and display the channels enlisted. */

#include <fstream>
#include <iostream>
#include "channels.h"
#include "fstr.h"
#include <string.h>

using namespace std;

int main(int argc, char* argv[]){
  // Read input from input.csv
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

  if ( (argc == 2) && (strcmp(argv[1], "-h") == 0) ) {
    cout << "Usage: listalpha filename\n";
    return 0;
  }
  if ( (argc == 1) || (argv[1][0] == '\0') )
  {
    if (! inputTable.cookInput("alphas.in")) return 0;
  }
  else
  {
    if (! inputTable.cookInput(argv[1])) return 0;
  }

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
  // Print out the channels the user have seleted
  chsTable.outChannels(channels, nChannels, labelChs);
  printf("You have selected channels:\n");
  printf("%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\n", "alpha", "l", "s", "j", "lambda", "I", "T",  "Coupling", "opT");
  for(int i = 0; i < nChannels; i++){
    printf("%-10d", labelChs[i]);
    for(int j = 0; j < _N_QMN_CHN_; j++)
      printf("%-10.2f", channels[i][j]);
    printf("\n");
  }

  for(int i = 0; i < nChannels; i++){
    delete[] channels[i];
  }
  delete[] channels;
}
