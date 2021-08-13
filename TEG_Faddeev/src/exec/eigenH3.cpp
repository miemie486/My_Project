/* Find the binding energy of a triton. */

#include <complex>
#include <fstream>
#include <iostream>
#include <ctime>
// #include "timer.h"
#include "math.h"
#include <string.h>
#include <vector>

// #include "def.h"
#include "kMatGen.h"
#include "channels.h"
#include "fstr.h"
#include "msg.h"
// #include "mkl.h"

using namespace std;

string fname, default_fname {"eigenH3.in"};
bool chntable_requested = false, flag_nlo = false, flag_n2lo = false;
REAL x1 {0.0}, x2 {0.0};
int local_verbose_level {_VBSLVL_NORMAL_};
string whole_argline;
int parse_cmdline(int argc, char* argv[]);

/* Usage:
  eigenH3 [filename] [options]
(filename does not have to appear before options)

filename: the name of the input file. If filename is omitted,
eigenH3.in will be used.

options:
-c : print out a table for quantum numbers allowed by parameters in filename or eigenH3.in
or

-v vbl : verbose level, 10 >= vbl >= 0

-1 x1 : add x1*V_NLO to V_LO

-2 x2 : add x2*V_N2LO to V_LO

Remarks
* If there is x2 but x1 is absent, x1 will be set 0

*/
int main(int argc, char* argv[]) {

  REAL default_ratio_NQovrNP {1.0};

  // Read input from
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
  inputTable.BE01={-15.0,-8.0};
  inputTable.error={1.0E-2};
  inputTable.ratioNQP = {default_ratio_NQovrNP};

  if (parse_cmdline(argc, argv) == -1) {
    cout << "eigenH3: Something wrong with command line arguments." << endl;
    return -1;
  }
  set_verbose_level_fdv(local_verbose_level);

  // debug
  // cout << fname << " " << x1 << " " << x2 << " "
    // << flag_nlo << " " << flag_n2lo << endl;
  // return 0;

  //This is to judge whether the file reading successful or not
  bool condition=inputTable.cookInput(fname);

  // Generalize channel table
  Chs chsTable {};

  // Parameters needed for generating channels table
  chsTable.writeMaxL(inputTable.maxL[0]);
  chsTable.writeMaxOpT(inputTable.maxOpT[0]);
  chsTable.writeOpJ(inputTable.opJ[0]);
  chsTable.writeParity(inputTable.parity[0]);
  chsTable.genChs();

  if (chntable_requested) {
    chsTable.showChannels(); // Show channel table for user to select
    return 0;
  }

  // Initialize the selected channels
  size_t nChannels = inputTable.channels.size();
  int labelChs[nChannels];
  for (size_t i = 0; i < nChannels; i++) {
    labelChs[i] = inputTable.channels[i];
  }

  REAL **channels = alloc_3Nchn_table(nChannels);
  chsTable.outChannels(channels, nChannels, labelChs);

  // Print out the channels the user have seleted
  cout << endl << endl << "/////////////////////////////////" << endl;
  if(inputTable.channels_j)
  {
    printf("You have selected channels:\n");
  }
  else
  {
    printf("You have selected channels: (default)\n");
  }

  printf("%-8s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s\n",
    "alpha", "l", "s", "j", "lmbd", "I", "T",  "lp", "opT", "Pi", "opJ");
  for(size_t i = 0; i < nChannels; i++) {
    printf(" %4d   ", labelChs[i]);
    for(int j = 0; j < _N_QMN_CHN_; j++)
      printf("%-5.1f", channels[i][j]);
    printf("\n");
  }

  // Initialize k matrix, note that all configurations can be changed later.
  REAL mambda = inputTable.mambda[0]; // Cutoff
  REAL _CTanMesh = mambda * inputTable.CMRatio[0];
  // Define the mesh points' range
  REAL dPara = inputTable.NP[0],
    para0 = inputTable.NP[1],
    para1 = inputTable.NP[2];
  // size_t lenPara = (para1 - para0) / dPara + 1;

  cout << "/////////////////////////////////" << endl;

  show_message_fdv("* eigenH3: cmdline = " + whole_argline, _VBSLVL_LOW_);
  //print parameters this cacualating will be used
  if(condition)
  {
    cout << "Input file = " << fname << endl;
  }
  else
  {
    cout << "Input file = Nothing" << endl;
    cout << "Those parameters are default parameters." << endl;
  }
  inputTable.PrintParameter();

  int uptoQn {0}, numPara {1};
  REAL potpara[NUM_PARA_TWOBODY] {0.0};

  if (flag_n2lo) {
    numPara = 3;
    uptoQn = 2;
    if (flag_nlo) {
      potpara[1] = x1;
    } else {
      potpara[1] = 0.0;
      flag_nlo = true;
    }
    potpara[2] = x2;
  } else {
    if (flag_nlo) {
      numPara = 2;
      uptoQn = 1;
      potpara[1] = x1;
    }
  }
  potpara[0] = uptoQn;

  if (flag_nlo) {
    printf("* x1 = %-8.4e\n", x1);
  }
  if (flag_n2lo) {
    printf("* x2 = %-8.4e\n", x2);
  }

  cout << "/////////////////////////////////" << endl;
    // current date/time based on current system
  time_t now = time(0);
  cout << "==>>>> Job starts at " << ctime(&now) << endl;
  printf("%-20s%-20s\n", "Mesh Points", "Binding Energy (MeV)");
  fflush(stdout);

  REAL BEH3, local_eps {1.0E-6};
  int NP = 15; // number of p mesh pts

  // Initialize chengdu_nnscat global functions/variables
  init_nnscat_chengdu_C();

  // Initialize TEG_Faddeev global functions/variables
  facLabel();

  // Complex-valued pots contains "cmplx"
  if ( inputTable.potential[0].find("cmplx") != string::npos ) {

    kMatrix<COMPLEX> KMatC {nChannels, _CTanMesh, channels, inputTable.opJ[0]};
    REAL phi = inputTable.phi[0]/180.0*3.141593;
    const std::complex<double> ii(0.0, 1.0);
    KMatC.setFacContour(exp(-ii*phi));
    KMatC.setPotName(inputTable.potential[0]); // potential name
    KMatC.loadTwoBodyPara(numPara, potpara);
    KMatC.setMass(_AvgMassN_);
    KMatC.setNX(inputTable.NX[0]);

    for (NP = para0; NP < para1 + local_eps; NP += dPara){
      KMatC.setNP(NP);
      KMatC.setNQ(int(NP*inputTable.ratioNQP[0]));
      KMatC.genGTable();
      KMatC.allocKMat();
      KMatC.initKMatUnit();
      KMatC.genIntzTable_Hmgns();

      // debug
      REAL ImE0 {0.0}, ImE1 {0.01};
      COMPLEX detK, pole;

      printf("NP = %-5d\n", NP);
      pole = KMatC.findResPole(inputTable.BE01[0], inputTable.BE01[1], ImE0, ImE1,
        inputTable.error[0], inputTable.error[0]*0.1, detK);

      printf("* [Final] NP = %-5d", NP);
      printf("H3BE = (%-22.15e,  %-22.15e)\n", pole.real(), pole.imag());
      fflush(stdout);

      // debug
      // cout <<  "det(E) = " << KMatC.findKHmgnsDetFunc(BEH3) << endl;
      cout <<  "detK(E) = " << detK << endl;
      // cout <<  "Check: Epole = " << pole << ",  det(E) = "
      //   << KMatC.findKHmgnsDetFunc(pole) << endl;

      now = time(0);
      cout << endl << "------ Currently  at " << ctime(&now)  << endl;
      fflush(stdout);
    }
  } else {

    kMatrix<REAL> inst_kMatR {nChannels, _CTanMesh, channels, inputTable.opJ[0]};
    inst_kMatR.setFacContour(1.0);
    inst_kMatR.setPotName(inputTable.potential[0]); // potential
    inst_kMatR.loadTwoBodyPara(numPara, potpara);
    inst_kMatR.setMass(_AvgMassN_);
    inst_kMatR.setNX(inputTable.NX[0]);

    // Find the binding energy
    for(NP = para0; NP < para1 + local_eps; NP += dPara){
      inst_kMatR.setNP(NP);
      inst_kMatR.setNQ(int(NP*inputTable.ratioNQP[0]));
      inst_kMatR.genGTable();
      inst_kMatR.allocKMat();
      inst_kMatR.initKMatUnit();
      inst_kMatR.genIntzTable_Hmgns();

      printf("NP = %-5d\n", NP);
      BEH3 = inst_kMatR.findBE(inputTable.BE01[0], inputTable.BE01[1],
        inputTable.error[0], inputTable.RFMTHD[0]);

      printf("* [Final] NP = %-5d", NP);
      printf("H3BE = %-22.15e\n", BEH3);

      fflush(stdout);
      now = time(0);
      cout << endl << "------ Currently  at " << ctime(&now) << endl;
      fflush(stdout);
    }
  }

  now = time(0);
  cout << "<<<<<< Job ends at   " << ctime(&now)  << endl;

  // for(int i = 0; i < nChannels; i++){
  //   delete[] channels[i];
  // }
  // delete[] channels;
  dealloc_3Nchn_table(nChannels, channels);
  return 0;

}

// Return:
//   -1 -> there is something wrong in command line arguments
//   0  -> success
int parse_cmdline(int argc, char* argv[]) {

  if (argc < 2) {
    fname = default_fname; // eigenH3
    return 0;
  }

  vector<string> arg_lst, leftover;
  whole_argline = argv[0];
  for (unsigned i = 1; i < (unsigned) argc; i++) {
    arg_lst.push_back(argv[i]);
    whole_argline += " ";
    whole_argline += argv[i];
  }
  leftover = arg_lst;

  // for (vector<string>::iterator it = arg_lst.begin(); it != arg_lst.end(); ++it) {
  vector<string>::iterator it = leftover.begin();
  while (it != leftover.end()) {

    if (*it == "-c") {
      chntable_requested = true;
      leftover.erase(it);
      continue;
    }

    if (*it == "-v") {
      leftover.erase(it);
      if (it != leftover.end()) {
        local_verbose_level = (int) stod(*it);
        leftover.erase(it);
      } else {
        cout << "eigenH3: couldn\'t get verbose value. use default value" << endl;
        local_verbose_level = _VBSLVL_NORMAL_;
      }
      continue;
    }

    if (*it == "-1") {
      leftover.erase(it);
      if (it != leftover.end()) {
        x1 = stod(*it);
        flag_nlo = true;
        leftover.erase(it);
      } else {
        cout << "eigenH3: couldn\'t get NLO x-value" << endl;
        return -1;
      }
      continue;
    }

    if (*it == "-2") {
      leftover.erase(it);
      if (it != leftover.end()) {
        x2 = stod(*it);
        flag_n2lo = true;
        leftover.erase(it);
      } else {
        cout << "eigenH3: couldn\'t get N2LO x-value" << endl;
        return -1;
      }
      continue;
    }

    ++it;
  }

  // if nothing left, use default filename; what is left must be the filename
  if (leftover.empty()) {
    fname = default_fname;
    return 0;
  }

  // Use the first non-option string in the command line as filename
  fname = leftover.front();
  return 0;
}


