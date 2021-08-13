#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <string.h>
#include <vector>
#include <stdio.h>

#include "def.h"
#include "kMatGen.h"
#include "channels.h"
#include "iodrv.h"
#include "msg.h"
#include "delta.h"

// Usage:
// elasuidelta.exe elasuidelta.in -q 28.89 -NP 20 -v 10 -1 0.01

int parse_cmdline(int argc, char* argv[]);
int read_inputs(string &fname);
bool inputparas_are_ready(void);
void print_out_paras(void);

string fname, default_fname {"elasuidelta.in"};
REAL maxL, maxOpT, opJ, parity, mambda, cmratio, phi;
int NP, NX;

vector<double> alpha_lst_bare;
size_t nChannels;
REAL** channels;
bool chntable_requested {false};

string pot_name;
int uptoQn {0}, numPara {1}; // By default, LO and no paras
REAL potpara[NUM_PARA_TWOBODY] {0.0};
bool flag_nlo {false};
REAL x1 {0.0};
// bool flag_n2lo {false};
// REAL x2 {0.0};

int local_verbose_level {_VBSLVL_NORMAL_};
string whole_argline;
REAL q0;
bool flag_q0 {false}, flag_NP {false};

// This program is built to calculate S matrix, phase shifts and mixangles in nd scattering.
// Necessary parameters needed to run this program contains NP, q, and
// an input file which lists channels and gives parity, opJ, potential name, mambdad, CMRatio, phi, NX.

int main(int argc, char* argv[]) {

  if (parse_cmdline(argc, argv) == -1) {
    cout << "elasuidelta: Something wrong with command line arguments." << endl;
    return -1;
  }
  set_verbose_level_fdv(local_verbose_level);
  if (! inputparas_are_ready()) {
    return -1;
  }

  read_inputs(fname);
  if (local_verbose_level == _VBSLVL_HIGH_) {
    print_out_paras();
    fflush(stdout);
  }

  // prepare 3N channels
  Chs chsTable(maxL, maxOpT, opJ, parity);
  chsTable.genChs();
  if (chntable_requested) {
    chsTable.showChannels(); // Show channel table for user to select
    return 0;
  }
  vector<int> alpha_lst(alpha_lst_bare.begin(), alpha_lst_bare.end());
  nChannels = alpha_lst.size();
  channels = alloc_3Nchn_table(nChannels);
  cout << "/////////////////////////////////" << endl;
  cout << "You have selected channels:" << endl;
  chsTable.showChannels(alpha_lst);
  cout << "/////////////////////////////////" << endl;
  chsTable.outChannels(channels, nChannels, alpha_lst);

  // Initialize chengdu_nnscat global functions/variables
  init_nnscat_chengdu_C();
  // Initialize TEG_Faddeev global functions/variables
  facLabel();

  REAL _CTanMesh = mambda * cmratio;
  kMatrix<COMPLEX> kMatC(nChannels, _CTanMesh, channels, opJ);
  const std::complex<double> ii(0.0, 1.0);
  kMatC.setFacContour(exp(-ii*phi/180.0*M_PI));
  kMatC.setPotName(pot_name);
  kMatC.loadTwoBodyPara(numPara, potpara);
  kMatC.setMass(_AvgMassN_);
  kMatC.setNX(NX);
  kMatC.setNP(NP);
  kMatC.setNAuxMesh((int) NP*1.4);
  kMatC.setNQ((int) NP*1.0);

  kMatC.genGTable();
  kMatC.allocKMat();
  kMatC.initKMatUnit();
  get_Upara_list(opJ, parity);
  int Size_S = (opJ == 0.5)? 4:9;
  COMPLEX U_array[Size_S];
  kMatC.get_U_array(parity, q0, U_array);
  COMPLEX Smat[Size_S];
  SMatGen(Smat, Size_S, opJ, parity, q0, U_array);
  print_out_paras();
  cout.precision(10);
  cout<<"Matrix of S: "<<endl;
  for(int i = 0; i < Size_S; i++)
    {cout<<"        S["<<i<<"] = "<<Smat[i]<<endl;}

  int Size_delta = (opJ == 0.5)? 2:3;
  int Size_mixangle = (opJ == 0.5)? 1:3;
  COMPLEX delta[Size_delta], mixangle[Size_mixangle];
  Gen_deltaArray(delta, Size_delta, mixangle, Size_mixangle, Smat, Size_S);
  cout<<"phase shifts :"<<endl;
  for(int i = 0; i < Size_delta; i++)
    {cout<<"        delta["<<i<<"] = "<<delta[i]<<endl;}
  cout<<"mixing angles :"<<endl;
  for(int i = 0; i < Size_mixangle; i++)
    {cout<<"        mixangle ["<<i<<"] = "<<mixangle[i]<<endl;}
  dealloc_3Nchn_table(nChannels, channels);
  return 0;
}


// Return:
//   -1 -> reading fname went wrong
//   0  -> success
int read_inputs(string &fname) {

  ifstream input_file;
  input_file.open(fname.c_str());
  if (! input_file) {
    cout << "read_inputs: couldn\'t read " << fname << endl;
    return -1;
  }

  string comment_line, tmpbuff;
  getline(input_file, comment_line);
  input_file >> maxL >> maxOpT >> opJ >> parity;
  getline(input_file, tmpbuff);

  getline(input_file, comment_line);
  getline(input_file, tmpbuff);
  iodrv::parse_line_to_reals(tmpbuff, alpha_lst_bare, ',');

  getline(input_file, comment_line);
  getline(input_file, pot_name);

  getline(input_file, comment_line);
  input_file >> mambda >> cmratio >> phi;
  getline(input_file, tmpbuff);

  getline(input_file, comment_line);
  input_file >> NX;

  input_file.close();
  return 0;
}

// Return:
//   -1 -> there is something wrong in command line arguments
//   0  -> success
int parse_cmdline(int argc, char* argv[]) {

  if (argc < 2) {
    fname = default_fname;
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
        cout << "elasuidelta: couldn\'t get verbose value. use default value" << endl;
        local_verbose_level = _VBSLVL_NORMAL_;
      }
      continue;
    }

    if (*it == "-q") {
      leftover.erase(it);
      if (it != leftover.end()) {
        q0 = stod(*it);
        flag_q0 = true;
        leftover.erase(it);
      } else {
        cout << "elasuidelta: couldn\'t get value for q0" << endl;
        return -1;
      }
      continue;
    }

    if (*it == "-NP") {
      leftover.erase(it);
      if (it != leftover.end()) {
        NP = (int) stod(*it);
        flag_NP = true;
        leftover.erase(it);
      } else {
        cout << "elasuidelta: couldn\'t get value for NP" << endl;
        return -1;
      }
      continue;
    }

    if (*it == "-1") {
      leftover.erase(it);
      if (it != leftover.end()) {
        x1 = stod(*it);
        flag_nlo = true;
        numPara = 2;
        uptoQn = 1;
        potpara[1] = x1;
        potpara[0] = uptoQn;
        leftover.erase(it);
      } else {
        cout << "elasuidelta: couldn\'t get NLO x-value" << endl;
        return -1;
      }
      continue;
    }

    // if (*it == "-2") {
    //   leftover.erase(it);
    //   if (it != leftover.end()) {
    //     x2 = stod(*it);
    //     flag_n2lo = true;
    //     leftover.erase(it);
    //   } else {
    //     cout << "elasuidelta: couldn\'t get N2LO x-value" << endl;
    //     return -1;
    //   }
    //   continue;
    // }

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

bool inputparas_are_ready() {

  if (! flag_NP) {
    cout << "elasuidelta: NP is not ready" << endl;
    return false;
  }
  if (! flag_q0) {
    cout << "elasuidelta: q0 is not ready" << endl;
    return false;
  }

  return true;
}

void print_out_paras(void) {

  show_message_fdv("* elasuidelta: cmdline = " + whole_argline, _VBSLVL_LOW_);
  cout << "* NP = " << NP << endl;
  // debug
  // cout << "* NX = " << NX << endl;
  // cout << "* mambda = " << mambda << endl;
  // cout << maxL << endl;
  // cout << maxOpT << endl;
  // cout << alpha_lst_bare[0] << ", " << alpha_lst_bare[5] << endl;
  // cout << pot_name << endl;
  cout << "* OPJ = " << opJ << endl;
  cout << "* parity = " << parity << endl;
  cout << "* q0 = " << q0 << endl;
  if (flag_nlo) {
    printf("* x1 = %-8.4e\n", x1);
  }
  // if (flag_n2lo) {
  //   printf("* x2 = %-8.4e\n", x2);
  // }
}
