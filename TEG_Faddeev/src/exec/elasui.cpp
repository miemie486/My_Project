// Bingwei Long  09/03/2020

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

// Usage:
// elasui.exe elasui.in -NP 20 -v 10 -q 1.0 -l 0 -I 0.5

int parse_cmdline(int argc, char* argv[]);
int read_inputs(string &fname);
bool inputparas_are_ready(void);
void print_out_paras(void);

string fname, default_fname {"elas_ui.in"};
REAL maxL, maxOpT, opJ, parity, mambda, cmratio, phi;
int NP, NX;

vector<double> alpha_lst_bare;
size_t nChannels;
REAL** channels;
bool chntable_requested {false};

string pot_name;
int local_verbose_level {_VBSLVL_NORMAL_};

string whole_argline;
REAL q0, I_p, lambda_p, I, lambda;
bool flag_q0 {false}, flag_Ip {false}, flag_I {false}, flag_lambda {false},
  flag_lambdap {false}, flag_NP {false};

int main(int argc, char* argv[]) {

  if (parse_cmdline(argc, argv) == -1) {
    cout << "elasui: Something wrong with command line arguments." << endl;
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

  // debug
  // return 0;

  // Initialize chengdu_nnscat global functions/variables
  init_nnscat_chengdu_C();
  // Initialize TEG_Faddeev global functions/variables
  facLabel();

  REAL _CTanMesh = mambda * cmratio;
  kMatrix<COMPLEX> kMatC(nChannels, _CTanMesh, channels, opJ);
  const std::complex<double> ii(0.0, 1.0);
  kMatC.setFacContour(exp(-ii*phi/180.0*M_PI));
  kMatC.setPotName(pot_name);
  kMatC.setMass(_AvgMassN_);
  kMatC.setNX(NX);
  kMatC.setNP(NP);
  kMatC.setNAuxMesh((int) NP*1.4);
  kMatC.setNQ((int) NP*1.0);

  kMatC.genGTable();
  kMatC.allocKMat();
  COMPLEX Ukf = kMatC.get_U_elastic(q0, lambda_p, I_p, lambda, I);

  print_out_paras();
  printf("* U_lambda_I = (%-22.15e,  %-22.15e)\n", Ukf.real(), Ukf.imag());
  Ukf = Ukf*2.0/3.0*M_PI*_AvgMassN_*_HBARC_;
  printf("* 2pi/3*mN*U = (%-22.15e,  %-22.15e)\n", Ukf.real(), Ukf.imag());

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
        cout << "elas_ui: couldn\'t get verbose value. use default value" << endl;
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
        cout << "elas_ui: couldn\'t get value for q0" << endl;
        return -1;
      }
      continue;
    }

    if (*it == "-I") {
      leftover.erase(it);
      if (it != leftover.end()) {
        I = stod(*it);
        flag_I = true;
        leftover.erase(it);
      } else {
        cout << "elas_ui: couldn\'t get value for I" << endl;
        return -1;
      }
      continue;
    }

    if (*it == "-Ip") {
      leftover.erase(it);
      if (it != leftover.end()) {
        I_p = stod(*it);
        flag_Ip = true;
        leftover.erase(it);
      } else {
        cout << "elas_ui: couldn\'t get value for I_p" << endl;
        return -1;
      }
      continue;
    }

    if (*it == "-l") {
      leftover.erase(it);
      if (it != leftover.end()) {
        lambda = stod(*it);
        flag_lambda = true;
        leftover.erase(it);
      } else {
        cout << "elas_ui: couldn\'t get value for lambda" << endl;
        return -1;
      }
      continue;
    }

    if (*it == "-lp") {
      leftover.erase(it);
      if (it != leftover.end()) {
        lambda_p = stod(*it);
        flag_lambdap = true;
        leftover.erase(it);
      } else {
        cout << "elas_ui: couldn\'t get value for lambda_p" << endl;
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
        cout << "elas_ui: couldn\'t get value for NP" << endl;
        return -1;
      }
      continue;
    }

    // if (*it == "-1") {
    //   leftover.erase(it);
    //   if (it != leftover.end()) {
    //     x1 = stod(*it);
    //     flag_nlo = true;
    //     leftover.erase(it);
    //   } else {
    //     cout << "elas_ui: couldn\'t get NLO x-value" << endl;
    //     return -1;
    //   }
    //   continue;
    // }

    // if (*it == "-2") {
    //   leftover.erase(it);
    //   if (it != leftover.end()) {
    //     x2 = stod(*it);
    //     flag_n2lo = true;
    //     leftover.erase(it);
    //   } else {
    //     cout << "elas_ui: couldn\'t get N2LO x-value" << endl;
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
    cout << "elas_ui: NP is not ready" << endl;
    return false;
  }
  if (! flag_q0) {
    cout << "elas_ui: q0 is not ready" << endl;
    return false;
  }
  if (! flag_I) {
    cout << "elas_ui: I is not ready" << endl;
    return false;
  }
  if (! flag_Ip) {
    flag_Ip = true;
    I_p = I;
  }
  if (! flag_lambda) {
    cout << "elas_ui: lambda is not ready" << endl;
    return false;
  }
  if (! flag_lambdap) {
    flag_lambdap = true;
    lambda_p = lambda;
  }
  return true;
}

void print_out_paras(void) {

  show_message_fdv("* elasui: cmdline = " + whole_argline, _VBSLVL_LOW_);
  cout << "* NP = " << NP << endl;
  // debug
  // cout << "* NX = " << NX << endl;
  // cout << "* mambda = " << mambda << endl;
  // cout << maxL << endl;
  // cout << maxOpT << endl;
  // cout << alpha_lst_bare[0] << ", " << alpha_lst_bare[5] << endl;
  // cout << pot_name << endl;
  cout << "* lambda = " << lambda << ", lambda_p = " << lambda_p << endl;
  cout << "* I = " << I << ", I_p = " << I_p << endl;
  cout << "* q0 = " << q0 << endl;
}

