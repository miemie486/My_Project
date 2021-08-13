/* Find the binding energy of a triton. */

#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>
#include "math.h"
#include <string.h>
#include <vector>
#include <stdio.h>

#include "def.h"
#include "kMatGen.h"
#include "channels.h"

using namespace std;

void parse_line_to_reals(string &str_nums, vector<REAL> &num_array, char delim);
void parse_line_to_complex(string &str_nums, vector<COMPLEX> &num_array);
void read_a_real(REAL& num);
void gen_grid(REAL _xup, REAL _xdn, REAL _yup, REAL _ydn, int nx, int ny,
  vector<COMPLEX> &num_array);


int main(int argc, char* argv[]) {

  string fname, default_fname {"showdetk.in"};
  // bool verbose_requested = false;
  // if (argc < 2) {
  //   fname = default_fname;
  // } else {
  //   if (argv[1][0] == '\0') {
  //     fname = default_fname;
  //   } else {
  //     if (strcmp(argv[1], "-v") == 0) {
  //       verbose_requested = true;
  //       if (argv[2][0] == '\0') {
  //         fname = default_fname;
  //       } else {
  //         fname = argv[2];
  //       }
  //     } else {
  //       fname = argv[1];
  //     }
  //   }
  // }

  REAL maxL, maxOpT, opJ, parity, mambda, cmratio, phi;
  size_t NP, NX;

  // const int maxlength_line = 10240;
  // char cmtline[maxlength_line];
  string cmtstr, pot_name;

  getline(cin, cmtstr);
  scanf("%lf\n", &maxL);

  getline(cin, cmtstr);
  scanf("%lf\n", &maxOpT);

  getline(cin, cmtstr);
  scanf("%lf\n", &opJ);

  getline(cin, cmtstr);
  scanf("%lf\n", &parity);

  // read channel labels
  string labelstr;
  vector<REAL> chnlable_lst;
  getline(cin, cmtstr);
  getline(cin, labelstr);
  parse_line_to_reals(labelstr, chnlable_lst, ',');

  getline(cin, cmtstr);
  getline(cin, pot_name);

  getline(cin, cmtstr);
  scanf("%lf\n", &mambda);

  REAL npr;
  getline(cin, cmtstr);
  scanf("%lf\n", &npr);
  NP = (int) npr;

  getline(cin, cmtstr);
  scanf("%lf\n", &npr);
  NX = (int) npr;

  getline(cin, cmtstr);
  scanf("%lf\n", &cmratio);

  getline(cin, cmtstr);
  scanf("%lf\n", &phi);

  vector<COMPLEX> Elst;
  vector<REAL> cordn;
  getline(cin, cmtstr);
  getline(cin, cmtstr);
  parse_line_to_reals(cmtstr, cordn, ',');

  // debug
  // cout << cmtstr << endl;
  // cout << cordn.size() << endl;
  // cout << cordn[0] << cordn[1] << cordn[2] << cordn[3] << cordn[4] << cordn[5]
  //   << endl;

  gen_grid(cordn[0], cordn[1], cordn[2], cordn[3], (int) cordn[4],
    (int) cordn[5], Elst);

  // // Print out input parameters
  // cout << endl << "/////////////////////////////////" << endl;
  cout << "maxL = " << maxL << endl;
  cout << "maxOpT = " << maxOpT << endl;
  cout << "opJ = " << opJ << endl;
  cout << "parity = " << parity << endl;
  for (size_t i = 0; i < chnlable_lst.size()-1; i++) {
    cout << chnlable_lst[i] << ", ";
  }
  // cout << chnlable_lst[chnlable_lst.size()-1] << endl;
  cout << chnlable_lst.back() << endl;
  cout << pot_name << endl;
  cout << mambda << endl;
  cout << "NP, NX = " << NP << ", " << NX << endl;
  cout << "cmratio = " << cmratio << endl;
  cout << "phi = "<< phi << endl;

  cout << "Number of E's = " << Elst.size() << endl;

  for (size_t i = 0; i < Elst.size()-1; i++) {
    cout << Elst[i] << " ";
  }
  cout << Elst.back() << endl;

  // Generate large channel table allowed by MaxL and MaxOpT
  Chs chsTable {};
  chsTable.writeMaxL(maxL);
  chsTable.writeMaxOpT(maxOpT);
  chsTable.writeOpJ(opJ);
  chsTable.writeParity(parity);
  chsTable.genChs();

  // Initialize the selected channels
  size_t nChannels = chnlable_lst.size();
  int labelChs[nChannels];
  for (size_t i = 0; i < nChannels; i++) {
    labelChs[i] = chnlable_lst[i];
  }
  REAL **channels = new REAL*[nChannels];
  for(size_t i = 0; i < nChannels; i++) {
    channels[i] = new REAL[_N_QMN_CHN_];
  }
  chsTable.outChannels(channels, nChannels, labelChs);

  // Print out the channels the user have seleted
  cout << endl << "/////////////////////////////////" << endl;
  cout << "You have selected channels:\n";
  printf("%-8s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s\n",
    "alpha", "l", "s", "j", "lmbd", "I", "T",  "lp", "opT", "Pi", "opJ");
  for(size_t i = 0; i < nChannels; i++) {
    printf(" %4d   ", labelChs[i]);
    for(int j = 0; j < _N_QMN_CHN_; j++)
      printf("%-5.1f", channels[i][j]);
    printf("\n");
  }

  REAL _CTanMesh = mambda * cmratio;
  // , local_eps {1.0E-6};
  int uptoQn {0}, numPara {1};
  REAL potpara[NUM_PARA_TWOBODY] {0.0};
  potpara[0] = uptoQn; potpara[1] = 0.0;

  time_t now = time(0);
  cout << "==>>>> Job starts at " << ctime(&now) << endl << flush;

  // Initialize chengdu_nnscat global functions/variables
  init_nnscat_chengdu_C();
  // Initialize TEG_Faddeev global functions/variables
  facLabel();

  // Complex-valued pots contains "cmplx"
  if ( pot_name.find("cmplx") != string::npos ) {

    kMatrix<COMPLEX> kMatR {nChannels, _CTanMesh, channels, opJ};
    const COMPLEX ii(0.0, 1.0);
    kMatR.setFacContour(exp(-ii*phi/180.0*3.141593));
    kMatR.setPotName(pot_name); // potential name
    kMatR.loadTwoBodyPara(numPara, potpara);
    kMatR.setMass(_AvgMassN_);
    kMatR.setNX(NX);
    kMatR.setNP(NP);
    kMatR.setNQ(NP);
    kMatR.genGTable();
    kMatR.allocKMat();
    kMatR.initKMatUnit();
    kMatR.genIntzTable_Hmgns();
    for (size_t i = 0; i < Elst.size(); i++) {
      printf("--> E = (%-17.10e,  %-17.10e)\n", Elst[i].real(), Elst[i].imag());
      fflush(stdout);
      COMPLEX detK = kMatR.findKHmgnsDetFunc(Elst[i]);
      printf("====> detK = (%-17.10e,  %-17.10e)\n", detK.real(), detK.imag());
      now = time(0);
      cout << _TIMER_COL_ << ctime(&now);
      fflush(stdout);
    }
  } else {
    cout << pot_name << " is not complex " << endl;
  }

  for(size_t i = 0; i < nChannels; i++){
    delete[] channels[i];
  }
  delete[] channels;

}

void parse_line_to_reals(string &str_nums, vector<REAL> &num_array, char delim) {

  istringstream ss(str_nums);
  string token;
  while(getline(ss, token, delim)) {
    num_array.push_back(stod(token));
  }
}

void gen_grid(REAL _xup, REAL _xdn, REAL _yup, REAL _ydn, int nx, int ny,
  vector<COMPLEX> &num_array) {

  REAL xup {_xup}, xdn {_xdn}, yup {_yup}, ydn {_ydn}, xstep, ystep;
  if (_xup > _xdn) {
    xup = _xup; xdn = _xdn;
  } else {
    xup = _xdn; xdn = _xup;
  }
  if (_yup > _ydn) {
    yup = _yup; ydn = _ydn;
  } else {
    yup = _ydn; ydn = _yup;
  }

  xstep = (xup - xdn)/(nx-1);
  ystep = (yup - ydn)/(ny-1);
  int iter {0};

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      num_array.push_back(COMPLEX(xdn + xstep*i, ydn + ystep*j));
      iter++;
    }
  }
}

// void parse_line_to_complex(string &str_nums, vector<COMPLEX> &num_array) {

//   istringstream ss_line(str_nums);
//   string token;
//   while(getline(ss_line, token, ' ')) {
//     istringstream ss_num(token);
//     COMPLEX c;
//     ss_num >> c;
//     num_array.push_back(c);
//   }
// }

// void read_a_real(REAL& num) {

//   cin >> num;
//   string ignoredline;
//   getline(cin, ignoredline);
//   getline(cin, ignoredline);

// }

