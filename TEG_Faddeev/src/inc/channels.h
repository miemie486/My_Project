#ifndef CHANNELS
#define CHANNELS

#include "def.h"
#include <vector>

/* Rosetta stone of notations for QMN
  Stadler  Elster   codenote      Meaning
   L        l         l          Orbit momentem of the pairs
   T        T         T          Isospin of the pairs
   opT      opT       opT        Total isospin
   l        lambda    lambda     Orbit momentem of the third paticle
   S        s         s          Spin of the pair
   I        j         j          Total momentem of the pair
   j        J         I          Total momentem of the third paticle
                      -1         Coupled channels' pair's orbit angular momentem
                      opT        Total isospin
 */
// notation of codenote is followed here
// alpha[_CHNIDX_###_] is an REAL array that stores quantum numbers of a 3N
// channel
#define _CHNIDX_l_            0
#define _CHNIDX_s_            1
#define _CHNIDX_j_            2
#define _CHNIDX_lmbd_         3
#define _CHNIDX_I_            4
#define _CHNIDX_T_            5
#define _CHNIDX_lcpld_        6 // l of the partner partial wave. -1 if uncoupled
#define _CHNIDX_opT_          7
#define _CHNIDX_pt_           8
#define _CHNIDX_opJ_          9

// void iniChannel(REAL **channels, int row); // Obsolete

bool is2BChnCoupled(REAL* alpha1, REAL* alpha2);
void copy_3Nchn(REAL* alpha_src, REAL* alpha_dest);
// Generate alpha_Nd quantum numbers
void gen_alphaNd_3Nchn(REAL* alpha, REAL lambda, REAL I, REAL opJ);

void copy_3Nchn_table(int n_channels, REAL** channels_src, REAL** channels_dest);
REAL* alloc_3Nchn(void);
void release_3Nchn(REAL* alpha);
REAL** alloc_3Nchn_table(int n_channels);
void dealloc_3Nchn_table(int n_channels, REAL **channels);

inline double cpCoe(double& s,double& l, double& j);

// Chs generates a large table of three-nucleon channel labels allowed by
// {maxL, maxOpT, opJ, parity}
class Chs{
  // It will be better if you have maxJ as well
private:
  std::vector<double> channels;
  std::vector<double> channelsLong;
  int col {_N_QMN_CHN_}, row {1};
  /*
    parity of channls. 0 = even and 1 = odd.
     */
  // bool parity {0};
  bool parity;
  double maxL, maxOpT, opJ;
  // double maxL {2.0}, maxOpT {0.5};
  // double opJ {0.5}; // total angular momentem

public:
  // Constructor
  Chs(){};
  // Chs(bool _parity, double _maxL, double _opJ) : parity(_parity), maxL(_maxL), opJ(_opJ) {};
  // Chs(double _maxL, double _maxOpT, double _opJ, bool _parity) : maxL(_maxL),
  //   maxOpT(_maxOpT), opJ(_opJ), parity(_parity) {};
  Chs(double _maxL, double _maxOpT, double _opJ, bool _parity) {

    maxL = _maxL;
    maxOpT = _maxOpT;
    opJ = _opJ;
    parity = _parity;
  };

  // Rule of three
  ~Chs();
  Chs(const Chs& other){};
  Chs& operator=(const Chs& other){return *this;};

  // void iniChannelsOld(int _rowOld);
  /* initializing the channels with which in function iniChannel. usually five to six channels.

     in:
     _rowOld: how many channels you want. For more info, read comments in function iniChannel

     out:
     channels: two dimentional array
  */
  bool testTri(double a, double b, double c);
  /* Test if three numbers obey the triangular condition

     in:
     a, b, c: Quantum number a, b and c

     out:
     return: 1 -- they obey. 0 -- not.
  */
  double readChs(int _row, int _col, std::string method = "None"){
    /* read the quantum number in new version channels.

       input:
       method: If method is = none, than quantum numbers are stroed as what they are and labeled with the sequence,
       {l, s, j, lambda, I ,T, coupling state, opT}.
       If method is equal to "Long", than quantum is output as integers and with sequence which is described in note/notes.pdf.
    */
    if(method == "Long")
      return channelsLong[_row * col + col];
    else
      return channels[_row * col + col];
  }
  double readMaxL(){return maxL;}
  double readMaxOpT(){return maxOpT;}
  double readOpJ(){return opJ;}
  int readCol(){return col;}
  void writeMaxL(double _maxL){maxL = _maxL;}
  void writeMaxOpT(double _maxOpT){maxOpT = _maxOpT;}
  void writeOpJ(double _opJ){opJ = _opJ;}
  void writeParity(bool _parity){parity = _parity;}
  void genChs();
  /* By default, show channels in the format of Zeyuan*/
  void showChannels(std::string method = "Zeyuan");
  // Print out a subset of Chs according to alpha_array
  void showChannels(std::vector<int> &alpha_array);
  void outChannels(double **outChs, int nChannels, int *alphaChs);
  void outChannels(double **outChs, int nChannels, std::vector<int> &alpha_array);
};

#endif
