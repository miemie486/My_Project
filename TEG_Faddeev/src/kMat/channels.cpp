/* Input channels, please read ../../notes/channels.pdf for the detail
 */
#include "channels.h"
#include <iostream>

// Obsolete
// void iniChannel(REAL **channels, int row){ // Traditional 6 channel in the Elster's note
//     /* Initialize typical 5 channels

//        Input:
//        channels: array with size row
//        row: number of rows of channels

//        Output:
//        channels: channels with crrect order
//     */

//     //Plese input the channels here, with routine:
//     //l,  s,  j,  lambda, J,    T,   Coupled Channel
//     //e.g.:
//     //l,  s,  j,  lambda, J,    T    Coupled Channel (uncoupled)
//     //0,   0,  0,  0,      .5,   1.    -1
//     // if coupled channel >= 0, l_bar = coupled channel. if coupled channel = -1, this channel is uncoupled
//     REAL channelsTemp[6][_N_QMN_CHN_] {
//       {0,   0,  0,  0,      .5,   1., -1, 0.5},
//         {0,   1,  1,  0,      .5,   0., 2, 0.5},
//           {2,   1,  1,  0,      .5,   0., 0, 0.5},
//             {0,   1,  1,  2,      1.5,  0., 2, 0.5},
//               {2,   1,  1,  2,      1.5,  0., 0, 0.5},
//                 {1,  1,  0,  1,      .5,  1., -1, 0.5}
//     };

//     for(int i = 0; i < row; i++){
//       for(int j = 0; j < _N_QMN_CHN_; j++){
//         channels[i][j] = channelsTemp[i][j];
//       }
//     }
// }

// If two-body subsystems are coupled, return true; else return false
bool is2BChnCoupled(REAL* alpha1, REAL* alpha2) {

  REAL local_eps {1.0e-6};

  if (abs(alpha1[_CHNIDX_lcpld_] + 1.0) < local_eps
    || abs(alpha2[_CHNIDX_lcpld_] + 1.0) < local_eps) {
    return false;
  } else {
    if (abs(alpha1[_CHNIDX_l_] - alpha2[_CHNIDX_lcpld_]) > local_eps && abs(alpha2[_CHNIDX_l_] - alpha1[_CHNIDX_lcpld_]) > local_eps ) {
      return false;
    } else {
      if (abs(alpha1[_CHNIDX_s_]-alpha2[_CHNIDX_s_]) < local_eps
          && abs(alpha1[_CHNIDX_j_]-alpha2[_CHNIDX_j_]) < local_eps
          && abs(alpha1[_CHNIDX_lmbd_]-alpha2[_CHNIDX_lmbd_]) < local_eps
          && abs(alpha1[_CHNIDX_I_]-alpha2[_CHNIDX_I_]) < local_eps
          && abs(alpha1[_CHNIDX_T_]-alpha2[_CHNIDX_T_]) < local_eps
          && abs(alpha1[_CHNIDX_opT_]-alpha2[_CHNIDX_opT_]) < local_eps
          && abs(alpha1[_CHNIDX_pt_]-alpha2[_CHNIDX_pt_]) < local_eps
          && abs(alpha1[_CHNIDX_opJ_]-alpha2[_CHNIDX_opJ_]) < local_eps) {
        return true;
      } else {
        return false;
      }
    }
  }
}

void copy_3Nchn(REAL* alpha_src, REAL* alpha_dest) {
  for (int i = 0; i < _N_QMN_CHN_; i++) {
    alpha_dest[i] = alpha_src[i];
  }
}

void gen_alphaNd_3Nchn(REAL* alpha, REAL lambda, REAL I, REAL opJ) {

  alpha[_CHNIDX_l_]    = 0.0;
  alpha[_CHNIDX_s_]    = 1.0;  // pair spin of 3S1 - 3D1
  alpha[_CHNIDX_j_]    = _SPIN_D_;
  alpha[_CHNIDX_lmbd_] = lambda;
  alpha[_CHNIDX_I_]    = I;
  alpha[_CHNIDX_T_]    = _ISOSPIN_D_;
  alpha[_CHNIDX_lcpld_] = 2.0;
  alpha[_CHNIDX_opT_]  = 0.5;
  alpha[_CHNIDX_pt_]   = ((int) lambda % 2 == 0 ? 0.0 : 1.0) ;
  alpha[_CHNIDX_opJ_]  = opJ;

};

REAL* alloc_3Nchn(void) {

  REAL *alpha_local = new REAL [_N_QMN_CHN_];
  return alpha_local;
}

REAL** alloc_3Nchn_table(int n_channels) {

  REAL **channels_local = new REAL* [n_channels];

  for(int i = 0; i < n_channels; i++) {
    channels_local[i] = alloc_3Nchn();
  }

  return channels_local;

}

void copy_3Nchn_table(int n_channels, REAL** channels_src, REAL** channels_dest) {

  for(int i = 0; i < n_channels; i++){
    copy_3Nchn(channels_src[i], channels_dest[i]);
  }
}

void release_3Nchn(REAL* alpha) {

  if (alpha) {
    delete[] alpha;
  }
}

void dealloc_3Nchn_table(int n_channels, REAL **channels) {

  if (channels) {
    for(int i = 0; i < n_channels; i++){
      release_3Nchn(channels[i]);
    };
    delete[] channels;
  }
}

Chs::~Chs(){
}

bool Chs::testTri(double a, double b, double c){
  /* Test if three numbers obey the triangular condition

     in:
     a, b, c: Quantum number a, b and c

     out:
     return: 1 -- they obey. 0 -- not.
  */
  if(a<= b+c && b<= a+c && c<= a+b && a+b+c - int(a+b+c) < _3NJ_EPS_)
    return 1;
  else
    return 0;
}

void Chs::genChs(){
  /* Generate all possible channels.
   */
  double maxS = 1, maxJ = maxL + 1, maxI = opJ + maxJ, maxLambda = maxI + 0.5, maxT = 1;
  double pt, coupling;
  //int rowRaw = (maxS + maxJ + maxL + maxI + maxLambda + maxT + maxOpT) / 0.5;
  //std::vector<vector<double>> channels;

  channels.clear();
  channelsLong.clear();

  // Load legal quantum numbers
  // Too ugly here, you should rewrite them as vector
  for(double s = maxS; s >= 0; s -= 1)
    for(double l = maxL; l >= 0; l -= 1)
      for(double j = maxJ; j >= 0; j -= 1)
        if(testTri(s, l, j))
          for(double lambda = maxLambda; lambda >= 0; lambda -= 1)
            if(int(lambda + l) % 2 == parity)
              for(double I = maxI; I >= 0; I -= 1)
                if(testTri(lambda, 0.5, I) && testTri(j, I, opJ))
                  for(double T = maxT; T >= 0; T -= 1)
                    if(int(s + l + T) % 2 == 1)
                      for(double opT = maxOpT; opT >= 0; opT -= 1){
                        pt = (int(lambda + l) % 2);
                        coupling = cpCoe(s, l, j);

                        /* Order of push_back conforms to the table defined
                        in channels.h */
                        channels.push_back(l);        //  0
                        channels.push_back(s);        //  1
                        channels.push_back(j);        //  2
                        channels.push_back(lambda);   //  3
                        channels.push_back(I);        //  4
                        channels.push_back(T);        //  5
                        channels.push_back(coupling); //  6
                        channels.push_back(opT);      //  7
                        channels.push_back(pt);       //  8
                        channels.push_back(opJ);      //  9

                        channelsLong.push_back(l);
                        channelsLong.push_back(s);
                        channelsLong.push_back(j);
                        channelsLong.push_back(lambda);
                        channelsLong.push_back(I);
                        channelsLong.push_back(T);
                        channelsLong.push_back(opT);
                        channelsLong.push_back(coupling);
                        channelsLong.push_back(pt);
                        channelsLong.push_back(opJ);      //  9
                      }
  row = channels.size() / col;
  // Find the position of all channels
  int seq[row] {0};
  double cpare;
  for(int idx1 = 0; idx1 < row; idx1++) // First channel
    for(int idx2 = 0; idx2 < row; idx2++) // Second channel
      for(int jdx = 0; jdx < col - 1; jdx++){ // Comparation. The final col is coupling state
        cpare = channelsLong[idx1 * col + jdx] - channelsLong[idx2 * col + jdx]; // Compare
        if(cpare > 0){ // If first channel is bigger, seq += 1 and break
          seq[idx1] += 1;
          break;
        }
        else if(cpare < 0) // If smaller, break
          break;
      }
  // Ordering, long's channels
  std::vector<double> channelsTemp {channelsLong};
  for(int i = 0; i < row; i++)
    for(int j = 0; j < col; j++){
      channelsLong[seq[i] * col + j] = channelsTemp[i * col + j];
    }
  // Ordering, my channels, making which has same sequence as long's channel.
  channelsTemp = channels;
  for(int i = 0; i < row; i++)
    for(int j = 0; j < col; j++){
      channels[seq[i] * col + j] = channelsTemp[i * col + j];
    }
}

void Chs::showChannels(std::string method){
  int idx = 0;

  std::cout << std::endl;
  std::cout << "J^pi = " << opJ << "^" << parity << "\n" << row << " chennels in total\n";

  std::vector<double> channelsTemp;
  if(method == "Long"){
    //printf("l\ts\tj\tlambda\tI\tT\topT\tCoupling\tParity\n");
    printf("%-10s%-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s\n",
      "alpha", "l", "s", "j", "lmbd", "I", "T", "opT", "lcpld", "Pi", "opJ");
    channelsTemp = channelsLong;
  }
  else{
    //printf("l\ts\tj\tlambda\tI\tT\tCoupling\topT\tParity\n");
    printf("%-10s%-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s%-7s\n",
      "alpha", "l", "s", "j", "lmbd", "I", "T",  "lp", "opT", "Pi", "opJ");
    channelsTemp = channels;
  }
  for(int i = 0; i < row; i++){
    printf("%-10d", i);
    for(int j = 0; j < col; j++){
      printf("%-7.1f", channelsTemp[idx]);
      idx += 1;
    }
    printf("\n");
  }
}

void Chs::showChannels(std::vector<int> &alpha_array) {

  printf("%-8s%-6s%-6s%-6s%-6s%-6s%-6s%-6s%-6s%-6s%-6s\n",
    "alpha", "l", "s", "j", "lmbd", "I", "T", "lp", "opT", "Pi", "opJ");

  for (unsigned i = 0; i < alpha_array.size(); i++){
    int alpha = alpha_array[i];
    printf(" %4d   ", alpha);
    for(int j = 0; j < col; j++){
      printf("%-6.1f", channels[alpha * col + j]);
    }
    printf("\n");
  }
}

inline double cpCoe(double& s,double& l, double& j){
  if( (int(j) != int(l)) && (int(j) != 0) )
    return (int(l) < int(j)) ? l+2 : l-2;
  else
    return -1;
}

void Chs::outChannels(double **outChs, int nChannels, int *alphaChs){
  // Output the existed channels to a two dimentional array, with size nChannels * 8. Channel label alpha 0 means the first row of channels.
  for(int i = 0; i < nChannels; i++){
    if(i > row - 1)
      printf("The required alpha label exceeds the exited one, please try to increase the maximum l\n");
    for(int j = 0; j < _N_QMN_CHN_; j++){
      outChs[i][j] = channels[alphaChs[i] * col + j];
    }
  }
}

void Chs::outChannels(double **outChs, int nChannels, std::vector<int> &alpha_array) {

  int *alphaChs = new int[nChannels];
  for (int i = 0; i < nChannels; i++) {
    alphaChs[i] = (int) alpha_array[i];
  }
  this->outChannels(outChs, nChannels, alphaChs);
  delete[] alphaChs;
}
