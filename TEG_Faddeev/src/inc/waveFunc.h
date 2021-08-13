#ifndef WAVEFUNC
#define WAVEFUNC

#include "kMatGen.h"
#include "mkl.h"

class waveFunc{
  /* Solve the waveFunc.
     wave(nChannels * pMesh * qMesh): wave function with different channel.
  */
 protected:
  double *wave;
 public:
  //double wavePart(int idxChannel, double p, double q);
  void genWave(kMatrix);
}
