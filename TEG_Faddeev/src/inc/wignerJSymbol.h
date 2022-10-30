#ifndef WIGNERJ
#define WIGNERJ

#include "iostream"
#include "math.h"
#include "algorithm"
#include "mathFunc.h"

class wignerJSymbol {
 private:
  long facArr[21] {0};
 public:
  // Constructor
  wignerJSymbol(){
    // Generating factals from 0! to 20!
    facArr[0]=1;
    for(short i=1; i<21; i++)
      facArr[i]=facArr[i-1]*i;
  }
  // Rule of three. This is just a collection of functions. Just remember don't use coping and =.

  long facto(auto i);
  bool ifTri(double, double, double);
  double triCoe(double, double, double);
  double threeJ(double j1, double j2, double j3, double m1, double m2, double m3);
  //  double sixJ(double j1, double j2, double j12, double j3, double j, double j23){return sixJ(j1, j2, j12, j3, j, j23);};
  //  double nineJ(double j1, double j2, double j12, double j3, double j4, double j34, double j13, double j24, double j) {return nineJ( j1,  j2,  j12,  j3,  j4,  j34,  j13,  j24,  j);}
  double CG(double j1, double j2, double j3, double m1, double m2, double m3);
};

#endif /* WIGNERJ */
