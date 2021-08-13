#ifndef THREE_NJ
#define THREE_NJ

#include "iostream"
#include "math.h"
#include "algorithm"
#include "def.h"

void facLabel();

REAL sixJ(REAL j1, REAL j2, REAL j12, REAL j3, REAL j, REAL j23);

REAL nineJ(REAL j1, REAL j2, REAL j12, REAL j3, REAL j4, REAL j34, REAL j13, REAL j24, REAL j);

REAL threeJ(REAL j1, REAL j2, REAL j3, REAL m1=0, REAL m2=0, REAL m3=0);

REAL CGcoef(double j1, double j2,double j3, double m1=0, double m2=0, double m3=0);

REAL threeJm0(REAL j1, REAL j2, REAL j3);

REAL CGm0(double j1, double j2,double j3);

#endif /* THREE_NJ */
