#ifndef UMATRIX
#define UMATRIX
 /* TMAT */
#include <stdio.h>
#include "kMatGen.h"
#include <fstream>
#include <ctime>
// #include "timer.h"
#include "fstr.h"
#include "mkl.h"
#endif
class UMAT
{
    private:
    REAL q0;
    REAL q0_p;
    REAL energyLoc;
    REAL m;
    REAL jd=1;
    REAL opJ;
    public:
    std::vector<COMPLEX> S_Matrix;
    std::vector<COMPLEX> delta_Matrix;
    std::vector<COMPLEX> U_Matrix;
    std::vector<COMPLEX> Mixing_Parameters;
    void SetMass(REAL _m){m=_m;};
    void SetEnergyLoc(REAL _energyLoc){energyLoc=_energyLoc;};
    void Setq0(REAL _q0){q0=_q0;};
    void Setq0_p(REAL _q0_p){q0_p=_q0_p;};
    void SetopJ(REAL _opJ){opJ=_opJ;};

    REAL GetMass(){return m;};

    COMPLEX U_sigma(REAL lambda_,REAL Sigma_, REAL lambda, REAL Sigma, string path);

    REAL hat(REAL x);

    REAL delta(REAL i,REAL j);

    COMPLEX S(REAL lambda_,REAL Sigma_, REAL lambda, REAL Sigma, string path);

    void Build_S_Matrix(REAL parity,string path);

    void Build_delta_Matrix_and_U_Matrix(REAL parity,string path);

    void Caculate_Mixing_Parameters(REAL parity,string path);
};
