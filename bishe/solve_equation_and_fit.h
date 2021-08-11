
#include<iostream>
#include<cmath>
#include<fstream>

class Caculate_Q
{
  private:
  double miu0=4*M_PI*1e-7; //Vacuum permeability 
  double H0=27500.0;  //surface magnatic field intensity
  double N=1e7;   //the number of the layer
  double d=3.0e-3;   //the thickness of the caveity  
  double Tc=9.23; //critical temperature  
  double rou_ph=1.45e-7; //rou_ph(300K)
  double frequency=1.3e9; //frequency
  double R_res=1.0e-9; //Residual resistance
  //double alpha=2.7*10e-9; //RRR - mean free path conversion
  double kB=1.38065e-23; //Boltzmann constant
  double Delta_0=3.53*kB*Tc/2; //superconducting energy gap at 0K
  double hbar=1.0546e-34; 
  double omega=2*M_PI*frequency;

  //parameters that need to be fitted

  double thickness_a=1.0e-7; //Dirty layer thickness
  double lambda_L=23.0e-9; //London's penetration depth
  double Hc=0.2/miu0; //critical field
  double xi_0=40.0e-9; //coherence length
  double l_1=3.0e-9; //Dirty layer's mean free path
  double l_2=20.0e-9; //bulk Nb's mean free path

  double g=270;
  double beta=1.0e-6;


}