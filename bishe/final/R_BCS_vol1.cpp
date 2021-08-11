

#include<iostream>
#include<cmath>
#include<ctime>
#include<assert.h>
#include<fstream>
#include<vector>

using namespace std;

class two_layer_model
{
  protected:
    double hbar=1.0546e-34;
    double miu0=4*M_PI*1e-7;      //Vacuum permeability 
    double kB=1.38065e-23;        //Boltzmann constant

    double Tc=9.23;               //critical temperature
    double rou_ph=1.45e-7;        //rou_ph(300K)
    double R_res=1.0e-9;          //Residual resistance
    double Delta_0=3.85*kB*Tc/2;  //superconducting energy gap at 0K
    double frequency=1.3e9;       //frequency
    double omega=2*M_PI*frequency;
    double g=270;                 //Geometric factor 

  public:
    double thickness_a=100.0e-9;    //Dirty layer thickness
    double lambda_L=32.0e-9;      //London's penetration depth
    double Hc=0.2/miu0;           //critical field
    double xi_0=40.0e-9;          //coherence length
    double l_1=3.0e-9;            //Dirty layer's mean free path
    double l_2=20.0e-9;           //bulk Nb's mean free path
    double x1=10.0;
    double x2=10.0;
    
    double Hc_T(double T){return (Hc*(1-pow(T/Tc,2)));};
    double Delta(double T){return (Delta_0*sqrt(cos(M_PI_2*(T/Tc)*(T/Tc))));};
    double lambda_L_T(double T){return (lambda_L*pow((1-pow(T/Tc,4)),0.5));};
    double lambda_L_T_H(double T, double H);
    double alpha(double T, double H);
    double C_dirty(double T, double H);
    double Rs1(double T, double H);
    double Rs2(double T, double H);
    double R_BCS(double T, double H);
    

};

double two_layer_model::lambda_L_T_H(double T, double H)
{
  double aa;
  aa=lambda_L_T(T)*(1+x1*(H/Hc)+x2*pow(H/Hc,2));
  return aa;
}

//在这里，我们假设alpha为定值。

double two_layer_model::alpha(double T, double H)
{
  double aa;
  aa=omega*miu0*xi_0*lambda_L_T_H(T,H)*lambda_L_T_H(T,H)/
  rou_ph*(M_PI*Delta(T)/(hbar*omega)*tanh(Delta(T)/(2*kB*T)));
  return aa;
}

double two_layer_model::C_dirty(double T,double H)
{
  double aa;
  aa=0.5*(Delta(T)/(kB*T))*log(Delta(T)/(hbar*omega))*
  (pow(omega,2)*pow(miu0,2)*pow(lambda_L_T_H(T,H),3))
  /(alpha(T,H)*rou_ph);
  return aa;
}

double two_layer_model::Rs1(double T, double H)  //surface resistence which mfp=l1
{
  double aa;
  aa=C_dirty(T,H)*l_1*pow((1+xi_0/l_1),1.5)*exp(-Delta(T)/(kB*T));
  return aa+R_res;
}

double two_layer_model::Rs2(double T, double H) 
{
  double aa;
  aa=C_dirty(T,H)*l_2*pow((1+xi_0/l_2),1.5)*exp(-Delta(T)/(kB*T));
  return aa+R_res;
}

double two_layer_model::R_BCS(double T, double H)
{
  double aa;
  aa=(Rs1(T,H)*(1-exp(-2*thickness_a/(lambda_L_T_H(T,H)*pow((1+xi_0/l_1),0.5))))+
  Rs2(T,H)*exp(-2*thickness_a/(lambda_L_T_H(T,H)*pow((1+xi_0/l_1),0.5))));
  return aa;
}
/*
int main()
{
  two_layer_model TLM;
  cout<<TLM.Rs1(2,0)<<endl;
}
*/