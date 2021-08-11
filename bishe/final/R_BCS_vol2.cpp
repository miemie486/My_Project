

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
    double Delta_0=3.70*kB*Tc/2;  //superconducting energy gap at 0K
    double frequency=1.3e9;       //frequency
    double omega=2*M_PI*frequency;
    double g=270;                 //Geometric factor 

  public:
    double thickness_a=100.0e-9;    //Dirty layer thickness
    double lambda_L=32.0e-9;      //London's penetration depth
    double Hc=0.18/miu0;           //critical field
    double xi_0=40.0e-9;          //coherence length
    double l_1=3.0e-9;            //Dirty layer's mean free path
    double l_2=20.0e-9;           //bulk Nb's mean free path
    double x1=0.0;
    double x2=0.0;
    
    double Hc_T(double T){return (Hc*(1-pow(T/Tc,2)));};
    double Delta(double T){return (Delta_0*sqrt(cos(M_PI_2*(T/Tc)*(T/Tc))));};
    double lambda_L_T(double T){return (lambda_L*pow((1-pow(T/Tc,4)),0.5));};
    double lambda_L_T_H(double T, double H, double l);
    double alpha(double T, double H, double l);
    double C_dirty(double T, double H);
    double Rs1(double T, double H);
    double Rs2(double T, double H);
    double R_BCS(double T, double H);
    double Delta_T_H(double T,double H);
};

double two_layer_model::Delta_T_H(double T, double H)
{
  double aa;
  aa=Delta(T)*(1+x1*(H/Hc)+x2*(H/Hc)*(H/Hc));
  return aa;
}

double two_layer_model::lambda_L_T_H(double T, double H,double l)
{
  double aa;
  aa=M_PI*miu0*Delta_T_H(T,H)/(hbar*alpha(T,H,l)*rou_ph/l)*
  tanh(Delta_T_H(T,H)/(2*kB*T));
  return sqrt(1/aa);
}



double two_layer_model::alpha(double T, double H, double l)
{
  double aa;
  /*
  aa=omega*miu0*xi_0*lambda_L_T_H(T,H,l)*lambda_L_T_H(T,H,l)/
  rou_ph*(M_PI*Delta(T)/(hbar*omega)*tanh(Delta(T)/(2*kB*T)));
  */
  return 2.7e-9;
}

double two_layer_model::C_dirty(double T,double H)
{
  double aa;
  aa=0.5*(Delta_T_H(T,H)/(kB*T))*log(Delta_T_H(T,H)/(hbar*omega));
  return aa;
}

double two_layer_model::Rs1(double T, double H)  //surface resistence which mfp=l1
{
  double aa;
  aa=C_dirty(T,H)*pow(omega*miu0*alpha(T,H,l_1)*rou_ph/l_1,0.5)*
  pow(omega*miu0*lambda_L_T_H(T,H,l_1)*lambda_L_T_H(T,H,l_1)/(alpha(T,H,l_1)*rou_ph/l_1),1.5)
  *exp(-Delta_T_H(T,H)/(kB*T));
  return aa+R_res;
}

double two_layer_model::Rs2(double T, double H) 
{
  double aa;
  aa=C_dirty(T,H)*pow(omega*miu0*alpha(T,H,l_2)*rou_ph/l_2,2)*
  pow(omega*miu0*lambda_L_T_H(T,H,l_2)*lambda_L_T_H(T,H,l_2)/(alpha(T,H,l_2)*rou_ph/l_2),1.5)
  *exp(-Delta_T_H(T,H)/(kB*T));
  return aa+R_res;
}

double two_layer_model::R_BCS(double T, double H)
{
  double aa;
  aa=(Rs1(T,H)*(1-exp(-2*thickness_a/(lambda_L_T_H(T,H,l_1))))+
  Rs2(T,H)*exp(-2*thickness_a/(lambda_L_T_H(T,H,l_1))));
  return aa;
}

int main()
{
  two_layer_model TLM;
  double miu0=4*M_PI*1e-7;
  double kB=1.38065e-23; 
  double Tc=9.23;
  /*
  cout<<TLM.Delta_T_H(2,0.1/miu0)/(kB*Tc)<<endl;
  cout<<TLM.lambda_L_T_H(2,0.1/miu0,3.0e-9)<<endl;
  cout<<TLM.lambda_L*pow(1+TLM.xi_0/3.0e-9,0.5)<<endl;
  cout<<TLM.Rs1(2,0)<<endl;
  */
  for(int i=1;i<=100;i++)
  {
    TLM.l_1=1e-9*pow(10,i/10.0);
    cout<<TLM.l_1<<" "<<TLM.Rs1(2,0)<<endl;
  }

}
