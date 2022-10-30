#include<iostream>
#include<cmath>

using namespace std;

double miu0=4*M_PI*1e-7; //Vacuum permeability 
double H0=27500.0;  //surface magnatic field intensity
double N=1e7;   //the number of the layer
double d=3.0e-3;   //the thickness of the caveity
      double lambda_L=23.0e-9; //London's penetration depth
double Tc=9.23; //critical temperature
      double Hc=0.2/miu0; //critical field
      double xi_0=40.0e-9; //coherence length
double rou_ph=1.45e-7; //rou_ph(300K)
double frequency=1.3e9; //frequency
double l_1=3.0e-9; //Dirty layer's mean free path
double l_2=20.0e-9; //bulk Nb's mean free path
double R_res=1.0e-9; //Residual resistance
      double thickness_a=1.0e-7; //Dirty layer thickness
//double alpha=2.7*10e-9; //RRR - mean free path conversion
double kB=1.38065e-23; //Boltzmann constant
double Delta_0=3.53*kB*Tc/2; //superconducting energy gap at 0K
double hbar=1.0546e-34; 
double omega=2*M_PI*frequency;


double g=270;
double beta=1.0e-6;



double Hc_T(double T)
{
  double a;
  a=Hc*(1-pow(T/Tc,2));
  return a;
}

double Delta(double T)  //superconducting energy gap at T
{
  double a;
  a=Delta_0*sqrt(cos(M_PI_2*(T/Tc)*(T/Tc)));
  return a;
}

double lambda_L_T(double T)
{
  double a;
  a=lambda_L*pow((1-pow(T/Tc,4)),-0.5);
  return a;
}

double lambda_L_T_H(double T,double H)
{
  double a;
  a=lambda_L_T(T)*sqrt(2.0)/(sqrt(1+sqrt(1-pow(H/Hc_T(T),2))));
  return a;
}

double alpha(double T)
{
  double a;
  a=omega*miu0*xi_0*lambda_L*lambda_L/rou_ph*(M_PI*Delta(T)/(hbar*omega)*tanh(Delta(T)/(2*kB*T)));
  return a=2.7e-9;
}

double C(double T)
{
  double a;
  a=0.5*(Delta(T)/(kB*T))*log(Delta(T)/(hbar*omega))*(pow(omega,2)*pow(miu0,2)*pow(lambda_L,3))/(alpha(T)*rou_ph);
  /*
  cout<< Delta(T)/(kB*T)<<endl;
  cout<< log(Delta(T)/(hbar*omega))<<endl;
  cout<< pow(omega,2)*pow(miu0,2)*pow(lambda_L,3)<<endl;
  cout<< alpha(T)*rou_ph<<endl;
  */
  return a;
}

double Rs1(double T)  //surface resistence which mfp=l1
{
  double a;
  a=C(T)*l_1*pow((1+xi_0/l_1),1.5)*exp(-Delta(T)/(kB*T));
  //cout<<exp(-Delta(T)/(kB*T))<<endl;
  //cout<<C(T)<<" "<< exp(-Delta(T)/(kB*T))<<endl;
  return a+R_res;
}

double Rs2(double T)  //surface resistence which mfp=l2
{
  double a;
  a=C(T)*l_2*pow((1+xi_0/l_2),1.5)*exp(-Delta(T)/(kB*T));
  return a+R_res;
}

double H(double z)  //magnatic field intensity
{
  double a;
  a=H0*exp(-z/lambda_L);
  return a;
}

double P(double T, int k)    //dissipate power
{
  double a;
  if(k*d/N<=thickness_a)
  {
    a=0.5*(Rs1(T)+R_res)*(H((k-1)*d/N)*H((k-1)*d/N)-H(k*d/N)*H(k*d/N));
    cout<<Rs1(T)<<endl;
  }
  else
  {
    a=0.5*(Rs2(T)+R_res)*(H((k-1)*d/N)*H((k-1)*d/N)-H(k*d/N)*H(k*d/N));
  }
  return a;
}


double Q(double H,double T)
{
  double a;
  a=(g/(Rs1(T)))+(g/Rs2(T)-g/Rs1(T))*exp(-thickness_a/(lambda_L+beta*(H/Hc_T(T))));
  return a;
}

double deltaQ(double H,double T)
{
  double a;
  //a=(g/Rs2(T)-g/Rs1(T))*exp(-thickness_a/(lambda_L+beta*(H/Hc_T(T))));
  a=(g/Rs2(T)-g/Rs1(T));
  return a;
}

double asd(double H,double T)
{
  double a;
  a=exp(-thickness_a/(lambda_L+beta*(H/Hc_T(T))));
  return a;
}

//test code
/*
double Rs(double T,double l)  //surface resistence which mfp=l1
{
  double a;
  a=C(T)*l*pow((1+xi_0/l),1.5)*exp(-Delta(T)/(kB*T));
  //cout<<exp(-Delta(T)/(kB*T))<<endl;
  //cout<<C(T)<<" "<< exp(-Delta(T)/(kB*T))<<endl;
  return a+R_res;
}

*/

int main()
{
  
  for(int i=0;i<=100;i++)
  {
    cout<<i*Hc_T(2)*1000.0*miu0/100.0<<" "<<Q(i*Hc_T(2)/100.0,2)<<" "<<deltaQ(i*Hc_T(2)/100.0,2)<<endl;
  }
  
}