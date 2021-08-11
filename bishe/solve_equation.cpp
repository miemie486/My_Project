#include<iostream>
#include<cmath>
#include<ctime>
#include<assert.h>
#include<fstream>
#include<omp.h>

using namespace std;

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
double Delta_0=3.84*kB*Tc/2; //superconducting energy gap at 0K
double hbar=1.0546e-34; 
double omega=2*M_PI*frequency;


double thickness_a=1.0e-8; //Dirty layer thickness
double lambda_L=32.0e-9; //London's penetration depth
double Hc_1=0.2/miu0; //critical field 1
double Hc_2=0.2/miu0; //critical field 2
double xi_0=62.0e-9/M_PI*2; //coherence length
double l_1=3.0e-9; //Dirty layer's mean free path
double l_2=50.0e-9; //bulk Nb's mean free path


double g=270;
double beta=1.0e-6;




double Hc_T_1(double T)
{
  double a;
  a=Hc_1*(1-pow(T/Tc,2));
  return a;
}

double Hc_T_2(double T)
{
  double a;
  a=Hc_2*(1-pow(T/Tc,2));
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

double lambda_L_T_H_1(double T,double H)
{
  double a;
  a=lambda_L_T(T)*sqrt(2.0)/(sqrt(1+sqrt(1-pow(H/Hc_T_1(T),2))));
  return a;
}

double lambda_L_T_H_2(double T,double H)
{
  double a;
  a=lambda_L_T(T)*sqrt(2.0)/(sqrt(1+sqrt(1-pow(H/Hc_T_2(T),2))));
  return a;
}


double alpha_1(double T,double H)
{
  double a;
  a=omega*miu0*xi_0*lambda_L_T_H_1(T,H)*lambda_L_T_H_1(T,H)/rou_ph*(M_PI*Delta(T)/(hbar*omega)*tanh(Delta(T)/(2*kB*T)));
  return a;
}

double alpha_2(double T,double H)
{
  double a;
  a=omega*miu0*xi_0*lambda_L_T_H_2(T,H)*lambda_L_T_H_2(T,H)/rou_ph*(M_PI*Delta(T)/(hbar*omega)*tanh(Delta(T)/(2*kB*T)));
  return a;
}

double C_dirty_1(double T,double H)
{
  double a;
  a=0.5*(Delta(T)/(kB*T))*log(Delta(T)/(hbar*omega))*(pow(omega,2)*pow(miu0,2)*pow(lambda_L_T_H_1(T,H),3))/(alpha_1(T,H)*rou_ph);
  /*
  cout<< Delta(T)/(kB*T)<<endl;
  cout<< log(Delta(T)/(hbar*omega))<<endl;
  cout<< pow(omega,2)*pow(miu0,2)*pow(lambda_L,3)<<endl;
  cout<< alpha(T)*rou_ph<<endl;
  */
  return a;
}

double C_dirty_2(double T,double H)
{
  double a;
  a=0.5*(Delta(T)/(kB*T))*log(Delta(T)/(hbar*omega))*(pow(omega,2)*pow(miu0,2)*pow(lambda_L_T_H_2(T,H),3))/(alpha_2(T,H)*rou_ph);
  /*
  cout<< Delta(T)/(kB*T)<<endl;
  cout<< log(Delta(T)/(hbar*omega))<<endl;
  cout<< pow(omega,2)*pow(miu0,2)*pow(lambda_L,3)<<endl;
  cout<< alpha(T)*rou_ph<<endl;
  */
  return a;
}


double Rs1(double T,double H)  //surface resistence which mfp=l1
{
  double a;
  a=C_dirty_1(T,H)*l_1*pow((1+xi_0/l_1),1.5)*exp(-Delta(T)/(kB*T));
  //cout<<exp(-Delta(T)/(kB*T))<<endl;
  //cout<<C(T)<<" "<< exp(-Delta(T)/(kB*T))<<endl;
  return a+R_res;
}

double Rs2(double T,double H)  //surface resistence which mfp=l2
{
  double a;
  a=C_dirty_2(T,H)*l_2*pow((1+xi_0/l_2),1.5)*exp(-Delta(T)/(kB*T));
  return a+R_res;
}

double H(double z)  //magnatic field intensity
{
  double a;
  a=H0*exp(-z/lambda_L);
  return a;
}

double Q(double H,double T)
{
  double a;
  //cout<<thickness_a/lambda_L_T_H(T,H)<<" "<<thickness_a<<" "<<lambda_L_T_H(T,H)<<" "<<H/Hc_T(T)<<"      ";
  //a=(g/(Rs1(T,H)))+(g/Rs2(T,H)-g/Rs1(T,H))*exp(-thickness_a/lambda_L_T_H(T,H));
  //cout<<g/Rs1(T,H)<<endl;
  //a=(g/(Rs1(T,H)))+(g/Rs2(T,H)-g/Rs1(T,H))*exp(-thickness_a/lambda_L_T_H(T,H));
  //a=g/(Rs1(T,H)*(1-exp(-2*thickness_a/(lambda_L_T_H(T,H)*(1+xi_0/l_1))))+Rs2(T,H)*exp(-2*thickness_a/(lambda_L_T_H(T,H)*(1+xi_0/l_1))));
  a=g/(Rs1(T,0)*(1-exp(-2*thickness_a/(lambda_L_T_H_1(T,H)*(1+xi_0/l_1))))+Rs2(T,0)*exp(-2*thickness_a/(lambda_L_T_H_1(T,H)*(1+xi_0/l_1))));
  return a;
}

double fit(double length,double (*H_and_Q)[2],double T)
{
  double a=0;
  for(int i=0;i<length;i++)
  {
    a=a+pow((H_and_Q[i][1]-Q(H_and_Q[i][0],T))/1.0e10,2);
  }
  return a;
}

int main()
{
  
  for(int i=0;i<=100;i++)
  {
    //cout<<i*Hc_T(2)*1000.0*miu0/100.0<<" "<<Q(i*Hc_T(2)/100.0,2)<<" "<< asd(i*Hc_T(2)/100.0,2)<<endl;
    cout<<i*Hc_T_1(2)*miu0/100*1000<<" "<<Q(i*Hc_T_1(2)/100.0,2)<<endl;
  }
  cout<<Rs1(4,0)<<" "<<Rs2(4,0)<<" "<<(1+xi_0/l_2)*lambda_L_T_H_1(4,0)<<endl;
}