

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


double thickness_a=1.0e-7; //Dirty layer thickness
double lambda_L=23.0e-9; //London's penetration depth
double Hc=0.2/miu0; //critical field
double xi_0=40.0e-9; //coherence length
double l_1=3.0e-9; //Dirty layer's mean free path
double l_2=20.0e-9; //bulk Nb's mean free path

double x1=1.0;
double x2=1.0;


double g=270;
double beta=1.0e-6;

const double Maxiter=10.0;
const int dimension=8;
const double uppper_limit[dimension]={50.0e-9,1.0e-7,0.2/miu0,100.0e-9,10.0e-9,100.0e-9,100.0,100.0};
const double lower_limit[dimension]={1.0e-9,10.0e-9,0.15/miu0,10.0e-9,1.0e-9,10.0e-9,-100.0,-100.0};




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
  a=lambda_L_T(T)*(1+x1*(H/Hc)+x2*pow(H/Hc,2));
  return a;
}

double alpha(double T,double H)
{
  double a;
  a=omega*miu0*xi_0*lambda_L_T_H(T,H)*lambda_L_T_H(T,H)/rou_ph*(M_PI*Delta(T)/(hbar*omega)*tanh(Delta(T)/(2*kB*T)));
  return a;
}

double C_dirty(double T,double H)
{
  double a;
  a=0.5*(Delta(T)/(kB*T))*log(Delta(T)/(hbar*omega))*(pow(omega,2)*pow(miu0,2)*pow(lambda_L_T_H(T,H),3))/(alpha(T,H)*rou_ph);
  /*
  cout<< Delta(T)/(kB*T)<<endl;
  cout<< log(Delta(T)/(hbar*omega))<<endl;
  cout<< pow(omega,2)*pow(miu0,2)*pow(lambda_L,3)<<endl;
  cout<< alpha(T)*rou_ph<<endl;
  */
  return a;
}

double C(double T)
{
  double a;
  a=0.5*(Delta(T)/(kB*T))*log(Delta(T)/(hbar*omega))*(pow(omega,2)*pow(miu0,2)*pow(lambda_L,3))/(2.7e-9*rou_ph);
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
  a=C_dirty(T,H)*l_1*pow((1+xi_0/l_1),1.5)*exp(-Delta(T)/(kB*T));
  //cout<<exp(-Delta(T)/(kB*T))<<endl;
  //cout<<C(T)<<" "<< exp(-Delta(T)/(kB*T))<<endl;
  return a+R_res;
}

double Rs2(double T,double H)  //surface resistence which mfp=l2
{
  double a;
  a=C_dirty(T,H)*l_2*pow((1+xi_0/l_2),1.5)*exp(-Delta(T)/(kB*T));
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
  a=g/(Rs1(T,0)*(1-exp(-2*thickness_a/(lambda_L_T_H(T,H)*(1+xi_0/l_1))))+Rs2(T,0)*exp(-2*thickness_a/(lambda_L_T_H(T,H)*(1+xi_0/l_1))));
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
  //读取文件中的数据
  ifstream file;
  file.open("droping.txt");
  double data[27][2];
  for(int i=0;!file.eof();i++)
  {
    file>>data[i][0]>>data[i][1];
    data[i][0]=data[i][0]*4.25/1000/miu0;
  }
  double best_position[dimension]={10e-8,23.0e-9,0.2/miu0,40.0e-9,3.0e-9,20.0e-9,1.0,1.0};
  double best_fit=1000000;
  double length=27.0;
 
  #pragma omp parallel for
  
    for(int i1=0;i1<10000000;i1++)
    {
      srand(clock());
      thickness_a=lower_limit[0]+(uppper_limit[0]-lower_limit[0])*rand()*1.0/RAND_MAX;
      lambda_L=lower_limit[1]+(uppper_limit[1]-lower_limit[1])*rand()*1.0/RAND_MAX;
      //thickness_a=1000*lambda_L;
      Hc=lower_limit[2]+(uppper_limit[2]-lower_limit[2])*rand()*1.0/RAND_MAX;
      xi_0=lower_limit[3]+(uppper_limit[3]-lower_limit[3])*rand()*1.0/RAND_MAX;
      l_1=lower_limit[4]+(uppper_limit[4]-lower_limit[4])*rand()*1.0/RAND_MAX;
      //l_2=lower_limit[5]+(uppper_limit[5]-lower_limit[5])*rand()*1.0/RAND_MAX;
      l_2=xi_0/2;
      x1=lower_limit[6]+(uppper_limit[6]-lower_limit[6])*rand()*1.0/RAND_MAX;
      x2=lower_limit[7]+(uppper_limit[7]-lower_limit[7])*rand()*1.0/RAND_MAX;
      double fangcha=fit(length,data,2);
      //cout<<fangcha<<endl;
      if(fangcha<best_fit)
      {
        best_fit=fangcha;
        best_position[0]=thickness_a;
        best_position[1]=lambda_L;
        best_position[2]=Hc;
        best_position[3]=xi_0;
        best_position[4]=l_1;
        best_position[5]=l_2;
        best_position[6]=x1;
        best_position[7]=x2;
      }

    }
  

  thickness_a=best_position[0];
  lambda_L=best_position[1];
  Hc=best_position[2];
  xi_0=best_position[3];
  l_1=best_position[4];
  l_2=best_position[5];
  x1=best_position[6];
  x2=best_position[7];
  for(int j=0;j<length;j++)
  {
    cout<<data[j][0]*miu0*1000<<" "<<Q(data[j][0],2)<<endl;//" "<<exp(-2*thickness_a/(lambda_L_T_H(2,data[j][0])*(1+xi_0/l_1)))<<endl;
  }
  cout<<best_position[0]<<" "<<best_position[1]<<" "<<best_position[2]*miu0<<" "<<best_position[3]<<" "<<best_position[4]<<" "<<best_position[5]<<" "<<best_position[6]<<" "<<best_position[7]<<" "<<best_fit<<endl;
  cout<<g/Rs1(2,0)<<" "<<g/Rs2(2,0)<<" "<<exp(-2*thickness_a/(lambda_L_T_H(2,0)*(1+xi_0/l_1)))<<endl;
  return 0;
}