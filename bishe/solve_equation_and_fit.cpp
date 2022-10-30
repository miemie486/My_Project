

#include<iostream>
#include<cmath>
#include<ctime>
#include<assert.h>
#include<fstream>

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
double Delta_0=3.53*kB*Tc/2; //superconducting energy gap at 0K
double hbar=1.0546e-34; 
double omega=2*M_PI*frequency;


double thickness_a=1.0e-7; //Dirty layer thickness
double lambda_L=23.0e-9; //London's penetration depth
double Hc=0.2/miu0; //critical field
double xi_0=40.0e-9; //coherence length
double l_1=3.0e-9; //Dirty layer's mean free path
double l_2=20.0e-9; //bulk Nb's mean free path


double g=270;
double beta=1.0e-6;

const double C1=2.0,C2=2.0;
const int n=100;
const int Maxiter=1000;
const double Maxv[6]={1.0e-7,1.0e-7,0.01/miu0,1.0e-7,1.0e-9,10.0e-9};
const double Initial_omega=1.4;
const double Final_omega=0.4;
const int dimension=6;
const double uppper_limit[6]={1.0e-8,1.0e-7,0.3/miu0,1.0e-7,10.0e-9,100.0e-9};
const double lower_limit[6]={1.0e-10,10.0e-9,0.15/miu0,10.0e-9,1.0e-9,10.0e-9};

double _velocity[dimension];



double f(double *x)
{
	double a=0;
  for(int i=0;i<dimension;i++)
  {
    a=a+x[i]*sin(x[i]/730.0)+x[i]*cos(x[i]/260.0);
  }
	return a;
}


double* velocity(double *v_before,double *identity_best,double *identity_now,double *globle_best,double omega)
{
	srand(clock());
  for(int i=0;i<dimension;i++)
  {
    _velocity[i]=omega*v_before[i]+C1*rand()*1.0/RAND_MAX*(identity_best[i]-identity_now[i])+C2*rand()*1.0/RAND_MAX*(globle_best[i]-identity_now[i]);
  }
	
	return 0;
}

int Min(double *array)
{
	double min=array[0];
	int a=0;
	for(int i=0;i<n;i++)
	{
		if((min-array[i])<=0){}
		else
		{
			min=array[i];
			a=i;
		}
	}
	return a;
}




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

double alpha(double T,double H)
{
  double a;
  a=omega*miu0*xi_0*lambda_L_T_H(T,H)*lambda_L_T_H(T,H)/rou_ph*(M_PI*Delta(T)/(hbar*omega)*tanh(Delta(T)/(2*kB*T)));
  return a=2.7e-9;
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
  a=C(T)*l_2*pow((1+xi_0/l_2),1.5)*exp(-Delta(T)/(kB*T));
  return a+R_res;
}

double H(double z)  //magnatic field intensity
{
  double a;
  a=H0*exp(-z/lambda_L);
  return a;
}

double Q(double H,double T,double *position)
{
  thickness_a=position[0];
  lambda_L=position[1];
  Hc=position[2];
  xi_0=position[3];
  l_1=position[4];
  l_2=position[5];
  double a;
  //cout<<thickness_a/lambda_L_T_H(T,H)<<" "<<thickness_a<<" "<<lambda_L_T_H(T,H)<<" "<<H/Hc_T(T)<<"      ";
  //a=(g/(Rs1(T,H)))+(g/Rs2(T,H)-g/Rs1(T,H))*exp(-thickness_a/lambda_L_T_H(T,H));
  //cout<<g/Rs1(T,H)<<endl;
  a=(g/(Rs1(T,H)))+(g/Rs2(T,H)-g/Rs1(T,H))*exp(-thickness_a/lambda_L_T_H(T,H));
  return a;
}

double fit(double length,double (*H_and_Q)[2],double T, double *position)
{
  double a=0;
  for(int i=0;i<length;i++)
  {
    a=a+pow((H_and_Q[i][1]-Q(H_and_Q[i][0],T,position))/1.0e10,2);
  }
  return a;
}


int main()
{
  
  ifstream file;
  file.open("droping.txt");
  double data[27][2];
  for(int i=0;!file.eof();i++)
  {
    file>>data[i][0]>>data[i][1];
    data[i][0]=data[i][0]*4.25/1000/miu0;
  }
  double length=27.0;
  
  //初始化位置速度
	double position[n][dimension];
	double position_best_identity[n][dimension];
	double v[n][dimension];
	double value[n];
	double position_best_globle[dimension];
	ofstream mycout("data.csv");
  //初始化速度与位置
	for(int i=0;i<n;i++)
	{
    for(int j=0;j<dimension;j++)
    {
    srand(clock());
		position[i][j]=(uppper_limit[j]-lower_limit[j])*rand()*1.0/RAND_MAX;

    if((position[i][2]<lower_limit[2])||(position[i][2]>uppper_limit[2]))
    {
      position[i][2]=lower_limit[2]+(uppper_limit[2]-lower_limit[2])*rand()*1.0/RAND_MAX;
    }

    position_best_identity[i][j]=position[i][j];
		v[i][j]=Maxv[j]*rand()*1.0/RAND_MAX;
		//cout<<rand()<<" "<<clock()*1.0/CLOCKS_PER_SEC<<" "<<f(position[i])<<endl;
    }

    if((position[i][2]<lower_limit[2])||(position[i][2]>uppper_limit[2]))
    {
      position[i][2]=lower_limit[2]+(uppper_limit[2]-lower_limit[2])*rand()*1.0/RAND_MAX;
    }
    value[i]=fit(length,data,2,position[i]);
    for(int j=0;j<dimension;j++)
    {
		position_best_globle[j]=position[Min(value)][j];
		mycout<<position[i][j]<<",";
    }
    mycout<<value[i]<<endl;

	}

	//迭代
	for(int i=0;i<Maxiter;i++)
	{
		double omega;
		omega=(Initial_omega-Final_omega)*((Maxiter-i)/Maxiter)+Final_omega;
		for(int j=0;j<n;j++)
		{
			double _position[dimension],_value,_position_best_identity[dimension];//用于对比数值大小
      _value=value[j];
      velocity(v[j],position_best_identity[j],position[j],position_best_globle,omega);

      for(int k=0;k<dimension;k++)
      {
        _position[k]=position[j][k];
        _position_best_identity[k]=position_best_identity[j][k];
        mycout<<_position[k]<<",";
        v[j][k]=_velocity[k];
        //cout<<v[j][k]<<endl;
        position[j][k]=position[j][k]+v[j][k];
        //防止点跑到外面
        //srand(time(0));
        if((position[j][k]<lower_limit[k])||(position[j][k]>uppper_limit[k]))
        {
          position[j][k]=lower_limit[k]+(uppper_limit[k]-lower_limit[k])*rand()*1.0/RAND_MAX;
        }
        if((position_best_identity[j][k]<lower_limit[k])||(position_best_identity[j][k]>uppper_limit[k]))
        {
          position_best_identity[j][k]=lower_limit[k]+(uppper_limit[k]-lower_limit[k])*rand()*1.0/RAND_MAX;
        }

      }
      mycout<<_value<<endl;
      if((position[j][2]<lower_limit[2])||(position[j][2]>uppper_limit[2]))
      {
        position[j][2]=lower_limit[2]+(uppper_limit[2]-lower_limit[2])*rand()*1.0/RAND_MAX;
      }
      value[j]=fit(length,data,2,position[j]);
      //cout<<f(position_best_identity[j])-value[j]<<endl;

      //得到单个点最优位置
      if((fit(length,data,2,position_best_identity[j])-value[j])<=0){}
      else
      {
        for(int k=0;k<dimension;k++)
        {
          position_best_identity[j][k]=position[j][k];
        }
      }
		}



		//得到全局最优位置
		double _position_best_globle[dimension];
    for(int j=0;j<dimension;j++)
    {
      _position_best_globle[j]=position[Min(value)][j];
    }
		if(fit(length,data,2,position_best_globle)-fit(length,data,2,_position_best_globle)<=0){}
		else
		{
      for(int j=0;j<dimension;j++)
      {
        position_best_globle[j]=_position_best_globle[j];
      }
		}
    for(int j=0;j<dimension;j++)
    {
      if((position_best_globle[j]<lower_limit[j])||(position_best_globle[j]>uppper_limit[j]))
      {
        position_best_globle[j]=lower_limit[j]+(uppper_limit[j]-lower_limit[j])*rand()*1.0/RAND_MAX;
      }
    }
		
	}
  
  for(int j=0;j<dimension;j++)
  {
    cout<<position_best_globle[j]<<" ";
  }
	cout<<fit(length,data,2,position_best_globle)<<endl;
  
	mycout.close();

  for(int j=0;j<length;j++)
  {
    cout<<data[j][0]*miu0*1000<<" "<<Q(data[j][0],2,position_best_globle)<<endl;
  }
  cout<<thickness_a<<" "<<lambda_L<<" "<<Hc<<" "<<xi_0<<" "<<l_1<<" "<<l_2<<endl;
	return 0;
}