

#include<iostream>
#include<cmath>
#include<ctime>
#include<assert.h>
#include<fstream>

using namespace std;

const double C1=2.0,C2=2.0;
const int n=100;
const int Maxiter=400;
const double Maxv[6]={1000.0,1000.0,1000.0,1000.0,1000.0,1000.0};
const double Initial_omega=1.4;
const double Final_omega=0.4;
const int dimension=6;
const double uppper_limit[6]={8190.0,8190.0,8190.0,8190.0,8190.0,8190.0};
const double lower_limit[6]={1.0,1.0,1.0,1.0,1.0,1.0};

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

int main()
{
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
		position_best_identity[i][j]=position[i][j];
		v[i][j]=Maxv[j]*rand()*1.0/RAND_MAX;
		//cout<<rand()<<" "<<clock()*1.0/CLOCKS_PER_SEC<<" "<<f(position[i])<<endl;
    }
    value[i]=f(position[i]);
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
          position[j][k]=(uppper_limit[k]-lower_limit[k])*rand()*1.0/RAND_MAX;
        }
        if((position_best_identity[j][k]<lower_limit[k])||(position_best_identity[j][k]>uppper_limit[k]))
        {
          position_best_identity[j][k]=(uppper_limit[k]-lower_limit[k])*rand()*1.0/RAND_MAX;
        }

      }
      mycout<<_value<<endl;
      value[j]=f(position[j]);
      //cout<<f(position_best_identity[j])-value[j]<<endl;

      //得到单个点最优位置
      if((f(position_best_identity[j])-value[j])<=0){}
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
		if(f(position_best_globle)-f(_position_best_globle)<=0){}
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
        position_best_globle[j]=(uppper_limit[j]-lower_limit[j])*rand()*1.0/RAND_MAX;
      }
    }
		
	}
  for(int j=0;j<dimension;j++)
  {
    cout<<position_best_globle[j]<<" ";
  }
	cout<<f(position_best_globle)<<endl;
	mycout.close();
	return 0;
}

