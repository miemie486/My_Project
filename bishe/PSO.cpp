
#include<iostream>
#include<cmath>
#include<ctime>
#include<assert.h>
#include<fstream>

using namespace std;

const double C1=2.0,C2=2.0;
const int n=100;
const int Maxiter=200;
const double Maxv=1000.0;
const double Initial_omega=1.4;
const double Final_omega=0.4;


double f(double x)
{
	double a;
	a=x*sin(x/730.0)+x*cos(x/260.0);
	return a;
}

double velocity(double v_before,double identity_best,double identity_now,double globle_best,double omega)
{
	srand(clock());
	double v;
	v=omega*v_before+C1*rand()/RAND_MAX*(identity_best-identity_now)+C2*rand()/RAND_MAX*(globle_best-identity_now);
	return v/50.0;
}

int Max(double *array)
{
	double min=array[0];
	int a=0;
	for(int i=0;i<n;i++)
	{
		if((min-array[i])>=0){}
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
	double position[n];
	double position_best_identity[n];
	double v[n];
	double value[n];
	double position_best_globle;
	ofstream mycout("data.csv");
	for(int i=0;i<n;i++)
	{
		//srand(time(0));
		position[i]=8190.0*rand()/RAND_MAX;
		position_best_identity[i]=position[i];
		v[i]=Maxv*rand()/RAND_MAX;
		value[i]=f(position[i]);
		position_best_globle=position[Max(value)];
		mycout<<position[i]<<","<<value[i]<<endl;
		//cout<<8190.0*rand()/RAND_MAX<<""<<value[i]<<endl;
	}
	//迭代
	for(int i=0;i<Maxiter;i++)
	{
		double omega;
		omega=(Initial_omega-Final_omega)*((Maxiter-i)/Maxiter)+Final_omega;
		for(int j=0;j<n;j++)
		{
			double _position,_value,_position_best_identity;//用于对比数值大小
			_position=position[j];
			_value=value[j];
			_position_best_identity=position_best_identity[j];
			mycout<<_position<<","<<_value<<endl;
			v[j]=velocity(v[j],position_best_identity[j],position[j],position_best_globle,omega);

			position[j]=position[j]+v[j];
			//防止点跑到外面
			if((position[j]<1.0)||(position[j]>8191.0))
			{
				position[j]=8190.0*rand()/RAND_MAX;
			}
			if((position_best_identity[j]<1.0)||(position_best_identity[j]>8191.0))
			{
				position_best_identity[j]=8190.0*rand()/RAND_MAX;
			}
			value[j]=f(position[j]);
			//得到单个点最优位置
			if((f(position_best_identity[j])-value[j])>=0){}
			else
			{
				position_best_identity[j]=position[j];
			}
		}
		//得到全局最优位置
		double _position_best_globle;
		_position_best_globle=position[Max(value)];
		if(f(position_best_globle)-f(_position_best_globle)>=0){}
		else
		{
			position_best_globle=_position_best_globle;
		}
		if((position_best_globle<1.0)||(position_best_globle>8191.0))
		{
			position_best_globle=8190.0*rand()/RAND_MAX;
		}
	}
	cout<<position_best_globle<<" "<<f(position_best_globle)<<endl;
	mycout.close();
	return 0;
}

