
#include"thermal_feed_back.cpp"

using namespace std;

double miu0=4*M_PI*1e-7; //Vacuum permeability 

double g=270;
double beta=1.0e-6;

const double C1=2.0,C2=2.0;
const int n=20;
const int Maxiter=50;
const int dimension=8;
const double Initial_omega=1.4;
const double Final_omega=0.4;
const double Maxv[dimension]={1.0e-7,1.0e-7,0.01/miu0,1.0e-7,1.0e-9,10.0e-9,0.5,0.5};
const double uppper_limit[dimension]={150.0e-9,33.0e-9,0.18/miu0,41.0e-9,10.0e-9,50.0e-9,1.0,1.0};
const double lower_limit[dimension]={1.0e-9,31.0e-9,0.181/miu0,40.0e-9,1.0e-9,10.0e-9,-1.0,-1.0};
/*

position[i]
i=0  thickness_a
i=1  lambda_L
i=2  Hc
i=3  xi_0
i=4  l_1
i=5  l_2
i=6  x1
i=7  x2

*/

double _velocity[dimension];


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
class FIT: public thermal_feed_back
{
  public:
    double fit(double length,double (*H_and_Q)[2])
    {
      double a=0;
      for(int i=0;i<length;i++)
      {
        set_H(H_and_Q[i][0]);
        caculate_T_distribution();
        //cout<<H_and_Q[i][1]<<" "<<Q()<<endl;
        a=a+pow((H_and_Q[i][1]-Q())/1.0e10,2);
      }
      /*
      cout<<thickness_a<<" "<<lambda_L<<" "<<
      Hc*miu0<<" "<<xi_0<<" "<<l_1<<" "<<" "<<l_2<<
      " "<<x1<<" "<<x2<<endl;
      */
      return a;
    };
    void set_parameters(double *position)
    {
      thickness_a=position[0];
      lambda_L=position[1];
      Hc=position[2];
      xi_0=position[3];
      l_1=position[4];
      l_2=position[5];
      x1=position[6];
      x2=position[7];
    };
};



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
	//ofstream mycout("data.csv");
  //初始化thermal feed back
  FIT fitparameter;
  fitparameter.set_l_phonon(1e-4);
  fitparameter.set_d(3e-3);
  fitparameter.set_N(50);
  fitparameter.set_RRR(300.0);
  fitparameter.set_T_He(2.0);
  fitparameter.init_T();
  fitparameter.set_H(0.1/miu0);
  fitparameter.caculate_T_distribution();

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

    fitparameter.set_parameters(position[i]);
    value[i]=fitparameter.fit(length,data);

    for(int j=0;j<dimension;j++)
    {
		position_best_globle[j]=position[Min(value)][j];
		//mycout<<position[i][j]<<",";
    }
    //mycout<<value[i]<<endl;

	}

	//迭代
	for(int i=0;i<Maxiter;i++)
	{
    cout<<i<<endl;
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
        //mycout<<_position[k]<<",";
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
      //mycout<<_value<<endl;
      if((position[j][2]<lower_limit[2])||(position[j][2]>uppper_limit[2]))
      {
        position[j][2]=lower_limit[2]+(uppper_limit[2]-lower_limit[2])*rand()*1.0/RAND_MAX;
      }

      fitparameter.set_parameters(position[j]);
      value[j]=fitparameter.fit(length,data);

      //cout<<f(position_best_identity[j])-value[j]<<endl;

      //得到单个点最优位置
      fitparameter.set_parameters(position_best_identity[j]);
      if((fitparameter.fit(length,data)-value[j])<=0){}
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
    double position_best_globle_fit,_position_best_globle_fit;

    fitparameter.set_parameters(position_best_globle);
    position_best_globle_fit=fitparameter.fit(length,data);
    fitparameter.set_parameters(_position_best_globle);
    _position_best_globle_fit=fitparameter.fit(length,data);

		if(position_best_globle_fit-_position_best_globle_fit<=0){}
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
  fitparameter.set_parameters(position_best_globle);
	cout<<fitparameter.fit(length,data)<<endl;
  /*
	//mycout.close();
  double aaa[8]={7.29186e-08,3.35301e-08,150703,3.95602e-08,7.24981e-09,5.093e-08,-1.54024,1.84796};
  fitparameter.set_parameters(aaa);
  */
  for(int j=0;j<length;j++)
  {
    fitparameter.set_H(data[j][0]);
    fitparameter.caculate_T_distribution();
    cout<<data[j][0]*miu0*1000<<" "<<fitparameter.Q()<<" "<<fitparameter.Delta_T_H(fitparameter.T[0],data[j][0])/(9.23*1.38e-23)<<endl;
  }
  cout<<fitparameter.thickness_a<<" "<<fitparameter.lambda_L<<" "<<fitparameter.Hc<<" "<<
  fitparameter.xi_0<<" "<<fitparameter.l_1<<" "<<fitparameter.l_2<<" "<<fitparameter.x1<<" "<<fitparameter.x2<<" ";
  
  cout<<fitparameter.fit(length,data)<<" "<<fitparameter.Hc*miu0<<" ";
  double kappa1=fitparameter.lambda_L_T_H(fitparameter.T[0],0)/fitparameter.xi_0;
  double kappa2=fitparameter.lambda_L_T_H(fitparameter.T[0],0)/fitparameter.xi_0;
  cout<<kappa1<<" "<<kappa2<<" "<<
  pow(2/kappa1,0.5)*fitparameter.Hc*miu0<<" "<<pow(2/kappa2,0.5)*fitparameter.Hc*miu0<<endl;
	return 0;
}