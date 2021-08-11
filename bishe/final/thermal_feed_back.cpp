
#include"R_BCS.cpp"

using namespace std;

class thermal_feed_back:public two_layer_model
{
  private:
    double A=0.141;
    double a=7.52e-7;
    double B=4.34e3;
    double D=1/(2.34e2);
    double g=270;

    double l_phonon;
    double RRR;
    double d; //Material thickness
    double T_He;
    double H; //The intensity of the surface magenatic field, dimension A/m
    int N ;   //The number of the layer the material will be divided.
    
  public:
    void set_l_phonon(double l){l_phonon=l;};
    void set_RRR(double _RRR){RRR=_RRR;};
    void set_d(double _d){d=_d;};
    void set_T_He(double T){T_He=T;};
    void set_N(double _N){N=_N;};
    void set_H(double _H){H=_H;};
    void init_T(){T.resize(N);};
    double y(double T);
    double f_0(double p);
    double f(double _y);
    double R(double _y);
    double Kappa(double T, double RRR, double l);
    double H_Kapitza(double T1, double T2);
    void caculate_T_distribution_step1(double R_BCS, double T_0);
    void caculate_T_distribution();
    double Q(){return g/R_BCS(T[0],H);};
    vector<double> T;

};

double thermal_feed_back::y(double T)
{
  double aa=Delta_0/(kB*T)*pow(cos(M_PI/2*pow(T/Tc,2)),0.5);
  return aa;
}

double thermal_feed_back::f_0(double p)
{
  double aa;
  aa=-(1+p/2+p*p/3+p*p*p/4);
  return aa;
}

double thermal_feed_back::f(double _y)
{

  double t=-exp(-_y);
  double intg=0;
  double n=10;
  double h=t/n;
  intg=(h/3.0)*(f_0(__DBL_MIN__)+f_0(t));
  for(int i=1;i<=n/2;i++)
  {
    intg=intg+(4.0*h/3.0)*f_0((2.0*i-1.0)*h)+(2.0*h/3.0)*f_0(2.0*i*h);
  }
  intg=intg-2.0*h/3.0*f_0(t);

  return intg;
}

double thermal_feed_back::R(double _y)
{
  double aa;
  aa=(12/(M_PI*M_PI))*(f(_y)+_y*log(1+exp(-_y))+_y*_y/(2*(1+exp(_y))));
  return aa;
}

double thermal_feed_back::Kappa(double T, double RRR, double l)  // l is the phonon mean free path
{
  double aa;
  double _y=y(T);
  aa=R(_y)/(1/(A*RRR*T)+a*T*T)+1/(1/(D*T*T*exp(_y))+1/(B*l*T*T*T));
  return aa;
}

double thermal_feed_back::H_Kapitza(double T1, double T2)
{
  double t=(T1-T2)/T2;
  double aa;
  double f=1+3/2*t+t*t+1/4*t*t*t;
  aa=200*pow(T2,4.65)*f;
  return aa;
}

void thermal_feed_back::caculate_T_distribution_step1(double R_BCS, double T_0)
{
  //Initialize T
  T[0]=T_0;
  double q=0.5*R_BCS*H*H;
  for(int i=1;i<N;i++)
  {
    if(i<N-1)
    {
      T[i]=-q*d/(N*Kappa(T[i-1],RRR,l_phonon))+T[i-1];
    }
    else
    {
      double xn=1,xn_min=0,xn_plus=0;
      while(fabs(xn-xn_min)>1e-6)
      {
        xn_plus=xn-(xn-T_He-q/H_Kapitza(xn,T_He))*(xn-xn_min)/
        (xn-q/H_Kapitza(xn,T_He)-xn_min+q/H_Kapitza(xn_min,T_He));
        xn_min=xn;
        xn=xn_plus;
      }
      T[i]=xn;
    }
  }
}

void thermal_feed_back::caculate_T_distribution()
{
  
  double T0_n=T_He+2.0,T0_n_min=T_He+1.1,T0_n_plus;
  while(fabs(T0_n-T0_n_min)>1e-6)
  {
    double q=0.5*R_BCS(T0_n,H)*H*H;
    double fn,fn_min;
    caculate_T_distribution_step1(R_BCS(T0_n,H),T0_n);
    fn=T[N-2]-T[N-1]-q*d/(N*Kappa(T[N-2],RRR,l_phonon));
    caculate_T_distribution_step1(R_BCS(T0_n,H),T0_n_min);
    fn_min=T[N-2]-T[N-1]-q*d/(N*Kappa(T[N-2],RRR,l_phonon));
    T0_n_plus=T0_n-fn*(T0_n-T0_n_min)/(fn-fn_min);
    T0_n_min=T0_n;
    T0_n=T0_n_plus;
  }
}

/*
int main()
{

  double miu0=4*M_PI*1e-7;
  double N=50;
  double H=0.15/miu0;

  thermal_feed_back TFB;
  TFB.set_l_phonon(1e-4);
  TFB.set_d(3e-3);
  TFB.set_N(N);
  TFB.set_RRR(300.0);
  TFB.set_T_He(2.0);
  TFB.init_T();
  TFB.set_H(H);
  cout<<TFB.alpha(TFB.T[0],H,TFB.l_1)<<endl;
  TFB.caculate_T_distribution();
  
  for(int i=0;i<N;i++)
  {
    cout<<TFB.T[i]<<endl;
  }
  cout<<TFB.alpha(TFB.T[0],H,TFB.l_1)<<endl;
  cout<<270/TFB.R_BCS(TFB.T[0],H)<<endl;

  
}
*/