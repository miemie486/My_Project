
#include <iostream>

#include "kMatGen.h"
#include "uMat.h"
#include "msg.h"

COMPLEX scatlen(REAL q0, REAL n, string path)
{
  UMAT umatrix;
  REAL sigma=(n-1)/2;
  COMPLEX a_;
  umatrix.Setq0(q0);  //set q0
  REAL H2BE = -2.22107; // MMWLY_400
  // REAL H2BE = -2.22884; // MMWLY_600
  umatrix.SetEnergyLoc(H2BE + 0.75*pow(q0,2)/_AvgMassN_);
  a_=(2.0*M_PI*umatrix.GetMass()/3.0)*umatrix.U_sigma(0,sigma,0,sigma,path);
  //std::cout<<umatrix.U_sigma(0,sigma,0,sigma,path)<<std::endl;

  /*
  int N=int(sqrt(umatrix.S_Matrix.size()));
  std::cout<<"S_Matrix"<<endl;
  for(int i=0;i<N;i++)
  {
    for(int j;j<N;j++)
    {
      std::cout<<umatrix.S_Matrix[N*i+j]<<" ";
    }
    std::cout<<endl;
  }
  std::cout<<"U_Matrix"<<endl;
  for(int i=0;i<N;i++)
  {
    for(int j;j<N;j++)
    {
      std::cout<<umatrix.U_Matrix[N*i+j]<<" ";
    }
    std::cout<<endl;
  }
  std::cout<<"delta_Matrix"<<endl;
  for(int i=0;i<N;i++)
  {
    for(int j;j<N;j++)
    {
      std::cout<<umatrix.delta_Matrix[N*i+j]<<" ";
    }
    std::cout<<endl;
  }
  */
  return a_;
  //return umatrix.U_sigma(0,sigma,0,sigma,path);
}

int main(int argc, char* argv[])
{
  string fname, default_fname = "scat_400.in";
  REAL q0, n;
  COMPLEX rslt;

  set_verbose_level_fdv(_VBSLVL_HIGH_);
  // set_verbose_level_fdv(_VBSLVL_LOW_);
  show_message_timestamp_fdv("\n\nscatlen: starts", _VBSLVL_LOW_);

  if (argc < 2)
  {
    fname = default_fname;
  }
  else
  {
    if (argv[1][0] == '\0')
    {
      fname = default_fname;
    }
    else
    {
      fname = argv[1];
    }
  }
  if(argc < 4)
  {
    q0=std::stod(argv[1]);
    n=std::stod(argv[2]);
  }
  else
  {
    q0=std::stod(argv[2]);
    n=std::stod(argv[3]);
  }

  rslt = scatlen(q0,n,fname) * _HBARC_;
  std::cout<<"////////////////////////"<<std::endl;
  std::cout << "* q0 = " << q0 << endl;
  // std::cout<<"* [Final] a_" << n << " = " << rslt << std::endl;
  printf("* [Final] a_%d = (%-22.15e,  %-22.15e)\n", (int) n, rslt.real(), rslt.imag());

  show_message_timestamp_fdv("scatlen: ends", _VBSLVL_LOW_);
  return 0;
}
