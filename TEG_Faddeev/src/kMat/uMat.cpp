
#include <complex>
#include "uMat.h"

/*caculate (208)-(214), ref. Glockle_3Nreview_96*/
//All of the comment in this file are using for testing code. If not necessary, don't delete them.

REAL UMAT::hat(REAL x)   //hat(x) is a function to caculate \hat{x}
{
  REAL a=2*x+1;
  return a;
}

REAL UMAT::delta(REAL i,REAL j)
{

  if (abs(i-j) < 1.0e-6)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

COMPLEX UMAT::U_sigma(REAL lambda_,REAL Sigma_, REAL lambda, REAL Sigma, string path) //Huber_1995 2.26
{
  using namespace std;
  REAL ratio_NQovrNP {1.0};
  Fstr inputTable; // Input parameters
  //Set default parameters
  inputTable.maxL={2};
  inputTable.maxOpT={0.5};
  inputTable.opJ={0.5};
  inputTable.parity={0};
  inputTable.channels={0,1,12,2,13,5};
  inputTable.mambda={400};
  inputTable.NP={3,8,11};
  inputTable.potential={"chengdu_MMWLY_cmplx_400"};
  inputTable.NX={7};
  inputTable.CMRatio={0.6};
  inputTable.RFMTHD={"Scant"};
  inputTable.BE01={-15.0,-3.0};
  inputTable.error={1.0E-2};
  inputTable.phi={0.5};

  // bool condition=inputTable.cookInput(path);
  inputTable.cookInput(path);

  // Generalize channel table
  Chs chsTable {};

  // Parameters needed for generating channels table
  chsTable.writeMaxL(inputTable.maxL[0]);
  chsTable.writeMaxOpT(inputTable.maxOpT[0]);
  chsTable.writeOpJ(inputTable.opJ[0]);
  chsTable.writeParity(inputTable.parity[0]);
  chsTable.genChs();
  /*
  chsTable.showChannels(); // Show channel table for user to select
  */
  // Initialize the selected channels
  size_t nChannels = inputTable.channels.size();
  int labelChs[nChannels];
  for (size_t i = 0; i < nChannels; i++) {
    labelChs[i] = inputTable.channels[i];
  }
  REAL **channels = new REAL*[nChannels];
  for(size_t i = 0; i < nChannels; i++) {
    channels[i] = new REAL[_N_QMN_CHN_];
  }
  chsTable.outChannels(channels, nChannels, labelChs);

  // Print out the channels the user have seleted
  cout << "/////////////////////////////////" << endl;
  if(inputTable.channels_j)
  {
    printf("You have selected channels:\n");
  }
  else
  {
    printf("You have selected channels: (default)\n");
  }

  printf("%-8s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s%-5s\n",
    "alpha", "l", "s", "j", "lmbd", "I", "T",  "lp", "opT", "Pi", "opJ");
  for(size_t i = 0; i < nChannels; i++) {
    printf(" %4d   ", labelChs[i]);
    for(int j = 0; j < _N_QMN_CHN_; j++)
    printf("%-5.1f", channels[i][j]);
    printf("\n");
  }
  cout << "/////////////////////////////////" << endl;

  // Initialize k matrix, note that all configurations can be changed later.
  REAL mambda = inputTable.mambda[0]; // Cutoff
  REAL _CTanMesh = mambda * inputTable.CMRatio[0];

  // debug
  cout << "mambda = " << mambda << endl;
  cout << "CTanMesh = " << _CTanMesh << endl;

  // Define the mesh points' range
  /*
  REAL dPara = inputTable.NP[0],
    para0 = inputTable.NP[1],
    para1 = inputTable.NP[2];
  */
  // size_t lenPara = (para1 - para0) / dPara + 1;
  /*
  cout << "/////////////////////////////////" << endl;

  //print parameters this cacualating will be used
  if(condition)
  {
    cout << "Input file = " << path << endl;
  }
  else
  {
    cout << "Input file = Nothing" << endl;
    cout << "Those parameters are default parameters." << endl;
  }
  inputTable.PrintParameter();
  cout << "/////////////////////////////////" << endl;
    // current date/time based on current system
  */

  // time_t now = time(0);

  /*
  cout << ">>>>>> Job starts at " << ctime(&now) << endl;
  printf( "Mesh Points\n");
  */
  fflush(stdout);
  int NP = inputTable.NP[1]; // numer of p mesh

  // Initialize chengdu_nnscat global functions/variables
  init_nnscat_chengdu_C();
  // Initialize TEG_Faddeev global functions/variables
  facLabel();


  /* Complex-valued pots contains "cmplx" */
  if ( inputTable.potential[0].find("cmplx") != string::npos )
  {

    kMatrix<COMPLEX> kMatC {nChannels, _CTanMesh, channels, inputTable.opJ[0]};
    REAL phi = inputTable.phi[0]/180.0*3.141593;
    const std::complex<double> ii(0.0, 1.0);
    kMatC.setFacContour(exp(-ii*phi));
    kMatC.setPotName(inputTable.potential[0]); // potential
    kMatC.setMass(_AvgMassN_);
    kMatC.setNX(inputTable.NX[0]);
    // printf("NP = %-20d\n", NP);
    fflush(stdout);
    kMatC.setNP(NP);

    // debug
    kMatC.setNAuxMesh((int) NP*1.4);
    // kMatC.setNAuxMesh(100);

    kMatC.setNQ(int(NP*ratio_NQovrNP));

    // debug
    // cout << "uMat.cpp:: 1" << endl;
    // time_t now;
    // now = time(0);
    // cout << _TIMER_COL_ << ctime(&now) << flush;

    kMatC.genGTable();
    kMatC.allocKMat();
    SetMass(kMatC.getMass());
    SetopJ(kMatC.getopJ());
    REAL J,J_;
    COMPLEX Usigma(0.0,0.0);
    for(J_=fabs(lambda_-1.0/2);J_<=lambda_+1.0/2;J_++)
    {
      for(J=fabs(lambda-1.0/2);J<=lambda+1.0/2;J++)
      {
        //All of these comment are test code. If not necessary, don't delete them.
        //printf("lambda_=%f 1/2 J_=%f jd=%f opJ=%f Sigma_=%f\n",lambda_,J_,jd,opJ,Sigma_);
        //printf("lambda=%f 1/2 J=%f jd=%f opJ=%f Sigma=%f\n",lambda,J,jd,opJ,Sigma);
        // COMPLEX Ukf = kMatC.get_U(energyLoc, q0, lambda_, J_, q0+1.0e-9, lambda, J);
        // COMPLEX Ukf = kMatC.get_U(energyLoc, q0, lambda_, J_, q0, lambda, J);
        COMPLEX Ukf = kMatC.get_U_elastic(q0, lambda_, J_, lambda, J);
        // cout << "kMatC.get_U_elastic = " << " " << Ukf << endl;

        // debug
        cout << Ukf*2.0*M_PI*_AvgMassN_*_HBARC_ << endl;

        cout << "* NP = " << NP << endl;

        Usigma=Usigma+sqrt(hat(J_)*hat(Sigma_))*pow(-1,opJ-J_)*sixJ(lambda_,1.0/2,J_,jd,opJ,Sigma_)
          *sqrt(hat(J)*hat(Sigma))*pow(-1,opJ-J)*sixJ(lambda,1.0/2,J,jd,opJ,Sigma)
          *Ukf;
        // Usigma=Usigma+sqrt(hat(J_)*hat(Sigma_))*pow(-1,opJ-J_)*sixJ(lambda_,1.0/2,J_,jd,opJ,Sigma_)
        //   *sqrt(hat(J)*hat(Sigma))*pow(-1,opJ-J)*sixJ(lambda,1.0/2,J,jd,opJ,Sigma)
        //   *kMatC.get_U(energyLoc,q0,lambda_,J_,q0+1.0e-9,lambda,J);
        //cout<<sqrt(hat(J_)*hat(Sigma_))<<endl;
        //cout<<pow(-1,opJ-J_)<<endl;
        //cout<<sixJ(lambda_,1.0/2,J_,jd,opJ,Sigma_)<<endl;
        //cout<<sqrt(hat(J)*hat(Sigma))<<endl;
        //cout<<pow(-1,opJ-J)<<endl;
        //cout<<sixJ(lambda,1.0/2,J,jd,opJ,Sigma)<<endl;
        //cout<<kMatC.get_U(energyLoc,q0,lambda_,J_,q0+pow(10.0,-9),lambda,J)<<endl;
        //cout<<Usigma<<endl;
      }
    }


    // now = time(0);
    //cout << endl << "------ Currently  at " << ctime(&now)  << endl;
    fflush(stdout);
    return Usigma;

  }
  // now = time(0);
  //cout << "<<<<<< Job ends at   " << ctime(&now)  << endl;
  for(size_t i = 0; i < nChannels; i++)
  {
    delete[] channels[i];
  }
  delete[] channels;
  return 0;

}



COMPLEX UMAT::S(REAL lambda_,REAL Sigma_, REAL lambda, REAL Sigma,string path)// Huber_1995 2.23
{
  using namespace std;
  COMPLEX a,i={0,1};
  a=delta(lambda,lambda_)*delta(Sigma_,Sigma)-(4.0/3)*M_PI*i*q0*m*
  pow(i,(lambda_-lambda))*U_sigma(lambda_,Sigma_,lambda,Sigma,path);
  return a;
}

void UMAT::Build_S_Matrix(REAL parity,string path) //Parity here is different from other files. All of the parity must be 1 or -1
{                                                   //Huber_1995 2.37
  std::cout<<"S_Matrix"<<endl;
  if(opJ-0.5<0.00000001)
  {
    std::cout<<"1/2"<<endl;
    std::cout<<opJ<<endl;
    REAL i1,i2,j1,j2;
    for(int i=0;i<=1;i++)
    {
      switch(i)
      {
        case 0:
          i1=parity*1/2.0;
          i2=1/2.0;
          break;
        case 1:
          i1=parity*1/2.0;
          i2=3/2.0;
          break;
      }
      for(int j=0;j<=1;j++)
      {
        switch(j)
        {
          case 0:
            j1=parity*1/2.0;
            j2=1/2.0;
            break;
          case 1:
            j1=parity*1/2.0;
            j2=3/2.0;
            break;
        }
        S_Matrix.push_back(S(opJ+i1,i2,opJ+j1,j2,path));
        std::cout<<S_Matrix[2*i+j]<<" ";
      }
      std::cout<<endl;
    }
  }
  else if(opJ-1.5<0.00000001)
  {
    std::cout<<"3/2"<<endl;
    std::cout<<opJ<<endl;
    REAL i1,i2,j1,j2;
    for(int i=0;i<=2;i++)
    {
      switch(i)
      {
        case 0:
          i1=-3/2.0*parity;
          i2=3/2.0;
          break;
        case 1:
          i1=1/2.0*parity;
          i2=1/2.0;
          break;
        case 2:
          i1=1/2.0*parity;
          i2=3/2.0;
          break;
      }
      for(int j=0;j<=2;j++)
      {
        switch(j)
        {
          case 0:
            j1=-3/2.0*parity;
            j2=3/2.0;
            break;
          case 1:
            j1=1/2.0*parity;
            j2=1/2.0;
            break;
          case 2:
            j1=1/2.0*parity;
            j2=3/2.0;
            break;
        }
        S_Matrix.push_back(S(opJ+i1,i2,opJ+j1,j2,path));
        std::cout<<S_Matrix[3*i+j]<<" ";
      }
      std::cout<<endl;
    }
  }
  else
  {
    std::cout<<"error"<<endl;
  }

}

void UMAT::Build_delta_Matrix_and_U_Matrix(REAL parity,string path) //Huber_1995 2.38 Delta matrix and U matrix.
{
  Build_S_Matrix(parity,path);
  COMPLEX ii={0,1};
  int N;
  N=int(sqrt(S_Matrix.size()));
  MKL_Complex16 SMatrix[N*N];
  for(long unsigned int i=0;i<S_Matrix.size();i++)
  {
    SMatrix[i]=S_Matrix[i];
  }
  int LDA=N,LDVL=N,LDVR=N;
  COMPLEX w[N],vl[LDVL*N], vr[LDVR*N];
  LAPACKE_zgeev(LAPACK_COL_MAJOR,'V', 'V', N, SMatrix, LDA, w, vl, LDVL, vr, LDVR);
  std::cout<<"delta_Matrix"<<endl;
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      if(i==j)
      {
        delta_Matrix.push_back(-ii*(log(pow(w[i].imag(),2)+pow(w[i].real(),2))/4)+0.5*atan(w[i].imag()/w[i].real()));
        std::cout<<1<<delta_Matrix[N*i+j]<<" ";
      }
      else
      {
        delta_Matrix.push_back(0);
        std::cout<<1<<delta_Matrix[N*i+j]<<" ";
      }
      std::cout<<endl;
    }
  }
  std::cout<<"U_Matrix"<<endl;
  for(int i=0;i<N;i++)
  {
    for(int j=0;j<N;j++)
    {
      U_Matrix.push_back(vl[N*i+j]);
      std::cout<<2<<U_Matrix[N*i+j]<<" ";
    }
    std::cout<<endl;
  }
}
/*
void UMAT::Caculate_Mixing_Parameters(REAL parity,string path)//Huber_1995 2.40 2.41 2.42
{
  Build_delta_Matrix_and_U_Matrix(parity,path);
  //epslion,xi,eta;
  if(U_Matrix.size()==9)
  {
    Mixing_Parameters.push_back(-std::atan2(U_Matrix[7],U_Matrix[8]));
    Mixing_Parameters.push_back(atan2(-U_Matrix[6],std::sqrt(std::pow(U_Matrix[7],2)+std::pow(U_Matrix[7],2))));
    Mixing_Parameters.push_back(-std::atan2(U_Matrix[0],U_Matrix[3]));
  }
  if(U_Matrix.size()==4)
  {
    Mixing_Parameters.push_back(std::acos(U_Matrix[0]));
  }
}
*/
