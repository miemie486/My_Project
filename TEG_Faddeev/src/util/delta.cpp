#include "delta.h"
#include "def.h"

using namespace std;

REAL Upara_list[_MAXN_UparaL_] {0};//store all possible of parameters for U_lambda_I
bool Ulist_initialized {false};


//list the parameters(lambda_p, I_p, lambda, I) of all needed U_lambda_I, 
//store them orderly in  Upara_list[], and output them on the screen.
void get_Upara_list(REAL opJ_, REAL parity_){
  REAL parity = abs(parity_)<_momentum_EPS_?1.0:-1.0;
  std::cout << "the parameters for U_lambda_I \n" ;
  printf("%-10s%-7s%-7s%-7s%-7s\n", "idx", "lmbd_p", "I_p", "lmbd", "I");

  int idx = 0;
  for(REAL lambda_p = min(abs(opJ_-1.5),abs(opJ_-0.5)); lambda_p <= opJ_+1.5; lambda_p+=1.0){
    if(abs(pow(-1.0,lambda_p)-parity) > _momentum_EPS_ )
      continue;
    REAL minI_p = max(abs(lambda_p-0.5) , abs(1.0-opJ_));
    REAL maxI_p = min(abs(lambda_p+0.5), abs(1.0+opJ_));
    for(REAL I_p = minI_p; I_p <= maxI_p; I_p+=1.0){
      for(REAL lambda = min(abs(opJ_-1.5),abs(opJ_-0.5)); lambda <= opJ_+1.5; lambda+=1.0){
        if(abs(pow(-1.0,lambda)-parity) > _momentum_EPS_ )
          continue;
        REAL minI = max(abs(lambda-0.5) , abs(1.0-opJ_));
        REAL maxI = min(abs(lambda+0.5) , abs(1.0+opJ_));
        for(REAL I = minI; I <= maxI; I+=1.0){
          printf("%-10d%-7.1f%-7.1f%-7.1f%-7.1f\n", idx, lambda_p, I_p, lambda, I);
          Upara_list[idx*4] = lambda_p; Upara_list[idx*4+1] = I_p; Upara_list[idx*4+2] = lambda; Upara_list[idx*4+3] = I;
          idx++;
        }
      }
    }
  }
  Ulist_initialized = true;
}

//return the idx of U given parameters: lambda_p, I_p, lambda, I
//the idx of certain U is decided by Upara_list
int SerialofIdx_U_list(REAL opJ_, REAL lambda_p, REAL I_p, REAL lambda, REAL I){
  if(!Ulist_initialized){
    cout<<"Upara_list is not ready!"<<endl;
    exit(-1);
  }

  int N = (opJ_ == 0.5)? 4:9;
  for(int i = 0; i < N; i++){
    if(abs(lambda_p - Upara_list[i*4+0]) < _momentum_EPS_ && abs(I_p - Upara_list[i*4+1]) < _momentum_EPS_ &&
        abs(lambda - Upara_list[i*4+2]) < _momentum_EPS_ && abs(I - Upara_list[i*4+3]) < _momentum_EPS_)
      return i;
  }
  cout<<"can't find the correlated index of U"<<endl;
  exit(-1);
}

//equ.(40) in codenote
COMPLEX U_lambda_sigma(REAL q0, REAL opJ, REAL lambda_prime, REAL sigma_prime, REAL lambda, REAL sigma, COMPLEX * U_array){
  facLabel();
  COMPLEX result ={0.0};
  REAL minI_prime = max(abs(lambda_prime-0.5) , abs(1.0-opJ));
  REAL maxI_prime = min(abs(lambda_prime+0.5), abs(1.0+opJ));
  REAL minI = max(abs(lambda-0.5) , abs(1.0-opJ));
  REAL maxI = min(abs(lambda+0.5) , abs(1.0+opJ));
  for(REAL I_prime = minI_prime; I_prime <= maxI_prime; I_prime+=1.0){
    for(REAL I = minI; I <= maxI; I+=1.0){
      result=result+sqrt((2.0*I_prime+1.0)*(2.0*sigma_prime+1.0)*(2.0*I+1.0)*(2.0*sigma+1.0))*pow(-1.0,int(opJ-I_prime+opJ-I))*
                sixJ(lambda_prime,0.5,I_prime,1.0,opJ,sigma_prime)*sixJ(lambda,0.5,I,1.0,opJ,sigma)*
                U_array[SerialofIdx_U_list(opJ, lambda_prime, I_prime, lambda, I)];
    }
  }
  return result;
}

//compute related S_lambda_sigma given U_lambda_sigma .
COMPLEX S_lambda_sigma(REAL q0, REAL opJ, REAL lambda_prime, REAL sigma_prime, REAL lambda, REAL sigma, COMPLEX * U_array){
  REAL delta_= 0.0;
  if(abs(lambda_prime-lambda)<_momentum_EPS_&&abs(sigma_prime-sigma)<_momentum_EPS_)
    delta_ = 1.0;
  COMPLEX ii = {0.0,1.0};
  return delta_-4.0/3.0*M_PI*ii*q0*_AvgMassN_*pow(ii,int(lambda_prime-lambda))*
          U_lambda_sigma(q0,opJ,lambda_prime,sigma_prime,lambda,sigma, U_array);
}


//generate the matrix of S_lambda_sigma
void SMatGen(COMPLEX * Smat, int N_Smat, REAL opJ, REAL _parity_, REAL q0, COMPLEX * U_array){// N_Smat is the size of Smat
  REAL parity_ =abs(_parity_)<_momentum_EPS_?1.0:-1.0; //_parity_ is the parameter input from argumentlist(0 means even and 1 means old), parity_ is the real parity
  int N = (opJ == 0.5)?4:9;
  if(N != N_Smat){
    cout<<"the size of Smat must be in accord with opJ"<<endl;
    exit(-1);
  }
  if(N_Smat == 4){
    Smat[0] = S_lambda_sigma(q0, opJ, abs(opJ+parity_*1.5), 1.5, abs(opJ+parity_*1.5), 1.5, U_array);
    Smat[1] = S_lambda_sigma(q0, opJ, abs(opJ+parity_*1.5), 1.5, opJ-parity_*0.5, 0.5, U_array);
    Smat[2] = S_lambda_sigma(q0, opJ, opJ-parity_*0.5, 0.5, abs(opJ+parity_*1.5), 1.5, U_array);
    Smat[3] = S_lambda_sigma(q0, opJ, opJ-parity_*0.5, 0.5, opJ-parity_*0.5, 0.5, U_array);
  }
  else{
    int id = int(opJ+0.5);
    REAL sign = (id%2==0?1.0:-1.0);
      Smat[0] = S_lambda_sigma(q0, opJ, opJ-sign*parity_*1.5, 1.5, opJ-sign*parity_*1.5, 1.5, U_array);
      Smat[1] = S_lambda_sigma(q0, opJ, opJ-sign*parity_*1.5, 1.5, opJ+sign*parity_*0.5, 0.5, U_array);
      Smat[2] = S_lambda_sigma(q0, opJ, opJ-sign*parity_*1.5, 1.5, opJ+sign*parity_*0.5, 1.5, U_array);
      Smat[3] = S_lambda_sigma(q0, opJ, opJ+sign*parity_*0.5, 0.5, opJ-sign*parity_*1.5, 1.5, U_array);
      Smat[4] = S_lambda_sigma(q0, opJ, opJ+sign*parity_*0.5, 0.5, opJ+sign*parity_*0.5, 0.5, U_array);
      Smat[5] = S_lambda_sigma(q0, opJ, opJ+sign*parity_*0.5, 0.5, opJ+sign*parity_*0.5, 1.5, U_array);
      Smat[6] = S_lambda_sigma(q0, opJ, opJ+sign*parity_*0.5, 1.5, opJ-sign*parity_*1.5, 1.5, U_array);
      Smat[7] = S_lambda_sigma(q0, opJ, opJ+sign*parity_*0.5, 1.5, opJ+sign*parity_*0.5, 0.5, U_array);
      Smat[8] = S_lambda_sigma(q0, opJ, opJ+sign*parity_*0.5, 1.5, opJ+sign*parity_*0.5, 1.5, U_array);
  }
}

//compute phase shifts and mixangles
void Gen_deltaArray(COMPLEX * deltaArray_, int N_delta, COMPLEX * Mixangle, int N_mixangle, COMPLEX * Smat_, int N_Smat){
  //deltaArray_ is the matrix storing phase shift values, N_delta is the size of deltaArray_
  //N_Smat is the size of Smat_
  //Note that the mixangle[] must be COMPLEX type;
  int N = sqrt(N_Smat); //N is the dimension of Smat

  if(N_delta != N){
    std::cout<<"The size of deltaArray must be eaqul to the dimension of Smat!"<<endl;
    exit(-1);
  }
  if((N == 2 && N_mixangle != 1) || (N == 3 && N_mixangle != 3)){
    std::cout<<"The size of mixangle-array must be in accord with the dimension of Smat!"<<endl;
    exit(-1);
  }

  COMPLEX * Smat_copy = new COMPLEX[N*N];
  for(int i = 0; i < N * N; i++)
  {
    Smat_copy[i] = Smat_[i];
  }
  int LDA = N, LDVL = N, LDVR = N;
  COMPLEX w[N], vl[LDVL*N], vr[LDVR*N];
  LAPACKE_zgeev(LAPACK_ROW_MAJOR,'V', 'N', N, Smat_copy, LDA, w, vl, LDVL, vr, LDVR);

  Reorganize_U(N, vl, w);

  //note: vl is the transpose matrix of U(mixangle matrix)
  cout<<"U(the matrix of mixing angles): "<<endl;
  for(int i = 0; i < N; i++){
    for(int r = 0; r < N; r++){
      cout<<"        U["<< i * N + r <<"] = "<<vl[r * N + i]<<endl;
    }
  }

  COMPLEX ii = {0.0,1.0};
  for(int i = 0; i < N; i++){
    deltaArray_[i] = -ii / 2.0 * log(w[i]) * 180.0 /M_PI;
  }

  if(N == 2){
    Mixangle[0] = acos(vl[0]) * 180.0 / M_PI;
  } else{
    Mixangle[0] = atan(vl[7] / vl[8]) * 180.0 /M_PI;
    Mixangle[1] = asin(vl[6]) * 180.0 / M_PI;
    Mixangle[2] = atan(vl[3] / vl[0]) * 180.0 /M_PI;
  }
  delete [] Smat_copy;
}

//Reorganize U matrix and the matrix of eigenvalue to satisfy the convention given by Hubber.
//see Hubber 1995, page: 1104, the paragraph above RESULTS.
void Reorganize_U(int N_line, COMPLEX * U, COMPLEX * Ev){

  int maxLine[N_line] = {0};
  for(int i = 0; i < N_line; i++){
    for(int r = 0; r < N_line; r++){
      if(U[r * N_line + i].real() > U[maxLine[i] * N_line + i].real())
      maxLine[i] = r;
    }
  }
  bool cc = true;
  for(int i = 0; i < N_line; i++){
    if(maxLine[i] != i)
      cc = false;
  }
  if(!cc){
    if((N_line == 2 && maxLine[0] == maxLine[1]) || 
      (N_line ==3 && (maxLine[0] == maxLine[1] || maxLine[0] == maxLine[2] || maxLine[1] == maxLine[2]))){
        cout<<"U(The matrix of mixangles) can not be reorganized."<<endl;
        for(int i =0; i < N_line * N_line; i++){
          cout<<"        U["<<i<<"] = "<<U[i]<<endl;
        }
        exit(-1);
    }
    COMPLEX U_[N_line*N_line], Ev_[N_line];
    for(int i = 0; i < N_line; i++){
      Ev_[i] = Ev[i];
      for(int r = 0; r < N_line; r++)
        U_[i * N_line + r] = U[i * N_line +r];
    }
    for(int i = 0; i < N_line; i++){
      Ev[maxLine[i]] = Ev_[i]; 
      for(int r = 0; r < N_line; r++){
        U[r * N_line + maxLine[i]] = U_[r * N_line + i];
      }
    }
  }
}
