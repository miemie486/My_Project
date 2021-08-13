// Calculation of wignerJSymbol (3, 6, 9 j symbols)
// The source of 6-j's exact formula is lost. Yet all 3nj simbols have been checked by http://dentonwoods.nfshost.com/3j6j9j.html, thus follow the convention in Wiki. CG coefficient's formula comes from formula (7) in http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html. The formula is same as which in Wiki.

#include "3nj.h"
#include "def.h"
//
using namespace std;

unsigned long long int facArray[_MAXN_FACTBL_] {0};

/* This is potentially dangerous. For in multi-threading environment, the
initialization may be done twice.*/
void facLabel()
{
  // Generating factals from 0! to 20!
  static bool initialized {false};
  if (initialized) return;
  facArray[0]=1;
  for(short i=1; i<_MAXN_FACTBL_; i++){
    facArray[i]=facArray[i-1]*i;
  }
  initialized = true;
}

void ifTri(REAL a, REAL b, REAL c){
  // Testing if it is a triangle. a can be 0 if fabs(a+b+c - int(a+b+c) < _3NJ_EPS_, a + b + c are legal as a set of quantum numbers.
  if(a<= b+c && b<= a+c && c<= a+b && fabs(a+b+c - int(a+b+c)) < _3NJ_EPS_){}
  else {
    printf("SixJ: This is not a triangle\n%f, %f, %f\n", a, b, c);
    abort();
  }
}

REAL delta(REAL a, REAL b, REAL c){
  // delta here is sqrt(delta) in http://mathworld.wolfram.com/TriangleCoefficient.html
  return sqrt(1.*facArray[int(a+c-b)]/facArray[int(a+b+c+1)]*facArray[int(a+b-c)]*facArray[int(b+c-a)]);
}


REAL threeJ(REAL j1, REAL j2, REAL j3, REAL m1, REAL m2, REAL m3){
  /* This in fact is threeJ symbol in special case which is convenient for calculating CG coefficient and hence the G factor. Therefore you might see more restriction than usual.

     error = 0 ==> no error
     error = 1 ==> error
  */
  // formulars (7) in http://mathworld.wolfram.com/Wigner3j-Symbol.html
  // Conditions https://en.wikipedia.org/wiki/3-j_symbol selection rule
  ifTri(j1, j2, j3);
  if(m1>=-fabs(j1) && m1<= fabs(j1) && m2>=-fabs(j2) && m2<=fabs(j2) && m3>=-fabs(j3) && m3<=fabs(j3) && fabs(m1+m2+m3)<_3NJ_EPS_ && fabs(j1+m1-int(j1+m1))<_3NJ_EPS_ && fabs(j2+m2-int(j2+m2))<_3NJ_EPS_ && fabs(j3+m3-int(j3+m3))<_3NJ_EPS_ ){}
  else {
    std::cout<<"incorrect m: each m must be set within [-j,j]) and m1+m2+m3 must be equal to 0 and j-m must be equal to an integer."<<std::endl;
    abort();
       };
    REAL sum1=delta(j1,j2,j3)*sqrt(facArray[int(j1+m1)]*facArray[int(j1-m1)]*facArray[int(j2+m2)]*facArray[int(j2-m2)]*facArray[int(j3+m3)]*facArray[int(j3-m3)]);
    REAL sum2 {0};
    REAL mint=max(REAL(max(REAL(0),j2-j3-m1)),j1-j3+m2);
    REAL maxt=min(min(j1+j2-j3,j1-m1),j2+m2);
    for(REAL t=mint;t<=maxt;t++)
    sum2+=1.*pow(-1,t)/(facArray[int(t)]*facArray[int(j3-j2+t+m1)]*facArray[int(j3-j1+t-m2)]*facArray[int(j1+j2-j3-t)]*facArray[int(j1-t-m1)]*facArray[int(j2-t+m2)]);
    if(int(fabs(j1-j2-m3))%2==0)
     {return sum1*sum2;}
    else
      return -sum1*sum2;
}


double CGcoef(double j1, double j2,double j3, double m1,double m2,double m3){
  // Formula (7) in http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html
  double sum {-1};
  if(m1>=-fabs(j1) && m1<= fabs(j1) && m2>=-fabs(j2) && m2<=fabs(j2) && m3>=-fabs(j3) && m3<=fabs(j3) && fabs(m1+m2-m3)<_3NJ_EPS_ && fabs(j1+m1-int(j1+m1))<_3NJ_EPS_ && fabs(j2+m2-int(j2+m2))<_3NJ_EPS_ && fabs(j3+m3-int(j3+m3))<_3NJ_EPS_ ){}
  else {
    std::cout<<"incorrect m: each m must be set within [-j,j]) and m1+m2 must be equal to m3 and j-m must equal an integer."<<std::endl;
    abort();
       };
  if(int(fabs(j1-j2+m3)) % 2==0){sum=1;}
  return sum*threeJ(j1, j2, j3, m1, m2, -m3)*sqrt(2*j3+1);
}

REAL threeJm0(REAL j1, REAL j2, REAL j3){
  /* This in fact is threeJ symbol in special case which is convenient for calculating CG coefficient and hence the G factor. Therefore you might see more restriction than usual.

     error = 0 ==> no error
     error = 1 ==> error
  */
  // formulars (15) in http://mathworld.wolfram.com/Wigner3j-Symbol.html
  // Conditions https://en.wikipedia.org/wiki/3-j_symbol selection rule
  int J=j1+j2+j3;
  int halfJ= J/2;
  REAL s {0};

  ifTri(j1, j2, j3);
  if(fabs(int(j1) - j1) > _3NJ_EPS_ || fabs(int(j2) - j2) > _3NJ_EPS_ || fabs(int(j3) - j3) > _3NJ_EPS_) printf("error in 3j");
  if(J==halfJ*2){
    s=sqrt(1.*facArray[int(J-2*j1)]/facArray[J+1]*facArray[int(J-2*j2)]*facArray[int(J-2*j3)])*facArray[halfJ]/facArray[int(halfJ-j1)]/facArray[int(halfJ-j2)]/facArray[int(halfJ-j3)];
    if(halfJ %2==0){return s;}
    else {return -s;}
  }
  else {return 0;}
}
double CGm0(double j1, double j2, double j3){
  // Formula (7) in http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html, note that m = 0
  double sum {-1};

  if(int(abs(j1-j2)) % 2==0){sum=1;}
  return sum*threeJm0(j1, j2, j3)*sqrt(2*j3+1);
}




////////////////////////// 6j symbol//////////////////////////
////////////////////////// Check input ///////////////////////
void sixJChecker(REAL a, REAL b, REAL c, REAL d, REAL e, REAL f){
  // Check if a sixJ matrix is legal. P324
  // {{a, b, c}, {d, e, f}}
  ifTri(a, b, c);
  ifTri(d, b, f);
  ifTri(d, e, c);
  ifTri(a, e, f);
  if(a + b + c - (int) (a + b + c) > _3NJ_EPS_ || b + d + f - (int) (b + d + f) > _3NJ_EPS_ || a + e + f - (int) (a + e + f) > _3NJ_EPS_ || d + e + c - (int) (d + e + c) > _3NJ_EPS_) {cout<< "sixj: violate with the int rule"<< endl; abort();}
    // int rule see: http://mathworld.wolfram.com/Wigner6j-Symbol.html
}


///////////////////////// Evaluate ///////////////////////////
REAL sixJ(REAL j1, REAL j2, REAL j3, REAL J1, REAL J2, REAL J3){
//See formula (4) in http://mathworld.wolfram.com/Wigner6j-Symbol.html
  sixJChecker(j1, j2, j3, J1, J2, J3);

  REAL zmin=max(max(max(j1+j2+j3,j1+J2+J3),J1+j2+J3),J1+J2+j3);

  REAL zmax=min(min(j1+j2+J1+J2,j2+j3+J2+J3),j3+j1+J3+J1);

  REAL sum1 {0}, sum2 {0},ft {0};

  sum1=delta(j1, j2, j3)*delta(j1, J2, J3)*delta(J1, j2, J3)*delta(J1, J2, j3);

  for(REAL t=zmin; t<=zmax; t++){

  ft=facArray[int(t-j1-j2-j3)]*facArray[int(t-j1-J2-J3)]*facArray[int(t-J1-j2-J3)]*facArray[int(t-J1-J2-j3)]*facArray[int(j1+j2+J1+J2-t)]*facArray[int(j2+j3+J2+J3-t)]*facArray[int(j3+j1+J3+J1-t)];
  sum2+=1.*pow(-1,t)*facArray[int(t+1)]/ft;

  }
   return sum2*sum1;

}

///////////////////////// Evaluate ///////////////////////////
REAL oldsixJ(REAL e, REAL a, REAL f, REAL b, REAL d, REAL c){
//This is a sixJ function in old version whose source could not find.
  sixJChecker(e, a, f, b, d, c);
  // {{e, a, f}, {b, d, c}}={{j1, j2, j12}, {j3, j, j23}}
  REAL zMax=min(min(2*b, -a+b+c), b-d+f), z {0}, index {1};
  REAL sum1 {0}, sum2 {0}, an {0};

  sum1=delta(a, b, c)*delta(a, e, f)*delta(c, d, e)*delta(b, d, f)*facArray[int(a+b+c+1)]/facArray[int(a+b-c)]*facArray[int(b+d+f+1)]/(facArray[int(c-d+e)]*facArray[int(c+d-e)]*facArray[int(a-e+f)]*facArray[int(-a+e+f)]*facArray[int(b+d-f)]);

  for(z=0; z<zMax+1; z++){
  an=1.*facArray[int(2*b-z)]/facArray[int(-a+b+c-z)]*facArray[int(b+c-e+f-z)]/facArray[int(b-d+f-z)]*facArray[int(b+c+e+f+1-z)]/facArray[int(a+b+c+1-z)]/facArray[int(b+d+f+1-z)]/facArray[int(z)];
    sum2+=index*an;
    index*=-1;
  }
  if(int(b+c+e+f)%2==0) return sum2*sum1;
  else return -sum2*sum1;
}

///////////////////////// 9j symbol////////////////////////

REAL nineJ(REAL j1, REAL j2, REAL j12, REAL j3, REAL j4, REAL j34, REAL j13, REAL j24, REAL j){
  // GaoDengLiangZiLiXue 2nd Edition (XingLin Ke) P325 (23.84), which is same as formula (3) in http://mathworld.wolfram.com/Wigner9j-Symbol.html
  REAL s {0}, sign {1};
  REAL min[4] {3}, max[3] {0};
  REAL minb, maxb;

  min[0]=abs(j1-j);
  min[1]=abs(j3-j24);
  min[2]=abs(j2-j34);
  minb = min[0];
  for(int idx = 1; idx < 3; idx++)
    if(minb < min[idx])
      minb = min[idx];

  max[0]=j1+j;
  max[1]=j3+j24;
  max[2]=j2+j34;
  maxb = max[0];
  for(int idx = 1; idx < 3; idx++)
    if(maxb > max[idx])
      maxb = max[idx];

  if(minb != int(minb)){sign=-1;}
  for(REAL jp=minb; jp<=maxb; jp++){
    s+=(2*jp+1)*sixJ(j1, j2, j12, j34, j, jp)*sixJ(j3, j4, j34, j2, jp, j24)*sixJ(j13, j24, j, jp, j1, j3);
  }

  return s*sign;
}

REAL test9j23_87(REAL a, REAL b, REAL c, REAL d, REAL e, REAL f){
  REAL sign {1};

  if(int(b+c+e+f) % 2!=0){sign=-1;}
  return sign/sqrt((2*e+1)*(2*f+1))*sixJ(a, b, e, d, c, f);
}
