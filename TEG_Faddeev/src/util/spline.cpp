//Method from
//Numerical Treatment of Few Body Equations in Momentum Space by the Spline Method
//W. Glockle, G. Hasberg and A.R. Neghabian
//Z.Phy. A - Atoms and Nuclei 305,217-221(1982)

#include "spline.h"

using namespace std;

inline REAL delta(int i, int j)
{
  if (i == j)
    return (REAL)1;
  else
    return (REAL)0;
}

inline REAL B(int j, int i, REAL *h){
  return 6 * delta(i, j + 1) / (h[j] + h[j + 1]) / h[j + 1]
    - 6 * delta(i, j) / h[j] / h[j + 1]
    + 6 * delta(i, j - 1) / (h[j] + h[j + 1]) / h[j];
}

inline REAL A(int j, int i, REAL *p, REAL *h, REAL *mu){
  if (j == 0) //using A[0,i]=0;
    return 0;
  else
    return B(j, i, h) / p[j] - mu[j] * A(j - 1, i, p, h, mu) / p[j];
}

inline REAL c(int j, int i, REAL *q, REAL *mu, int meshNum, REAL *p, REAL *h)
{
  if (j == meshNum-1)
    return 0;
  else
    return q[j] * c(j + 1, i, q, mu, meshNum, p, h) + A(j, i, p, h, mu);
}

//using initial value ini_q=0
inline void ini_p_q(REAL *p, REAL *q, REAL *mu, REAL *lam, REAL meshNum){
  p[0] = 2;
  q[0] = 0;
  for (int i(1); i < meshNum - 1; i++) {
    p[i] = mu[i] * q[i - 1] + 2;
    q[i] = -lam[i] / p[i];
  }
}

inline void ini_h(REAL *h, REAL *x, int meshNum)
{
  for (int i(1); i < meshNum; i++) {
    h[i] = x[i] - x[i - 1];
  }
}

inline void ini_lam(REAL *lam, REAL *h, int meshNum)
{
  lam[0] = 0;
  for (int i (1); i < meshNum - 1; i++) {
    lam[i] = h[i + 1] / (h[i] + h[i + 1]);
  }

}

inline void ini_mu(REAL *mu, REAL *lam, int meshNum)
{
  for (int i(1); i < meshNum - 1; i++) {
    mu[i] = 1 - lam[i];
  }
}

void genSplineCoefs(REAL *coe1, REAL *coe2, REAL *coe3, REAL *x, int meshNum){
  /* Generate coefficients for calculating spline
     coe1, coe2, coe3: arrays with length (meshNum - 1) * meshNum
     x: mesh points
     meshNum: number of mesh points
   */

  REAL *h = new REAL[meshNum];
  REAL *lam = new REAL[meshNum - 1];
  REAL *mu = new REAL[meshNum - 1];
  REAL *p = new REAL[meshNum - 1];
  REAL *q = new REAL[meshNum - 1];

  ini_h(h, x, meshNum);
  ini_lam(lam, h, meshNum);
  ini_mu(mu, lam, meshNum);
  ini_p_q(p, q, mu, lam, meshNum);

  REAL cji {0}, cj1i {0};

  for (int i = 0; i < meshNum; i++){
      for (int j = 0; j < meshNum - 1; j++){
        cji = c(j, i, q, mu, meshNum, p, h);
        cj1i = c(j + 1, i, q, mu, meshNum, p, h);

        coe1[i * (meshNum - 1) + j] = (delta(i, j + 1) - delta(i, j)) / h[j + 1] - h[j + 1] / 6 * (2 * cji + cj1i);
        coe2[i * (meshNum - 1) + j] = cji / 2;
        coe3[i * (meshNum - 1) + j] = (cj1i - cji )/ h[j + 1] / 6;
      }
    }
}

void splineArray(REAL *x, REAL var, int meshNum, REAL *Si, REAL *coe1, REAL *coe2, REAL *coe3){

  //if (var<x[0] || var>x[meshNum-1]) {
  //  std::cout << "my_error" << std::endl;
  //}
  if (var<x[0] || var>x[meshNum-1]){
    for (int i = 0; i < meshNum; i++) {
        Si[i] = 0;
      }
    return;
  }

  int j {0};

  for (int i = 0; i < meshNum - 1; i++){
    if (var <= x[i + 1] && var >= x[i] && i != 0)
      {
        j = i;
        break;
      }
  }

  REAL con = var - x[j];
  int idx = j;

  for (int i = 0; i < meshNum; i++)
    {
      Si[i] = delta(i, j) + con * (coe1[idx] + con * (coe2[idx] + con * coe3[idx]));
      idx += (meshNum - 1);
    }
}

REAL spline(REAL *x, REAL var, int meshNum, REAL *coe1, REAL *coe2, REAL *coe3, int iRequire){

  int j {0};

  if (var < x[0] || var > x[meshNum-1]) {
    return 0;
  }

  for (int i = 0; i < meshNum - 1; i++){
    if (var <= x[i + 1] && var >= x[i] && i != 0) {
        j = i;
        break;
      }
  }

  REAL con = var - x[j];
  int idx = j + iRequire * (meshNum - 1);
  REAL ans = delta(iRequire, j) + con * (coe1[idx] + con * (coe2[idx] + con * coe3[idx]));

  return ans;
}

splineFunc::splineFunc(size_t lenx, size_t lenxi){
  n = lenx;
  ni = lenxi;
  coe1 = new REAL[n * (n - 1)];
  coe2 = new REAL[n * (n - 1)];
  coe3 = new REAL[n * (n - 1)];
  xi = new REAL[ni];
  yi = new REAL[ni];
  x = new REAL[n];
  y = new REAL[n];
  for(size_t idx = 0; idx < n; idx++){
    x[idx] = 0;
    y[idx] = 0;
  }
  for(size_t idx = 0; idx < n * (n - 1); idx++){
    coe1[idx] = 0;
    coe2[idx] = 0;
    coe3[idx] = 0;
  }
  for(size_t idx = 0; idx < ni; idx++){
    xi[idx] = 0;
    yi[idx] = 0;
  }
}

splineFunc::~splineFunc(){
  delete[] coe1;
  delete[] coe2;
  delete[] coe3;
  delete[] x;
  delete[] y;
  delete[] xi;
  delete[] yi;
};

void splineFunc::printArr(REAL*array, size_t len){
  for(size_t idx = 0; idx < len; idx++) printf("%.3e\t", array[idx]);
  printf("\n");
};

void splineFunc::genYFun(){
  for(size_t idx = 0; idx < n; idx++)y[idx] = func(x[idx]);
}

void splineFunc::genYIFun(void){
  REAL sj {0};

  genSplineCoefs(coe1, coe2, coe3, x, n);
  for (size_t i = 0; i < ni; i++) {
  	for (size_t j = 0; j < n; j++) {
      sj = spline(x, xi[i], n, coe1, coe2, coe3, j);
  		yi[i] = yi[i] + sj * y[j];
  	}
  }
}

void splineFunc::plotTable(ofstream& myfile){
  /* output all data into file
     input:
     str: name and location of the file. .csv file is recommanded.
     output:
     file with name str. The first row is x, second row is y, third row is xi, final row is yi.
  */
  for(size_t idx = 0; idx < n - 1; idx++) myfile << x[idx] << ",";
  myfile << x[n - 1] << endl;
  for(size_t idx = 0; idx < n - 1; idx++) myfile << y[idx] << ",";
  myfile << y[n - 1] << endl;
  for(size_t idx = 0; idx < ni - 1; idx++) myfile << xi[idx] << ",";
  myfile << xi[ni - 1] << endl;
  for(size_t idx = 0; idx < ni - 1; idx++) myfile << yi[idx] << ",";
  myfile << yi[ni - 1] << endl;
}
