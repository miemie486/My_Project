/*
  Generating w_n and x_n in Gauss Quadrature.
  Ref: https://en.wikipedia.org/wiki/Gaussian_quadrature
*/

#include "gaussIntegral.h"

void pl(int n, REAL x, REAL p_array[2]){
  /* Evaluate the value of Legendre polynomials with order n and n-1 by using recursion relations[1].
     Ref:
     [1] First formular in Recursion, Wikipedia https://en.wikipedia.org/wiki/Legendre_polynomials.

     Input:
     n: order of Legendre polynomial.
     x: p_(n)(x).

     Output:
     p_array[0]=p_(n-1)(x): Legendre polynomials with order n-1.
     p_array[1]=p_(n)(x): Legendre polynomials with order n.

  */
  REAL t;  // A temperory variable to restore value of p_n(x).
  int j;     // Recursional index.

  if(n<2) {  // p_0(x)= 1, p_1(x)=x.
      p_array[0]=1;
      p_array[1]=x;
    }
  else {
      p_array[0]=1;
      p_array[1]=x;
      for(j=2;j<=n;j++) {
          t=p_array[1];
          p_array[1]=((2*j-1)*x*p_array[1]-(j-1)*p_array[0])/j;
          p_array[0]=t;
        }
    }
}

REAL root(int n, int i){
  /* Find ith root of function p_n(x).

     Input:
     n: order of Legendre polynomial.
     i: ith root of function p_n(x)

     Output:
     return: ith root of function p_n(x)
   */

  REAL x0,x1,dp, p_array[2] {0};

  x1=cos(M_PI*(i-0.25)/(n+0.5));
  x0=0;
  while(fabs(x1-x0) > _GAULEG_EPS_)
    {
      x0=x1;
      pl(n,x1,p_array);
      dp=n*(x1*p_array[1]-p_array[0])/(x1*x1-1);
      x1=x1-p_array[1]/dp;
    }
  return x1;
}

void createGauss(REAL xn[], REAL wn[], int n, REAL xmin, REAL xmax){
  /*
    Creat x_n and w_n in first table, Wikipedia https://en.wikipedia.org/wiki/Legendre_polynomials. http://mathfaculty.fullerton.edu/mathews/n2003/GaussianQuadMod.html. If xmax > xmin, x_n and w_n are mesh points and weight in interval (xmin, xmax). If xmax == xmin, they are in interval (0, +inf) with C = xmax (normally take C = 1000 MeV).

    Input:
    n: number of mesh points. if(n == odd) n+=1;
    xmin, xmax: x belongs to [xmin, xmax].

    Output:
    xn[]: x_n.
    wn[]: w_n.

  */

  REAL dp[n], p_array[2] {0};
  REAL d1, d2;
  int i;

  //Find roots of Legendre Poly with order n
  for(i=0;i<n/2;i++) {
    xn[i]=root(n,i+1);
    pl(n,xn[i],p_array);
    dp[i]=-n*p_array[0]/(xn[i]*xn[i]-1);
    wn[i]=2./(1-xn[i]*xn[i])/(dp[i]*dp[i]);
    wn[n-i-1]=wn[i];
    xn[n-i-1]=-xn[i];
  }

  if(n % 2 != 0){
    i = n/2;
    xn[i]=0;
    pl(n,0,p_array);
    dp[i]=n*p_array[0];
    wn[i]=2./(dp[i]*dp[i]);
  }

  d1=(xmax-xmin)/2;
  d2=(xmax+xmin)/2;

  for(int i=0; i< n; i++){
    xn[i]=xn[i]*d1+d2;
    wn[i]=wn[i]*d1;
  }

  REAL tmp;

  /* Reverse the order of xn[] and wn[] */
  for(int i=0; i< n/2; i++){
    tmp = xn[i];
    xn[i] = xn[n-1-i];
    xn[n-1-i] = tmp;
    tmp = wn[i];
    wn[i] = wn[n-1-i];
    wn[n-1-i] = tmp;
  }

}

/* From Numerical Recipes */
// void createGauss(REAL xn[], REAL wn[], int n, REAL xmin, REAL xmax)
// {
//   int m,j,i;
//   REAL z1,z,xm,xl,pp,p3,p2,p1;

//   m=(n+1)/2;
//   xm=0.5*(xmax+xmin);
//   xl=0.5*(xmax-xmin);
//   for (i=0;i<m;i++) {
//     z=cos(M_PI*(i+0.75)/(n+0.5));
//     do {
//       p1=1.0;
//       p2=0.0;
//       for (j=0;j<n;j++) {
//         p3=p2;
//         p2=p1;
//         p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
//       }
//       pp=n*(z*p1-p2)/(z*z-1.0);
//       z1=z;
//       z=z1-p1/pp;
//     } while (abs(z-z1) > _GAULEG_EPS_);
//     xn[i]=xm-xl*z;
//     xn[n-1-i]=xm+xl*z;
//     wn[i]=2.0*xl/((1.0-z*z)*pp*pp);
//     wn[n-1-i]=wn[i];
//   }
// }

void createGaussInf(REAL xn[], REAL wn[], int n, REAL tanC){
  /*
    Creat gaussian mesh points from 0 to inf

    Input:
    n: number of mesh points. if(n == odd) n+=1;
    tanC: para C

    Output:
    xn[]: x_n.
    wn[]: w_n.

  */

  REAL aux_xn[n], aux_wn[n];
  createGauss(aux_xn, aux_wn, n, -1.0, 1.0);

  for(int i=0; i< n; i++){
    wn[i] = tanC * M_PI * aux_wn[i] / (4.0 * pow(cos(0.25*M_PI*(aux_xn[i] + 1.0)), 2));
    xn[i] = tanC * tan(0.25*M_PI*(aux_xn[i] + 1.0));
  }

}
