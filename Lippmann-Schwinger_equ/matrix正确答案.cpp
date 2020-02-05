#include<stdio.h>
#include<math.h>
#include "lapacke.h"
//各种常数
int N=400;
double g=10.0;
double M0=1.0;
double C=30.0;
double lambda=35.0;
double v(double q_out,double q_in)
{
    double a;
    a=g*exp(-pow(q_out/lambda,4))*exp(-pow(q_in/lambda,4));
    return a;
}
double k(double i,double q)
{
    double a;
    if(i==N+1)
    {
        a=q;
    }
    else
    {
        double x;
        x=i*2.0/N-1;
        a=C*tan(M_PI*(x+1)/4);
    }
    return a;
}
double s(double i)
{
    double a;
    double x,w;
    x=i*2.0/N-1;
    w=2.0/N;
    a=C*(M_PI/4)*w/(pow(cos(M_PI*(x+1)/4),2));
    return a;
}
double delta(double i,double j)
{
    double a;
    if(i==j)
    {
        a=1;
    }
    else
    {
        a=0;
    }
    return a;
}
double u(double j,double q)
{
    double a=0;
    int i;
    if(j==N+1)
    {
        for(i=1;i<=N;i++)
        {   
            a=a+s(i)/(pow(k(i,q),2)-pow(q,2));
        }
        a=-a*M0*pow(q,2);
    }
    else
    {
        a=M0*pow(k(j,q),2)*s(j)/(pow(k(j,q),2)-pow(q,2));
    }
    return a;
}
double A(double i,double j,double q)
{
    double a;
    a=delta(i,j)+u(j,q)*v(k(i,q),k(j,q));
    return a;
}

int main()
{
    lapack_int n=N+1;
    lapack_int ipiv[N+1];
    int numa,numb;
    double q=2.0;      //这是我们要算的R(q,q)的q值
    double i,j;
    double a[(N+1)*(N+1)],b[N+1];
    numa=0;
    for(i=1;i<=N+1;i++)
    {
        for(j=1;j<=N+1;j++)
        {
            a[numa]=A(i,j,q);
            if(j==N+1)
            {
                //printf("%e %e",A(i,j,q),v(k(i,q),v(j,q)));
            }
            else
            {
                //printf("%e ",A(i,j,q));
            }
            
            numa=numa+1;
        }
        //printf("\n");
    }
    numb=0;
    for(i=1;i<=N+1;i++)
    {
        b[numb]=v(k(i,q),q);
        numb=numb+1;
        //printf("%e %e\n",v(k(i,q),q),k(i,q));
    }
    LAPACKE_dgesv(LAPACK_ROW_MAJOR,n,1,a,n,ipiv,b,1);
    printf("%.10e\n",b[N]);

}