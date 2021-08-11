#include "mkl.h"
#include<iostream>
#include<cmath>
#include<omp.h>

#define m 1.0
#define V 2.0
#define P 3.0
#define L 4.0
#define N 5

double delta(double i, double j)
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

double A(double i, double j, double sigma)
{
    double sum=0,a;
    for(int k=1;k<=4*N;k++)
    {
        sum=sum+1/(pow(sigma,2)-pow((2*M_PI*k/L-P/3),2)-(2*M_PI*k/L-P/3)*(2*M_PI*i/L-P/3)-pow((2*M_PI*i/L-P/3),2));
    }
    a=delta(i,j)*(1-m*V/L*sum)-(2*m*V/L)/(pow(sigma,2)-pow((2*M_PI*j/L-P/3),2)-(2*M_PI*j/L-P/3)*(2*M_PI*i/L-P/3)-pow((2*M_PI*i/L-P/3),2));
    return a;
}

double det(double sigma)
{
    MKL_INT n=N,info;  //设置矩阵大小
    //MKL_INT info[N];
    MKL_INT ipiv[N];
    int numa=0,numb=0;
    double a[N*N],b=1;
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=N;j++)
        {
            a[numa]=A(i,j,sigma);
            //std::cout<<a[numa]<<" ";
            numa=numa+1;
        }
        //std::cout<<std::endl;
    }
    dgetrf(&n,&n,a,&n,ipiv,&info);
    for(int i=1;i<=N;i++)
    {
        for(int j=1;j<=N;j++)
        {
            std::cout<<" "<<a[numb];
            numb=numb+1;
        }
        std::cout<<std::endl;
    }
    
    return b;
}

int main()
{
    det(0.05);
}