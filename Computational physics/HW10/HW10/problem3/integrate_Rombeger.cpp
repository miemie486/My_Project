#include<stdio.h>
#include<math.h>
#include<time.h>
double f(double x)
{
    double a;
    a=1/(1+pow(x,2));
    return a;
}
int main()
{
    double start,end;
    start=clock();
    double a=0,b=1,R[255][255],h;
    int i=1,j=1,k;
    R[1][1]=0.5*(f(a)+f(b));
    while(fabs(R[i][i]-M_PI_4f64)>=pow(10,-10))
    {
        i=i+1;
        h=1/pow(2,i-1);
        R[i][1]=R[i-1][1];
        for(k=1;k<=pow(2,i-2);k++)
        {
            R[i][1]=R[i][1]+2*h*f(a+(2*k-1)*h);
        }
        R[i][1]=0.5*R[i][1];
        for(j=2;j<=i;j++)
        {
            R[i][j]=R[i][j-1]+(R[i][j-1]-R[i-1][j-1])/(pow(4,j-1)-1);
        }
        printf("%.16f %.16f\n",h,fabs(R[i][i]-M_PI_4f64));
    }


    end=clock();
    //printf("%f\n",(end-start)/(CLOCKS_PER_SEC));
}