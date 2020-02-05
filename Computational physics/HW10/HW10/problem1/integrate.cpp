#include<stdio.h>
#include<math.h>
double rho=6.022*pow(10,28),thetaD=428,kB=1.38065*pow(10,-23),V=pow(10,-3),e=M_E;
double f(double x)
{
    double a;
    if(x==0)
    {
        a=0;
    }
    else
    {
        a=pow(x,4)*pow(e,x)/(pow(pow(e,x)-1,2));
    }
    
    return a;
}
double C(double T)
{
    int i;
    double a,b,n=50,h;
    a=9*V*rho*kB*pow(T/thetaD,3);
    h=thetaD/(n*T);
    b=(h/3)*(f(0)+f(thetaD/T));
    for(i=1;i<=n/2;i++)
    {
        b=b+(4*h/3)*f((2*i-1)*h)+(2*h/3)*f(2*i*h);
    }
    b=b-f(thetaD/T)*(2*h/3);
    return a*b;
    
}
int main()
{
    int i;
    double x0=5,x1=500,n=4500,x;
    FILE *fp;
    fp=fopen("data.txt","w+");
    for(i=0;i<=n;i++)
    {
        x=5+i*(x1-x0)/n;
        fprintf(fp,"%f ",x);
    }
    fprintf(fp,"\n");
    for(i=0;i<=n;i++)
    {
        x=5+i*(x1-x0)/n;
        fprintf(fp,"%f ",C(x));
    }
    fcloseall();

}