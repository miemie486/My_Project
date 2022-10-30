#include <stdio.h>
#include <quadmath.h>
__float128 arctan(__float128 x)
{
    double i;
    __float128 n=50000;
    __float128 arc,l,o;
    l=(-1);
    o=x;
    arc=0.0q;
    char buf[256],buf1[256],buf2[256];
    for(i=1;i<=n;i++)
    {
       arc=arc+o*powq(10,100);
       o=-1*o*x*x*(2*i-1)/(2*i+1);
    }
    return arc;
}
int main()
{
    __float128 pi2,x1=0.2q,x2=1/239.0q,pi3;
    char bbuf[200];
    pi2=(4*arctan(x1)-arctan(x2))*4;
    pi3=(4*atanq(x1)-atanq(x2))*4;
    pi2=pi2-3*powq(10,100);
    pi2=rintq(pi2);
    quadmath_snprintf(bbuf,sizeof bbuf,"%.0Qf",100,pi2);
    printf("第二题的pi=3.%s\n",bbuf);
    printf("库里算的pi=%.100Qf\n",pi3);
}