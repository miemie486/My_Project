#include <stdio.h>
#include <math.h>
#include <quadmath.h>
int main()
{
    double n=500000000;
    double pi=0,pi1=0,l1=1,l=-1,i,ii,pi2,pi3;
    for(i=1;i<=n;i++)
    {
        l=l*(-1);
        pi=pi+l*(1/(2*i-1));
    }
    printf("正着算的pi=%.100f\n",pi*4);
    printf("与math.h内储存值之间的误差：%.100f\n",pi*4-M_PI);
    for(ii=n;ii>=1;ii--)
    {
        l1=l1*(-1);
        pi1=pi1+l1*(1/(2*ii-1));
    }
    printf("倒着算的pi=%.100f\n",pi1*4);
    printf("与math.h内储存值之间的误差：%.100f\n",pi1*4-M_PI);
}
