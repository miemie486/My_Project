#include<stdio.h>
#include<math.h>
#include<quadmath.h>

__float128 f(__float128 x)
{
    __float128 a;
    a=924*powq(x,6)-2772*powq(x,5)+3150*powq(x,4)-1680*powq(x,3)+420*x*x-42*x+1;
    return a;
}

int main()
{
    __float128 a,c,b,xn_2,xn_1,xn,xplus,xnn,n=1010.0q,xn1,i;
    __float128 fxn_2,fxn_1,fxn,fxnn,fxn1,fxplus,RR;
   for(i=2;i<=10000;i++)
   {
       fxn_2=f((i-2)/n);
       xn_2=(i-2)/n;
       fxn_1=f((i-1)/n);
       xn_1=(i-1)/n;
       fxn=f(i/n);
       xn=i/n;
       fxnn=fxn;
       xnn=xn;
       fxn1=f((i+1)/n);
       xn1=(i+1)/n;
       fxplus=fxn1;
       xplus=xn1;
       if(fxn*fxn1<0)
       {
           while(abs(xnn-xplus)>=powq(10,-15))
            {
                c=f(xn);
                a=((xn_1-xn)*(f(xn_2)-f(xn))-(xn_2-xn)*(f(xn_1)-f(xn)))/((xn_2-xn_1)*(xn_2-xn)*(xn_1-xn));
                b=((xn_2-xn)*(xn_2-xn)*(f(xn_1)-f(xn))-(xn_1-xn)*(xn_1-xn)*(f(xn_2)-f(xn)))/((xn_2-xn_1)*(xn_2-xn)*(xn_1-xn));

                if(b<0)
                {
                    xplus=xn+(2*c/(-b+sqrtq(powq(b,2)-4*a*c)));
                }
                else
                {
                    xplus=xn-(2*c/(b+sqrtq(powq(b,2)-4*a*c)));
                }
                xn_2=xn_1;
                xn_1=xn;
                xnn=xn;
                xn=xplus;
            }
            printf("%.15Qf\n",xplus);
       }
   }
}