#include<math.h>
#include<stdio.h>
#include<quadmath.h>
__float128 f(__float128 x)
{
    __float128 a,h=0.5q;
    a=x*tanq(x)-sqrtq(powq(h,2)-powq(x,2));
    a=a*cosq(x);
    return a;
}
__float128 g(__float128 x)
{
    __float128 a,h=1.0q;
    a=x*(1/tanq(x))+sqrtq(powq(h,2)-powq(x,2));
    a=a*sinq(x);
    return a;
}
int main()
{
    __float128 a,c,b,xn_2,xn_1,xn,xplus,xnn,n=1000.0q,xn1,i;
    __float128 fxn_2,fxn_1,fxn,fxnn,fxplus,fxn1;

   for(i=2;i<=1000;i++)
   {
       fxn_2=f(-0.5q+(i-2)/n);
       xn_2=-0.5q+(i-2)/n;
       fxn_1=f(-0.5q+(i-1)/n);
       xn_1=-0.5q+(i-1)/n;
       fxn=f(-0.5q+i/n);
       xn=-0.5q+i/n;
       fxnn=fxn;
       xnn=xn;
       fxn1=f(-0.5q+(i+1)/n);
       xn1=-0.5q+(i+1)/n;
       fxplus=fxn1;
       if(fxn*fxn1<0)
       {
           while(abs(xnn-xplus)>=powq(10,-11))
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
            printf("%.20Qf\n",xplus);
       }
   }
    
    for(i=2;i<=1000;i++)
   {
       xn_2=g(-1.0q+2*(i-2)/n);
       xn_1=g(-1.0q+2*(i-1)/n);
       xn=g(-1.0q+2*i/n);
       xnn=xn;
       xn1=g(-1.0q+2*(i+1)/n);
       xplus=xn1;
       if(xn*xn1<0)
       {
           while(abs(xnn-xplus)>=powq(10,-11))
            {
                c=g(xn);
                a=((xn_1-xn)*(g(xn_2)-g(xn))-(xn_2-xn)*(g(xn_1)-g(xn)))/((xn_2-xn_1)*(xn_2-xn)*(xn_1-xn));
                b=((xn_2-xn)*(xn_2-xn)*(g(xn_1)-g(xn))-(xn_1-xn)*(xn_1-xn)*(g(xn_2)-g(xn)))/((xn_2-xn_1)*(xn_2-xn)*(xn_1-xn));

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
            printf("%.20Qf\n",xplus);
       }
   }
}