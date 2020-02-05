#include<math.h>
#include<stdio.h>
#include<quadmath.h>

__float128 G=6.674q*powq(10,-11),M=5.974q*powq(10,24),m=7.348q*powq(10,22),R=3.844q*powq(10,8),omega=2.662q*powq(10,-6);

__float128 f(__float128 r)
{
    __float128 a;
    a=G*M/(r*r)-G*m/((R-r)*(R-r))-omega*omega*r;
    return a;
}

int main()
{
    __float128 a,c,b,xn_2,xn_1,xn,xplus,xnn,n=1010.0q,xn1,i;
    __float128 fxn_2,fxn_1,fxn,fxnn,fxn1,fxplus,RR;
    RR=10*R;
   for(i=2;i<=10000;i++)
   {
       fxn_2=f(1+RR*(i-2)/n);
       xn_2=1+RR*(i-2)/n;
       fxn_1=f(1+RR*(i-1)/n);
       xn_1=1+RR*(i-1)/n;
       fxn=f(1+RR*i/n);
       xn=1+RR*i/n;
       fxnn=fxn;
       xnn=xn;
       fxn1=f(1+RR*(i+1)/n);
       xn1=1+RR*(i+1)/n;
       fxplus=fxn1;
       xplus=xn1;
       //printf("%Qf\n",f(xn));
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
                printf("%Qf\n",xplus);
                xn_2=xn_1;
                xn_1=xn;
                xnn=xn;
                xn=xplus;
            }
            printf("%.20Qf\n",xplus);
       }
   }
}