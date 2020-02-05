#include<stdio.h>
#include<quadmath.h>
#include<time.h>
__float128 f(__float128 x)
{
    __float128 a;
    a=powq(sinq(sqrtq(100.0q*x)),2);
    return a;
}
int main()
{
    double start,end;
    __float128 n,k,i;
    __float128 Rk,rk1=0,h,error=1;
    k=1;
    start=clock();
    //while(error>powq(10.0q,-10.0q))
    while(fabsq(0.45583253230908513731596260144464q-Rk)>powq(10.0q,-10.0q))
    {
        h=1/powq(2.0q,k);
        Rk=(h/3.0q)*(f(0.0q)+f(1.0q));
        for(i=1;i<=powq(2.0q,k-1);i++)
        {
            Rk=Rk+(4.0q*h/3.0q)*f((2.0q*i-1.0q)*h)+(2.0q*h/3.0q)*f(2.0q*i*h);
        }
        Rk=Rk-(2.0q*h/3.0q)*f(1);
        //error=7000000/(180*powq(powq(2.0q,k),4));
        printf("%.0Qf\n",powq(2.0q,k));
        printf("error%.32Qf\n",0.45583253230908513731596260144464q-Rk);
        rk1=Rk;
        //printf("error=%.32Qf\n",error);
        printf("R%.0Qf=%.32Qf\n",k,Rk);
        k=k+1;
    }
    end=clock();
    printf("%f\n",(end-start)/(CLOCKS_PER_SEC));

}