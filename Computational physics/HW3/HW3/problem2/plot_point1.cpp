#include<stdio.h>
#include<math.h>
#include<quadmath.h>

float valuef1(float x)
{
    float b;
    b=pow(x-1,9);
    b=b;
    return b;
}
float valuef2(float x)
{
    float b;
    b=pow(x,9)-9*pow(x,8)+36*pow(x,7)-84*pow(x,6)+126*pow(x,5)-126*pow(x,4)+84*pow(x,3)-36*pow(x,2)+9*x-1;
    b=b;
    return b;
}
double valued1(double x)
{
    double b;
    b=pow(x-1,9);
    b=b;
    return b;
}
double valued2(double x)
{
    double b;
    b=pow(x,9)-9*pow(x,8)+36*pow(x,7)-84*pow(x,6)+126*pow(x,5)-126*pow(x,4)+84*pow(x,3)-36*pow(x,2)+9*x-1;
    b=b;
    return b;
}
__float128 valueq1(__float128 x)
{
    __float128 b;
    b=powq(x-1,9);
    b=b;
    return b;
}
__float128 valueq2(__float128 x)
{
    __float128 b;
    b=b=powq(x,9)-9*powq(x,8)+36*powq(x,7)-84*powq(x,6)+126*powq(x,5)-126*powq(x,4)+84*powq(x,3)-36*powq(x,2)+9*x-1;
    return b;
}
int main()
{
    int i;
    float x1,v1;
    double x2,v2;
    __float128 x3,v3;
    FILE *fp11,*fp12;
    FILE *fp21,*fp22;
    FILE *fp31,*fp32;
    fp11=fopen("datafloat1.txt","w+");
    fp12=fopen("datafloat2.txt","w+");
    fp21=fopen("datadouble1.txt","w+");
    fp22=fopen("datadouble2.txt","w+");
    fp31=fopen("dataquad1.txt","w+");
    fp32=fopen("dataquad2.txt","w+");
    float qwe1=1-0.3,ewq1=0.6,n1=2000.0;//计算范围为-0.01至0.3，计算2000个点。
    for(i=0;i<=2000;i++)//float合起来算的
    {
        x1=qwe1+i*ewq1/n1;
        v1=valuef1(x1);
        fprintf(fp11,"%.100f %.100f\n",x1,v1);
    };
    for(i=0;i<=2000;i++)//float分开算的
    {
        x1=qwe1+i*ewq1/n1;
        v1=valuef2(x1);
        fprintf(fp12,"%.100f %.100f\n",x1,v1);
    };
    //sdfgahkfbahj//
    //dsgsagsdgsdg//
    double qwe2=1-0.3,ewq2=0.6,n2=2000.0;
    for(i=0;i<=2000;i++)//double合起来算的
    {
        x2=qwe2+i*ewq2/n2;
        v2=valued1(x2);
        fprintf(fp21,"%.100f %.100f\n",x2,v2);
    };
    for(i=0;i<=2000;i++)//double分开算的
    {
        x2=qwe2+i*ewq2/n2;
        v2=valued2(x2);
        fprintf(fp22,"%.100f %.100f\n",x2,v2);
    };
    //fghjjfghjkl//
    //sdgsgldsjgg//
    __float128 qwe3=1.0q-0.003q,ewq3=0.006q,n3=2000.0q;
    for(i=0;i<=2000;i++)//quad合起来算的
    {
        x3=qwe3+i*ewq3/n3;
        v3=valueq1(x3);
        fprintf(fp31,"%.100Qf %.100Qf\n",x3,v3);
    };
    for(i=0;i<=2000;i++)//quad分开算的
    {
        x3=qwe3+i*ewq3/n3;
        v3=valueq2(x3);
        fprintf(fp32,"%.100Qf %.100Qf\n",x3,v3);
    };
    fclose(fp11);
    fclose(fp12);
    fclose(fp22);
    fclose(fp21);
    fclose(fp31);
    fclose(fp32);
}
