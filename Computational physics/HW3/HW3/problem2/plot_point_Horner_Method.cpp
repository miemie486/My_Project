#include<stdio.h>
#include<math.h>
#include<quadmath.h>

float f(float x)
{
    float b;
    b=((((((((x-9)*x+36)*x-84)*x+126)*x-126)*x+84)*x-36)*x+9)*x-1;
    return b;
}
int main()
{
    int i;
    float x1,v1,k;
    FILE *fp;
    fp=fopen("dataHorner.txt","w+");
    float qwe1=1-0.3,ewq1=0.6,n1=2000.0;
    for(i=0;i<=2000;i++)
    {
        x1=qwe1+i*ewq1/n1;
        v1=f(x1);
        fprintf(fp,"%.100f %.100f\n",x1,v1);
    }
    fclose(fp);
}