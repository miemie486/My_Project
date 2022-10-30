#include<stdio.h>
#include<math.h>
#include<quadmath.h>
int main()
{
    float r,b,i,x1,x2;
    for(i=2;i<=8;i++)
    {
        b=pow(10,i);
        printf("x1=%.20f\n",(b+sqrt(pow(b,2)-4))/2);
        printf("x2=%.20f\n",(b-sqrt(pow(b,2)-4))/2);
    }
    printf("\n\n");
    for(i=2;i<=8;i++)
    {
        b=pow(10,i);
        printf("x2'=%.20f\n",2/(b+sqrt(pow(b,2)-4)));
    }
}
