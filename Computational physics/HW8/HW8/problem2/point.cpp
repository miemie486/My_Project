#include<stdio.h>
#include<math.h>
#include<string.h>
#define n 4   //数据文件个数
double ff(double x)
{
    double a;
    a=x*sin(2*M_PI*x+1);
    return a;
}
int main()
{
    char path[n][255],a[255],b[255];
    FILE *f[n];
    int k[n]={9,17,19,21};
    int i,j;
    double x,y;
    for(i=0;i<n;i++)
    {
        strcpy(a,"data");
        sprintf(b, "%d", k[i]);
        strcpy(a,strcat(a,b));
        strcpy(a,strcat(a,".txt"));
        strcpy(path[i],a);
        f[i]=fopen(path[i],"w+");
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<k[i];j++)
        {
            x=-1+2*j/(k[i]-1.0);
            fprintf(f[i],"%f ",x);
        }
        fprintf(f[i],"\n");
        for(j=0;j<k[i];j++)
        {
            y=ff(-1+2*j/(k[i]-1.0));
            fprintf(f[i],"%f ",y);
        }
    }
    fcloseall();
}