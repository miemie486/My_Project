#include<stdio.h>
#include<math.h>
#include<time.h>
#define n 15
/*
double f(double x[n])
{
    int i;
    double a=0;
    for(i=0;i<n;i++)
    {
        a=a+x[i]*x[i];
    }
    return a;
}
*/
int main()
{
    int m=5,ii;
    double start,end;
    double x[n],i[n],k=0,k1=0,a;
    start=clock();
    for(i[0]=0;i[0]<m;i[0]++)
    {
        x[0]=(1.0*i[0])/m+0.5/m;
        for(i[1]=0;i[1]<m;i[1]++)
        {
            x[1]=(1.0*i[1])/m+0.5/m;
            for(i[2]=0;i[2]<m;i[2]++)
            {
                x[2]=(1.0*i[2])/m+0.5/m;
                for(i[3]=0;i[3]<m;i[3]++)
                {
                    x[3]=(1.0*i[3])/m+0.5/m;
                    for(i[4]=0;i[4]<m;i[4]++)
                    {
                        x[4]=(1.0*i[4])/m+0.5/m;
                        for(i[5]=0;i[5]<m;i[5]++)
                        {
                            x[5]=(1.0*i[5])/m+0.5/m;
                            for(i[6]=0;i[6]<m;i[6]++)
                            {
                                x[6]=(1.0*i[6])/m+0.5/m;
                                for(i[7]=0;i[7]<m;i[7]++)
                                {
                                    x[7]=(1.0*i[7])/m+0.5/m;
                                    for(i[8]=0;i[8]<m;i[8]++)
                                    {
                                        x[8]=(1.0*i[8])/m+0.5/m;
                                        for(i[9]=0;i[9]<m;i[9]++)
                                        {
                                            x[9]=(1.0*i[9])/m+0.5/m;
                                            for(i[10]=0;i[10]<m;i[10]++)
                                            {
                                                x[10]=(1.0*i[10])/m+0.5/m;
                                                for(i[11]=0;i[11]<m;i[11]++)
                                                {
                                                    x[11]=(1.0*i[11])/m+0.5/m;
                                                    for(i[12]=0;i[12]<m;i[12]++)
                                                    {
                                                        x[12]=(1.0*i[12])/m+0.5/m;
                                                        for(i[13]=0;i[13]<m;i[13]++)
                                                        {
                                                            x[13]=(1.0*i[13])/m+0.5/m;
                                                            for(i[14]=0;i[14]<m;i[14]++)
                                                            {
                                                                x[14]=(1.0*i[14])/m+0.5/m;
                                                                k1=k1+1;
                                                                a=0;
                                                                for(ii=0;ii<n;ii++)
                                                                {
                                                                    a=a+x[ii]*x[ii];
                                                                }
                                                                if(a<=1)
                                                                {
                                                                    k=k+1;
                                                                }
                                                            }            
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("%.16f\n",k*pow(2,n)/pow(m,n));
    printf("%.0f %.0f\n",k,k1);
    end=clock();
    printf("%f\n",(end-start)/(CLOCKS_PER_SEC));
}