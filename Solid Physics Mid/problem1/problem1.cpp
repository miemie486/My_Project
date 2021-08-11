#include<iostream>
#include<cmath>
#include<vector>
#include <algorithm>

#define N 800
double h=1.0/(2*N);
double f=1,M=1,a=1,S=sqrt(3)*pow(M_PI,2)/2;
double omega(double x, double y)
{
    double p0;
    p0=(2.0*f/M)*(3.0-cos((1.0/sqrt(3.0))*y*a-x*a)-cos(a*2.0*y/sqrt(3.0))-cos((-1.0/sqrt(3.0))*y*a-x*a));
    return sqrt(p0);
}

double del_omega(double x,double y)
{
    double p0,p1;
    p0=(omega(x+h,y)-omega(x-h,y))/(2*h);
    p1=(omega(x,y+h)-omega(x,y-h))/(2*h);
    return sqrt(pow(p0,2)+pow(p1,2));
}
int main()
{
    std::vector<double> density_of_states_omega;
    std::vector<double> density_of_states_y;
    for(double p=1;p<=N;p++)
    {
        std::vector<double> equal_line_up;
        std::vector<double> equal_line_down;
        std::vector<double> equal_line_x_up;
        std::vector<double> equal_line_x;
        std::vector<double> equal_line_y;
        // double omega0=p/N*sqrt(4*f/M*(1-cos(M_PI*a)));
        double omega0=p/N*3;
        for(int i=0;i<=N;i++)
        {
            
            double x=-M_PI+2*M_PI*i/N;
            if(x<-M_PI/2)
            {
               
                double y_up0,y_down0;
                y_up0=sqrt(3)*(x+M_PI);
                y_down0=-sqrt(3)*(x+M_PI);
                double lenth=(y_up0-y_down0)/N;
               
                for(int j=0;j<=N;j++)
                {
                    double number;
                    number=(omega(x,y_up0-j*lenth)-omega0)*(omega(x,y_up0-(j+1)*lenth)-omega0);
                    if(number<=0)
                    {
                        equal_line_up.push_back(fabs(y_up0-(j+0.5)*lenth));
                        equal_line_x_up.push_back(x);
                        break;
                    }
                    
                }
                for(int j=0;j<=N;j++)
                {
                    double number;
                    number=(omega(x,y_down0+j*lenth)-omega0)*(omega(x,y_down0+(j+1)*lenth)-omega0);
                    if(number<=0)
                    {
                        equal_line_down.push_back(y_down0+(j+0.5)*lenth);
                        break;
                    }
                }
            }
            if(x>=-M_PI/2&&x<=M_PI/2)
            {
               
                double y_up0,y_down0;
                y_up0=sqrt(3)*M_PI/2;
                y_down0=-sqrt(3)*M_PI/2;
                double lenth=(y_up0-y_down0)/N;

                for(int j=0;j<=N;j++)
                {
                    double number;
                    number=(omega(x,y_up0-j*lenth)-omega0)*(omega(x,y_up0-(j+1)*lenth)-omega0);
                    if(number<=0)
                    {
                        equal_line_up.push_back(fabs(y_up0-(j+0.5)*lenth));
                        equal_line_x_up.push_back(x);
                        break;
                    }
                    
                }
                for(int j=0;j<=N;j++)
                {
                    double number;
                    number=(omega(x,y_down0+j*lenth)-omega0)*(omega(x,y_down0+(j+1)*lenth)-omega0);
                    if(number<=0)
                    {
                        equal_line_down.push_back(y_down0+(j+0.5)*lenth);
                        break;
                    }
                }
            }
            if(x>M_PI/2)
            {
               
                double y_up0,y_down0;
                y_up0=-sqrt(3)*(x-M_PI);
                y_down0=sqrt(3)*(x-M_PI);
                double lenth=(y_up0-y_down0)/N;

                for(int j=0;j<=N;j++)
                {
                    double number;
                    number=(omega(x,y_up0-j*lenth)-omega0)*(omega(x,y_up0-(j+1)*lenth)-omega0);
                    if(number<=0)
                    {
                        equal_line_up.push_back(fabs(y_up0-(j+0.5)*lenth));
                        equal_line_x_up.push_back(x);
                        break;
                    }
                    
                }
                for(int j=0;j<=N;j++)
                {
                    double number;
                    number=(omega(x,y_down0+j*lenth)-omega0)*(omega(x,y_down0+(j+1)*lenth)-omega0);
                    if(number<=0)
                    {
                        equal_line_down.push_back(y_down0+(j+0.5)*lenth);
                        break;
                    }
                }
            } 
        }
        std::reverse(equal_line_down.begin(),equal_line_down.end());
        equal_line_y.insert(equal_line_y.end(),equal_line_up.begin(),equal_line_up.end());
        equal_line_y.insert(equal_line_y.end(),equal_line_down.begin(),equal_line_down.end());
        equal_line_x.insert(equal_line_x.end(),equal_line_x_up.begin(),equal_line_x_up.end());
        std::reverse(equal_line_x_up.begin(),equal_line_x_up.end());
        equal_line_x.insert(equal_line_x.end(),equal_line_x_up.begin(),equal_line_x_up.end());
        double y_order=0,x_order=0;
        double g_omega=0;
        for(auto i=equal_line_x.begin();i!=equal_line_x.end();i++)
        {
            double dl,del_omega_;
            dl=sqrt(pow(equal_line_x[x_order+1]-equal_line_x[x_order],2)+pow(equal_line_y[y_order+1]-equal_line_y[y_order],2));
            del_omega_=del_omega(equal_line_x[x_order],equal_line_y[y_order]);
            if(del_omega_==0)
            {
                g_omega=g_omega;
            }
            else
            {
                g_omega=g_omega+dl/del_omega_;
            }
            y_order=y_order+1;
            x_order=x_order+1;
        }
        density_of_states_y.push_back(S/pow(2*M_PI,2)*g_omega);
        density_of_states_omega.push_back(omega0);
        std::cout<<omega0<<" "<<g_omega/pow(2*M_PI,2)<<std::endl;
    }
}