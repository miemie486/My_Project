#include <stdio.h>
#include<iostream>
//lapacke headers
#include "lapacke.h"
#include<vector>

using namespace std;

class asd
{
    private:
        vector<double> T;
    public:
        void set_T()
        {
            for(int i=0;i<10;i++)
            {
                T.push_back(i);
            }
        };
        double get_T()
        {
            double aa;
            aa=T[3];
            return aa;
        };
        double kk=3;

};

class abc:public asd
{
    public:
        double kkk=4;
};

class abcd:public abc
{
    public:
        double kkkk=5;
};
int main()
{
    abcd kkkkk;
    cout<<kkkkk.kk<<endl;
}