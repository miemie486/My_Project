#include <stdio.h>
#include <math.h>
#include <functional>
#include <algorithm>

namespace mathFunc{
  //double rootFinder(std::function<double(double)> func, double a, double b, double error, std::string str);
  double rootFinder(double (*f) (double), double a, double b, double error);
  double max(double *a, int len);
  double min(double *a, int len);
}
