#ifndef TIMER
#define TIMER

#include <iostream>
#include <algorithm>
#include <omp.h>
#include <string>

class Timer{
 private:
  size_t _N, _index = 0;
  double *_clockPoint;
  std::string *_name;
 public:
 Timer(size_t N): _N(N), _clockPoint(new double[N]), _name(new std::string[N - 1]) {};
  ~Timer() {delete[] _clockPoint; delete[] _name;};
  Timer& operator= (Timer&);

  void record(std::string str = ""); // record current time
  void print() const; // print time using
  double recent() const;
  double interval(int) const;
};

#endif /* TIMER */
