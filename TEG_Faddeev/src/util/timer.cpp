// Record the exe time
#include "timer.h"

Timer& Timer::operator= (Timer &other){
  _N = other._N;
  _index = other._index;
  std::copy(other._clockPoint, other._clockPoint + _N, _clockPoint);
  std::copy(other._name, other._name + _N - 1, _name);

  return *this;
}

void Timer::record(std::string str){
  if(_index > _N){
    std::cout << "There is no enough clocks\n";
    abort();
  }
  _clockPoint[_index] = omp_get_wtime();
  if(_index > 0) _name[_index - 1] = str;
  _index += 1;
}

double Timer::recent() const {
  int idx = _index - 1;
  return _clockPoint[idx] - _clockPoint[idx - 1];
}

double Timer::interval(int i) const{
  return _clockPoint[_index-1] - _clockPoint[i];
}

void Timer::print() const {
  if(_index < 2){
    std::cout << "Time records are less than 2" << std::endl;
    abort();
  }
  for(size_t idx = 1; idx < _index; idx++)
    std::cout << _name[idx - 1] << "\t" << _clockPoint[idx] - _clockPoint[idx - 1] << "\n";
}

