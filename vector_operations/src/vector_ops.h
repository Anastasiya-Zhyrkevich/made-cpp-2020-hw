#pragma once
#include <vector>
#include <iostream>
#include <cassert>


namespace task {

std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2)
{
  assert(v1.size() == v2.size());
  std::vector<double> res(v1.size());

  for (size_t ind = 0; ind < v1.size(); ind++) {
    res[ind] = v1[ind] + v2[ind];
  } 
  return res;
}

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2)
{
  assert(v1.size() == v2.size());
  std::vector<double> res(v1.size());

  for (size_t ind = 0; ind < v1.size(); ind++) {
    res[ind] = v1[ind] - v2[ind];
  } 
  return res;
}

std::vector<double> operator+(const std::vector<double>& v) 
{
  return v;
}

std::vector<double> operator-(const std::vector<double>& v)
{
  std::vector<double> res(v);
  for (size_t ind = 0; ind < v.size(); ind++) {
    res[ind] *= -1.;
  } 
  return res; 
}

double operator*(const std::vector<double>& v1, const std::vector<double>& v2)
{
  assert(v1.size() == v2.size());
  double res = 0.0;

  for (size_t ind = 0; ind < v1.size(); ind++) {
    res += v1[ind] * v2[ind];
  } 
  return res;
}

std::vector<double> operator%(const std::vector<double>& v1, const std::vector<double>& v2)
{
  assert(v1.size() == v2.size());

  if (v1.size() == 2) {
    return std::vector<double>(v1.size(), 0.0);
  }

  assert(v1.size() == v2.size());
  std::vector<double> res(v1.size());

  for (size_t ind = 0; ind < v1.size(); ind++) {
    size_t i1 = (ind + 1) % 3;
    size_t i2 = (ind + 2) % 3;

    res[ind] = v1[i1] * v2[i2] - v1[i2] * v2[i1];
  } 
  return res;
}

bool operator||(const std::vector<double>& v1, const std::vector<double>& v2)
{
  const double EPS = 1e-9;
  assert(v1.size() == v2.size());

  for (size_t ind = 0; ind < v1.size(); ind++) {
    if (fabs(v1[0] * v2[ind] - v1[ind] * v2[0]) > EPS) {
      return false;
    }
  } 
  return true;
}


bool operator&&(const std::vector<double>& v1, const std::vector<double>& v2)
{
  return ( v1 || v2 ) && (v1 * v2 >= 0.0);
}

std::istream& operator>> (std::istream& is, std::vector<double>& v)
{
  size_t n;
  is >> n;

  v.resize(n);
  for (size_t ind = 0; ind < n; ind++) {
    is >> v[ind];
  }
  
  return is;
}

std::ostream& operator<<(std::ostream& os, const std::vector<double>& v)
{
  for (size_t ind = 0; ind < v.size(); ind++) {
    os << v[ind] << ' ';
  }
  os << std::endl;
  return os;
}

void reverse(std::vector<double>& v)
{
  for (size_t ind = 0; 2 * ind < v.size(); ind++) {
    double tmp = v[ind];
    v[ind] = v[v.size() - ind - 1];
    v[v.size() - ind - 1] = tmp;
  } 
}

std::vector<int> operator|(const std::vector<int>& v1, const std::vector<int>& v2)
{
  assert(v1.size() == v2.size());
  std::vector<int> res(v1.size());

  for (size_t ind = 0; ind < v1.size(); ind++) {
    res[ind] = v1[ind] | v2[ind];
  } 
  return res;
}

std::vector<int> operator&(const std::vector<int>& v1, const std::vector<int>& v2)
{
  assert(v1.size() == v2.size());
  std::vector<int> res(v1.size());

  for (size_t ind = 0; ind < v1.size(); ind++) {
    res[ind] = v1[ind] & v2[ind];
  } 
  return res;
}

}  // namespace task
