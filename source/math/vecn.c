#include "vecn.h"

void vecset(int n, double * r, double s){
  int i;
  for (i=0; i<n; i++){
    r[i] = s;
  }
  return;
}

void veccp(int n, double * u, double * v){
  int i;
  for (i=0; i<n; i++){
    u[i] = v[i];
  }
  return;
}

void vecscal(size_t n, double * u, double s){
  size_t i;
  for (i=0; i<n; i++){
    u[i] *= s;
  }
  return;
}

void vecadds(int n, double * u, double * v, double s){
  int i;
  for (i=0; i<n; i++){
    u[i] += v[i]*s;
  }
  return;
}

