#include "vecn.h"

void vecset(size_t n, double * u, double s){
  for(size_t i=0; i<n; i++){
    u[i] = s;
  }
  return;
}

void veccp(size_t n, double * u, double * v){
  for(size_t i=0; i<n; i++){
    u[i] = v[i];
  }
  return;
}

void vecsum(size_t n, double * w, double * u, double * v){
  for(size_t i=0; i<n; i++){
    w[i] = u[i]+v[i];
  }
  return;
}

void vecscal(size_t n, double * u, double s){
  for(size_t i=0; i<n; i++){
    u[i] *= s;
  }
  return;
}

void vecadd(size_t n, double * u, double * v){
  for(size_t i=0; i<n; i++){
    u[i] += v[i];
  }
  return;
}

void vecadds(size_t n, double * u, double * v, double s){
  for(size_t i=0; i<n; i++){
    u[i] += v[i]*s;
  }
  return;
}

double vecabsmax(size_t n, double * u){
  double s = fabs(u[0]);
  for(size_t i=1; i<n; i++){
    double t = fabs(u[i]);
    if(s<t){
      s = t;
    }
  }
  return s;
}

