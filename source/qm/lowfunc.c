#include "lowfunc.h"

static double u_eq44(int n, double r);

double B(int l1, int l2, int l3, int m1, int m2, int m3){

  if( ((l1+l2+l3)%2) || (l1+l2<l3) || (l1+l3<l2) || (l2+l3<l1) ){
    return 0.0;
  }

  int am1 = abs(m1);
  int am2 = abs(m2);
  int am3 = abs(m3);
  if( (( m1)&&(!m2)&&(!m3)) || ((!m1)&&( m2)&&(!m3)) || ((!m1)&&(!m2)&&( m3)) ){
    return 0.0;
  }
  else if( (!m1)&&( m2)&&( m3) ){
    if( (am2!=am3) || (m2==-m3) ){
      return 0.0;
    }
  }
  else if( ( m1)&&(!m2)&&( m3) ){
    if( (am1!=am3) || (m1==-m3) ){
      return 0.0;
    }
  }
  else if( ( m1)&&( m2)&&(!m3) ){
    if( (am1!=am2) || (m1==-m2) ){
      return 0.0;
    }
  }
  else{
    if( (am1+am2!=am3) && (am1+am3!=am2) && (am2+am3!=am1) ){
      return 0.0;
    }
  }

  int t;
  if (l2 > l1) { SWAP(l1, l2, t); SWAP(m1, m2, t);}
  if (l3 > l2) { SWAP(l2, l3, t); SWAP(m2, m3, t);}
  if (l2 > l1) { SWAP(l1, l2, t); SWAP(m1, m2, t);}
  if (l3 > l2) { SWAP(l2, l3, t); SWAP(m2, m3, t);}

#include "B.h"
  abort();

}

double f_eq42_n0(double r, double a, double c0){
  return exp(-a*r)*c0;
}

static double pow_pos(int n, double x){
  if(n==0){
    return 1.0;
  }
  else if(n==1){
    return x;
  }
  double s = x;
  for(int i=1; i<n; i++){
    s *= x;
  }
  return s;
}

double g_eq43_l0(double r, double q, double a, double c0){
  double ex = exp(-a*r);
  double r1 = 1.0/r;
  double u  = 1.0 - ex; // eq44
  double t1 = q * r1 * u;
  double t2 = a * ex * c0;
  return t1+t2;
}

double g_eq43_c0(int l, double r, double q, double a){
  double r1     = 1.0/r;
  double r1_l1  = pow_pos(l+1, r1);
  return q * r1_l1 * u_eq44(2*l, a*r);
}

static double u_eq44(int n, double r){
  double t = 1.0;
  double s = t;
  for(int m=1; m<=n; m++){
    t *= (r/m);    // TODO 1/m can be stored
    s += t;
  }
  return 1.0 - s*exp(-r);
}

double g6_eq45(double r, double a, double c){
  double r1 = 1.0/r;
  double r2 = r1*r1;
  double r4 = r2*r2;
  double r6 = r4*r2;
  return c*r6*u_eq44(9, a*r);
}

double fo_eq46(double r, double a, double b, double c){
  double t1 = exp(-a*r+b);
  double t2 = t1*t1;
  return c*t1/(1.0+t2);
}

