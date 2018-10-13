#include "lowfunc_grad.h"
#include "lowfunc.h"

fgpair_t f_eq42_n0_grad(double r, double a, double c0){
  double f = exp(-a*r)*c0;
  return (fgpair_t){ f, -a*f };
}

fgpair_t g_eq43_l0_grad(double r, double q, double a, double c0){
  double ex  = exp(-a*r);
  double r1 = 1.0/r;
  double u  = 1.0 - ex; // eq44
  double t1 = q * r1 * u;
  double t2 = a * ex * c0;
  double f  = t1+t2;
  double g  = a*(r1*q*ex-t2) - r1*t1;
  return (fgpair_t){ f, g };
}

fgpair_t g_eq43_c0_grad(int l, double r, double q, double a){
  double r1    = 1.0/r;
  double r1_l1 = pow_pos(l+1, r1);
  double l1r1  = (l+1)*r1;
  double ar = a * r;
  double t = 1.0;
  double s = t;
  for(int m=1; m<=2*l; m++){
    t *= (ar/m);    // part of u_eq44_grad
    s += t;         // part of u_eq44
  }
  double ex = exp(-ar);
  double qr = q * r1_l1;
  double u  = 1.0 - s*ex;
  double f  = qr * u;
  double g  = qr * a * t * ex - l1r1 * f;
  return (fgpair_t){ f, g };
}

double g6_eq45_grad(double r, double a, double c){
  double r1 = 1.0/r;
  double r2 = r1*r1;
  double r4 = r2*r2;
  double r6 = r4*r2;
  double ar = a * r;
  double t = 1.0;
  double s = t;
  for(int m=1; m<=9; m++){
    t *= (ar/m);    // part of u_eq44_grad
    s += t;         // part of u_eq44
  }
  return c*r6*((6.0*r1*s+a*t)*exp(-ar) - 6.0*r1);
}

fgpair_t fo_eq46_grad(double r, double a, double b, double c){
  double t1  = exp(-a*r+b);
  double t2  = t1*t1;
  double t21 = 1.0/(t2+1.0);
  double f   = c*t1*t21;
  double g   = f * a * t21 * (t2-1.0);
  return (fgpair_t){ f, g };
}

