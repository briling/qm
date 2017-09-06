#include <stdlib.h>
#include <math.h>

#define SWAP(X,Y,T) { T = X; X = Y; Y = T; }
#define SQRT3    1.73205080756887729352
#define SQRTPI   1.77245385090551602729

double B(int l1, int l2, int l3, int m1, int m2, int m3);
double f_eq42_n0(double r, double a, double c0);
double g_eq43_l0(double r, double q, double a, double c0);
double g_eq43_c0(int l, double r, double q, double a);
double g6_eq45(double r, double a, double c);
double fo_eq46(double r, double a, double b, double c);

