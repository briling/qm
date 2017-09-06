#include "3d.h"

void rotmx(double * rot, double * u, double phi){
  double cphi = cos(phi);
  double sphi = sin(phi);
  double R[9] = {cphi, 0.0,  0.0,
                 0.0,  cphi, 0.0,
                 0.0,  0.0,  cphi};
  double ux[9] = {0.0, -u[2], u[1],
                  u[2], 0.0, -u[0],
                 -u[1], u[0], 0.0};
  double uu[9] = {u[0]*u[0], u[0]*u[1], u[0]*u[2],
                  u[0]*u[1], u[1]*u[1], u[1]*u[2],
                  u[0]*u[2], u[1]*u[2], u[2]*u[2]};
  vecadds(9, R, ux, sphi);
  vecadds(9, R, uu, 1.0-cphi);
  veccp  (9, rot, R);
  return;
}

