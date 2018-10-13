#ifndef MOL_H
#define MOL_H

#include "common.h"

typedef struct {
  styp   * s;
  double * r;
  int    * q;
  double * m;
  int    * k;
  int    * b;
  int    * l;
  int      n;
  int      c;
  int      z;
  int      mult;
} mol;

/*==========================================================================================*/

mol  * mol_read    (FILE * f);
void   mol_print_m (mol  * m,  int bohr, const char s[], FILE * f);
void   mol_print2  (mol  * m,  FILE * f);

void eigensort(int n, double * val, double * vec);
int zmat2cart(int n,  double * mr,  double r[3], int a1, int a2,  int a3, double R, double phi, double theta);

#endif

