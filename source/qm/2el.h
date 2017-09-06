#include "mol.h"
#include "qm_t.h"
#define mpos4(i,j,k,l,M)  ((i)*(M)*(M)*(M)+(j)*(M)*(M)+(k)*(M)+(l))

double R_old(int u, int v, int u1, int v1, basis * bo, mol * m, qmdata * qmd);
double R2_old(int a, int v, int u1, int v1, basis * bo, basis * bv, mol * m, qmdata * qmd);

double * mmmm0_fill(int * alo, double * rij, euler * z, basis * bo, mol * m, qmdata * qmd);
void     mmmm6_add (int * alo, double * mmmm, double * rij, basis * bo, mol * m, qmdata * qmd);
double * pmmm_fill (int * alo, int * alv, euler * z, basis * bo, basis * bv, mol * m, qmdata * qmd);

void mmmm_check(double * mmmm, basis * bo, mol * m, qmdata * qmd);
void pmmm_check(double * pmmm, basis * bo, basis * bv, mol * m, qmdata * qmd);

double R (int u, int v, int u1, int v1, double * mmmm, basis * bo, mol * m, qmdata * qmd);
double R2(int a, int v, int u1, int v1, double * pmmm, basis * bo, basis * bv, mol * m, qmdata * qmd);

