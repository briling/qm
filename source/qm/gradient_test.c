#include "qm.h"
#include "matrix.h"
#include "eq.h"
#include "2el.h"
#include "gradient.h"

double calc_energy(int Na, int Nb, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd){
  int Mo = bo->M;
  int Mv = bv->M;
  double * Ca   = calloc(1,sizeof(double)*Mo*Mo);
  double * Cb   = calloc(1,sizeof(double)*Mo*Mo);
  double * Va   = calloc(1,sizeof(double)*Mo);
  double * Vb   = calloc(1,sizeof(double)*Mo);
  double * Da   = calloc(1,sizeof(double)*symsize(Mo));
  double * Db   = calloc(1,sizeof(double)*symsize(Mo));
  double * Dmp  = calloc(1,sizeof(double)*Mo*Mv);
  double * Fw   = calloc(1,sizeof(double)*symsize(Mo));
  double * f    = calloc(1,sizeof(double)*symsize(Mo));
  double * fmp  = calloc(1,sizeof(double)*Mo*Mv);
  double * H    = calloc(1,sizeof(double)*symsize(Mo));
  double * Hmp  = calloc(1,sizeof(double)*Mo*Mv);
  double * rij  = calloc(1,sizeof(double)*symsize(m->n));
  axis  * z    = calloc(1,sizeof(axis )*(m->n)*(m->n));
  distang(rij, z, m);
  double * mmmm = mmmm0_fill(alo, rij, z, bo, m, qmd);
  double * pmmm = pmmm_fill (alo, alv, z, bo, bv, m, qmd);
  f_eq25_mm(f,   z, alo,      bo,     m, qmd);
  f_eq25_mp(fmp, z, alo, alv, bo, bv, m, qmd);
  H_eq22_mm(f,   H,   alo,      mmmm, bo,     m, qmd);
  H_eq22_mp(fmp, Hmp, alo, alv, pmmm, bo, bv, m, qmd);
  mmmm6_add(alo, mmmm, rij, bo, m, qmd);
  mx_id(Mo, Ca);
  veccp(symsize(Mo), Fw, f);
  jacobi(Fw, Ca, Va, Mo, 1e-15, 20, NULL);
  eigensort(Mo, Va, Ca);
  veccp(Mo*Mo, Cb, Ca);
  double E = scf_diis(Na, Nb, E0_eq2(m, qmd), Ca, Cb, Va, Vb, Da, Db, Dmp, 32, 8, 1e-16, alo, alv, H, Hmp, mmmm, pmmm, bo, bv, m, qmd, NULL);
  free(Fw);
  free(Ca);
  free(Cb);
  free(Va);
  free(Vb);
  free(Da);
  free(Db);
  free(Dmp);
  free(mmmm);
  free(pmmm);
  free(H);
  free(Hmp);
  free(f);
  free(fmp);
  free(z);
  free(rij);
  return E;
}

void gradient_test(int Na, int Nb, double * Da, double * Db, double * Hmp, double * pmmm, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd){
  int N = 3*(m->n);
  double * g = malloc(N*sizeof(double));
  gradient(Da, Db, Hmp, pmmm, g, alo, alv, bo, bv, m, qmd);
  g_print(m->n, g, "", stdout);
#if 1
  double d1 = 1e-4;
  for(int i=0; i<N; i++){
    double ri = m->r[i];
    m->r[i] = ri + d1;
    double E1 = calc_energy(Na, Nb, alo, alv, bo, bv, m, qmd);
    m->r[i] = ri - d1;
    double E2 = calc_energy(Na, Nb, alo, alv, bo, bv, m, qmd);
    m->r[i] = ri;
    double gi = (E1-E2)*0.5/d1;
    printf("g%c(%2d) :  a=%20.15lf   n=%20.15lf   (%20.15lf)\n", 'x'+(i%3), i/3, g[i], gi, (g[i]-gi));
  }
#endif

  free(g);
  return;
}

