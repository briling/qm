#include "qm.h"
#include "vecn.h"
#include "eq.h"
#include "gradient.h"
#include "derivatives.h"

void gradient(double * Da, double * Db, double * H, double * pmmm,
    double * g, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd){

  vecset(3*m->n, g, 0.0);
  E0_eq2_grad   (g, m, qmd);
  mmmm6_grad    (g, Da, Db, alo, bo, m, qmd);
  f_eq25_mm_grad(g, Da, Db, alo, bo, m, qmd);
  R_mmmm_grad   (g, Da, Db, alo, bo, m, qmd);

  int Mo = bo->M;
  int Mv = bv->M;
  double * Fa   = malloc(sizeof(double)*Mo*Mv);
  double * Fb   = malloc(sizeof(double)*Mo*Mv);
  double * Xa   = malloc(sizeof(double)*Mo*Mv);
  double * Xb   = malloc(sizeof(double)*Mo*Mv);
  double * sa   = malloc(sizeof(double)*Mo);
  double * sb   = malloc(sizeof(double)*Mo);
  double * FaXa = malloc(sizeof(double)*Mo*Mo);
  double * FbXb = malloc(sizeof(double)*Mo*Mo);
  double * ja   = malloc(sizeof(double)*Mo);
  double * jb   = malloc(sizeof(double)*Mo);

  F_eq8  (Da, Db, H, Fa, Fb, alo, alv, pmmm, bo, bv, m, qmd);
  X_eq7  (Fa, Fb, Xa, Xb, bo, bv, qmd);
  s_eq15 (Mv, Xa, sa, alo, bo, m, qmd);
  s_eq15 (Mv, Xb, sb, alo, bo, m, qmd);
  FX(Mo, Mv, Fa, Fb, Xa, Xb, FaXa, FbXb);
  d2E2dF2_j(Da, Db, FaXa, FbXb, ja, jb, alo, bo);

  E2_grad(g, Da, Db, Fa, Fb, Xa, Xb, sa, sb, ja, jb, alo, alv, bo, bv, m, qmd);

  free(Fa);
  free(Fb);
  free(Xa);
  free(Xb);
  free(sa);
  free(sb);
  free(FaXa);
  free(FbXb);
  free(ja);
  free(jb);

  return;
}

void gradient_r(double * D, double * H, double * pmmm,
    double * g, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd){

  vecset(3*m->n, g, 0.0);
  E0_eq2_grad     (g, m, qmd);
  mmmm6_grad_r    (g, D, alo, bo, m, qmd);
  f_eq25_mm_grad_r(g, D, alo, bo, m, qmd);
  R_mmmm_grad_r   (g, D, alo, bo, m, qmd);

  int Mo = bo->M;
  int Mv = bv->M;

  double * Fmp = malloc(sizeof(double)*Mo*Mv);
  double * X   = malloc(sizeof(double)*Mo*Mv);
  double * FX  = malloc(sizeof(double)*Mo*Mo);
  double * s   = malloc(sizeof(double)*Mo);
  double * j   = malloc(sizeof(double)*Mo);
  F_eq8_r  (D, H, Fmp, alo, alv, pmmm, bo, bv, m, qmd);
  X_eq7_r  (Fmp, X, bo, bv, qmd);
  s_eq15 (Mv, X, s, alo, bo, m, qmd);
  FX_r (Mo, Mv, Fmp, X, FX);
  d2E2dF2_j_r(D, FX, j, alo, bo);

  E2_grad_r(g, D, Fmp, X, s, j, alo, alv, bo, bv, m, qmd);

  free(Fmp);
  free(X);
  free(s);
  free(FX);
  free(j);
  return;
}

/*---------------------------------------------------------------------------*/

