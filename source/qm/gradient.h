#ifndef GRADIENT_H
#define GRADIENT_H

#include "qm_t.h"
#include "mol.h"

typedef struct{
  double f;
  double g;
} fgpair_t;

void gradient(double * Da, double * Db, double * H, double * pmmm,
    double * g, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd);
void gradient_r(double * D, double * H, double * pmmm,
    double * g, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd);
void gradient_test(int Na, int Nb, double * Da, double * Db, double * Hmp, double * pmmm, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd);

void E0_eq2_grad(double * g, mol * m, qmdata * qmd);
void E0_ext_grad(double field[3], double * g, double * Da, double * Db, int * alo, mol * m, qmdata * qmd);
void f_eq25_mm_grad(double * g, double * Da, double * Db, int * alo, basis * bo, mol * m, qmdata * qmd);
void f_eq25_mm_grad_r(double * g, double * D, int * alo, basis * bo, mol * m, qmdata * qmd);
void R_mmmm_grad(double * g, double * Da, double * Db, int * alo, basis * bo, mol * m, qmdata * qmd);
void R_mmmm_grad_r(double * g, double * D, int * alo, basis * bo, mol * m, qmdata * qmd);
void mmmm6_grad(double * g, double * Da, double * Db, int * alo, basis * bo, mol * m, qmdata * qmd);
void mmmm6_grad_r(double * g, double * D, int * alo, basis * bo, mol * m, qmdata * qmd);
void E2_grad(double * g, double * Da, double * Db,
    double * Fa, double * Fb, double * Xa, double * Xb,
    double * sa, double * sb, double * ja, double * jb,
    int * alo, int * alv, basis * bo , basis * bv, mol * m, qmdata * qmd);
void E2_grad_r(double * g, double * D, double * F, double * X, double * s, double * j,
    int * alo, int * alv, basis * bo , basis * bv, mol * m, qmdata * qmd);

fgpair_t F_eq47_grad(int fbi,  int lu, int lv, int qu, int qv, double r, qmdata * qmd);
fgpair_t F_eq48_grad(int f1bi, int la, int lv, int qa, int qv, double r, qmdata * qmd);
fgpair_t V_eq49_mm_grad(int ubi,  int lu, int lv, int qu, int qk, double r, qmdata * qmd);
fgpair_t V_eq49_mp_grad(int u1bi, int lu, int lv, int qu, int qk, double r, qmdata * qmd);
fgpair_t S_eq50_grad(int foi, int lu, int ld, int qu, int qd, double r, qmdata * qmd);
double   V_eq51_grad(int lu, int lu1, int qu, int qk, double r, qmdata * qmd);
fgpair_t G_eq52_mmmm_grad(int m, int l, int l1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double r, qmdata * qmd);
fgpair_t G_eq52_mmmp_grad(int m, int l, int l1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double q1, double q2, double r, qmdata * qmd);
double   G6_eq53_grad(int lu, int lu1, int qu, int qu1, double r, qmdata * qmd);

#endif
