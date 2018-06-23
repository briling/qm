#include <math.h>
#include "mol.h"
#include "qm_t.h"
#include "lowfunc.h"

double A(int l, int m1, int m2, double z[3]);
double A_new(int l, int m1, int m2, euler * eu);
void distang(double * rij, euler * eu, mol * m);

double E0_eq2(mol * m, qmdata * qmd);
double E1_eq3(int Mo, double * H, double * Da, double * Db, double * Fa, double * Fb);
double E2_eq5(int Mo, double * Da, double * Db, double * F2a, double * F2b);

void F_eq4(double * Da, double * Db, double * H, double * Fa, double * Fb, int * alo, double * mmmm, basis * bo, mol * m, qmdata * qmd);
void F2_eq6(int Mo, double * FaXa, double * FbXb, double * sa, double * sb, double * F2a, double * F2b);
void X_eq7(double * Fa, double * Fb, double * Xa, double * Xb, basis * bo, basis * bv, qmdata * qmd);
void FX(int Mo, int Mv, double * Fa, double * Fb, double * Xa, double * Xb, double * FaXa, double * FbXb);
void F_eq8(double * Da, double * Db, double * H, double * Fa, double * Fb, int * alo, int * alv, double * pmmm, basis * bo, basis * bv, mol * m, qmdata * qmd);
void D_eq9(int N, int M, double * C, double * D);
void s_eq15(int Mv, double * X, double * s, int * alo, basis * bo, mol * m, qmdata * qmd);

void H_eq22_mm(double * f, double * H, int * alo, double * mmmm, basis * bo, mol * m, qmdata * qmd);
void H_eq22_mp(double * f, double * H, int * alo, int * alv, double * pmmm, basis * bo, basis * bv, mol * m, qmdata * qmd);
void f_eq25_mm(double * f, euler * z, int * alo, basis * bo, mol * m, qmdata * qmd);
void f_eq25_mp(double * f, euler * z, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd);

double R0_eq34_old(int mu, int mv, int mu1, int mv1, int lu, int lv, int lu1, int lv1, int ku, mol * m, qmdata * qmd);
double R0_eq34 (int u, int v, int u1, int v1, int ku, double * rij, basis * bo, mol * m, qmdata * qmd);
double R00_eq35(int mu, int mv, int mu1, int mv1, int lu, int lv, int lu1, int lv1, int q, qmdata * qmd);

double R0_eq39_mmmm_old(int mu, int mv, int mu1, int mv1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double ruu1[3], qmdata * qmd);
double R0_eq39_mmmp_old(int mu, int mv, int mu1, int mv1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double ruu1[3], qmdata * qmd);
double R0_eq39_mmmm(int u, int v, int u1, int v1, int qu, int qu1, euler * z, basis * bo, qmdata * qmd);
double R0_eq39_mmmp(int u, int v, int u1, int v1, int qu, int qu1, euler * z, basis * bo, basis * bv, qmdata * qmd);

void dip_mm(double d[3], int mu, int mv, int lu, int lv, int qu, qmdata * qmd);
void dip_pm(double d[3], int ma, int mv, int la, int lv, int qa, qmdata * qmd);
int qlll_mm(int qu, int lu, int lv, int l, double * q, qmdata * qmd);
int qlll_pm(int qu, int la, int lv, int l, double * q, qmdata * qmd);
double fundconst(int l, int l1, int m);
void dipole(double * Da, double * Db, double * Dmp, double dip[3], int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd);

double E0_ext (double field[3], mol * m, qmdata * qmd);
void   H_ext  (double field[3], double * F, int * alo, basis * bo, mol * m, qmdata * qmd);
void   Hmp_ext(double field[3], double * F, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd);

double F_eq47     (int fbi,  int lu,  int lv, int qu, int qv, double r, qmdata * qmd);
double F_eq48     (int f1bi, int la,  int lv, int qa, int qv, double r, qmdata * qmd);
double V_eq49_mm  (int ubi,  int lu,  int lv, int qu, int qk, double r, qmdata * qmd);
double V_eq49_mp  (int ubi,  int lu,  int lv, int qu, int qk, double r, qmdata * qmd);
double S_eq50     (int foi,  int lu,  int ld, int qu, int qd, double r, qmdata * qmd);
double V_eq51     (int lu,   int lu1, int qu, int qk,  double r, qmdata * qmd);
double G_eq52_mmmm(int m, int l, int l1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double r, qmdata * qmd);
double G_eq52_mmmp(int m, int l, int l1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double q1, double q2, double r, qmdata * qmd);
double G6_eq53    (int lu,   int lu1, int qu, int qu1, double r, qmdata * qmd);


static inline void F2_8_7_14_15_6(
    double * Da, double * Db,
    double * Hmp, double * Fmpa, double * Fmpb,
    double * Xa, double * Xb,
    double * FaXa, double * FbXb,
    double * sa, double * sb,
    double * F2a, double * F2b,
    int * alo, int * alv, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd){
  int Mo = bo->M;
  int Mv = bv->M;
  F_eq8  (Da, Db, Hmp, Fmpa, Fmpb, alo, alv, pmmm, bo, bv, m, qmd);
  X_eq7  (Fmpa, Fmpb, Xa, Xb, bo, bv, qmd);
  s_eq15 (Mv, Xa, sa, alo, bo, m, qmd);
  s_eq15 (Mv, Xb, sb, alo, bo, m, qmd);
  FX(Mo, Mv, Fmpa, Fmpb, Xa, Xb, FaXa, FbXb);
  F2_eq6 (Mo, FaXa, FbXb, sa, sb, F2a, F2b);
  return;
}

void dEdF(double * Da, double * Db,
    double * Xa,    double * Xb,
    double * FaXa,  double * FbXb,
    double * sa,    double * sb,
    double * Fmpa,  double * Fmpb,
    double * dEdFa, double * dEdFb,
    int    * alo, basis * bo, basis * bv, qmdata * qmd);

void Heff(double * dEdFa, double * dEdFb,
    double * Fa,  double * Fb,
    double * F2a, double * F2b,
    double * FA,  double * FB,
    int * alo, int * alv, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd);

void Heff_test(int Na, int Nb,
    double * Ca, double * Cb, double * H, double * Hmp,
    int * alo, int * alv, double * mmmm, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd);

