#include "eq.h"
#include "matrix.h"

void Deff(double * Da, double * Db,
    double * Xa,   double * Xb,
    double * FaXa, double * FbXb,
    double * sa,   double * sb,
    double * Fmpa, double * Fmpb,
    double * Dmp,
    int * alo, basis * bo, basis * bv, qmdata * qmd){

  int Mo = bo->M;
  int Mv = bv->M;

  for(int v1=0; v1<Mo; v1++){

    int kv1 = bo->k[v1];
    int lv1 = bo->l[v1];
    int bra = alo[kv1] +  lv1   * lv1;
    int ket = alo[kv1] + (lv1+1)*(lv1+1);
    double cv1 = 2.0/(1.0+2.0*lv1);
    double ssa = sa[v1]*sa[v1]*cv1;
    double ssb = sb[v1]*sb[v1]*cv1;
    double ev1 = qmd->fa [bo->Q[v1]*qmd->nLo+lv1];

    for(int a1=0; a1<Mv; a1++){
      int a1v1 = a1*Mo+v1;
      double ea1 = qmd->f2a[bv->Q[a1]*qmd->nLv+bv->l[a1]];
      double ev1a1 = 1.0/(ev1-ea1);
      double ssxea = ssa * Xa[a1v1] * ev1a1;
      double ssxeb = ssb * Xb[a1v1] * ev1a1;
      double sav1a1 = sa[v1]*ev1a1;
      double sbv1a1 = sb[v1]*ev1a1;

      double D  = 0.0;
      double da = 0.0;
      double db = 0.0;
      for(int u=0; u<Mo; u++){
        int lu  = bo->l[u];
        double eu = qmd->fa [bo->Q[u ]*qmd->nLo+lu ];
        double eua1 = 1.0/(eu-ea1);
        int uv1 = MPOSIF(u,v1);
        int a1u = a1*Mo+u;
        D += Da[uv1] * Fmpa[a1u] * (sav1a1 + sa[u]*eua1);
        D += Db[uv1] * Fmpb[a1u] * (sbv1a1 + sb[u]*eua1);

        for(int v=bra; v<ket; v++){  /* functions on atom kv1 with L==lv1 */
          if(bo->l[v] == bo->l[v1]){
            int i = MPOSIF(u,v);
            int j = u*Mo+v;
            da -= Da[i] * FaXa[j];
            db -= Db[i] * FbXb[j];
          }
        }
      }
      D += ssxea * da + ssxeb * db;

      Dmp[a1v1] = D;
    }
  }

  return;
}

void Deff_test(int Na, int Nb,
    double * Ca,  double * Cb,
    double * Hmp, double * Dmp,
    int * alo, int * alv, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd){
  int Mo = bo->M;
  int Mv = bv->M;
  double * Da   = malloc(sizeof(double)*symsize(Mo));
  double * Db   = malloc(sizeof(double)*symsize(Mo));
  double * Fmpa = malloc(sizeof(double)*Mo*Mv);
  double * Fmpb = malloc(sizeof(double)*Mo*Mv);
  double * Xa   = malloc(sizeof(double)*Mo*Mv);
  double * Xb   = malloc(sizeof(double)*Mo*Mv);
  double * FaXa = malloc(sizeof(double)*Mo*Mo);
  double * FbXb = malloc(sizeof(double)*Mo*Mo);
  double * sa   = malloc(sizeof(double)*Mo);
  double * sb   = malloc(sizeof(double)*Mo);
  double * F2a  = malloc(sizeof(double)*symsize(Mo));
  double * F2b  = malloc(sizeof(double)*symsize(Mo));
  D_eq9(Na, Mo, Ca, Da);
  D_eq9(Nb, Mo, Cb, Db);

  for(int a=0; a<Mv; a++){
    for(int v=0; v<Mo; v++){
      int i = a*Mo+v;
      double d = 1e-4;
      double H = Hmp[i];
      Hmp[i] = H+d;
      F2_8_7_14_15_6(Da, Db, Hmp, Fmpa, Fmpb, Xa, Xb, FaXa, FbXb, sa, sb, F2a, F2b, alo, alv, pmmm, bo, bv, m, qmd);
      double e1 = E2_eq5(Mo, Da, Db, F2a, F2b);
      Hmp[i] = H-d;
      F2_8_7_14_15_6(Da, Db, Hmp, Fmpa, Fmpb, Xa, Xb, FaXa, FbXb, sa, sb, F2a, F2b, alo, alv, pmmm, bo, bv, m, qmd);
      double e2 = E2_eq5(Mo, Da, Db, F2a, F2b);
      Hmp[i] = H;
      double de = (e1-e2)/d*0.5;
      printf("%2d,%2d:   a=  % lf  n=  % lf  a-n=  % .1e\n", a, v, Dmp[i], de, Dmp[i]-de);
    }
  }

  free(Da);
  free(Db);
  free(Fmpa);
  free(Fmpb);
  free(Xa);
  free(Xb);
  free(FaXa);
  free(FbXb);
  free(sa);
  free(sb);
  free(F2a);
  free(F2b);
  return;
}

