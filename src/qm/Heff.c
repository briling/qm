#include "eq.h"
#include "matrix.h"
#include "2el.h"

void Heff(double * Da, double * Db,
    double * Xa,   double * Xb,
    double * FaXa, double * FbXb,
    double * sa,   double * sb,
    double * Fa,   double * Fb,
    double * F2a,  double * F2b,
    double * Fmpa, double * Fmpb,
    double * FA,   double * FB,
    int    * alo,  int    * alv, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd){

  int Mo = bo->M;
  int Mv = bv->M;
  int MM = Mo*Mv;
  double * dEdFa = malloc(sizeof(double)*MM*2);
  double * dEdFb = dEdFa + MM;

  for(int u1=0; u1<Mo; u1++){
    int lu1 = bo->l[u1];
    int ku1 = bo->k[u1];
    double eu1 = qmd->fa [bo->Q[u1]*qmd->nLo+lu1];
    double cu1 = 2.0/(1.0+2.0*lu1);
    for(int a1=0; a1<Mv; a1++){
      int a1u1 = a1*Mo+u1;
      double ea1 = qmd->f2a[bv->Q[a1]*qmd->nLv+bv->l[a1]];
      double eu1a1 = 1.0/(eu1-ea1);
      double sa_u1a1 = sa[u1]*eu1a1;
      double sb_u1a1 = sb[u1]*eu1a1;

      double star1a  = 0.0;
      double star1b  = 0.0;
      double star2aX = 0.0;
      double star2bX = 0.0;
      double star2aF = 0.0;
      double star2bF = 0.0;

      for(int v=0; v<Mo; v++){
        int a1v  = a1*Mo+v;
        int i    = MPOSIF(u1,v);
        star2aX += Da[i] * Xa[a1v] * sa[v];
        star2bX += Db[i] * Xb[a1v] * sb[v];
        star2aF += Da[i] * Fmpa[a1v];
        star2bF += Db[i] * Fmpb[a1v];
        for(int u=alo[ku1]+lu1*lu1 ; u<alo[ku1]+(lu1+1)*(lu1+1); u++){ /* functions on atom ku1 with L==lu1 */
          int vu  = v*Mo+u;
          int vus = MPOSIF(u,v);
          star1a += Da[vus] * FaXa[vu];
          star1b += Db[vus] * FbXb[vu];
        }
      }

      dEdFa[a1*Mo+u1] = star2aX + ( star2aF - star1a * cu1*sa[u1]*Xa[a1u1] ) * sa_u1a1 ;
      dEdFb[a1*Mo+u1] = star2bX + ( star2bF - star1b * cu1*sb[u1]*Xb[a1u1] ) * sb_u1a1 ;
    }
  }

  for(int u1=0; u1<Mo; u1++){
    int ku1 = bo->k[u1];

    for(int v1=u1; v1<alo[ku1+1]; v1++){
      double d = 0.0;
      for(int a=0; a<Mv; a++){
        int ka = bv->k[a];
        if(ka==ku1){
          continue;
        }
        for(int u=alo[ka]; u<alo[ka+1]; u++){
          double R = R2(a,u,u1,v1,pmmm,bo,bv,m,qmd);
          int au = a*Mo+u;
          d += R* ( dEdFa[au] + dEdFb[au] );
        }
      }
      int i = mpos(u1,v1);
      FA[i] = Fa[i] + F2a[i] + d;
      FB[i] = Fb[i] + F2b[i] + d;
    }

    for(int kv1=ku1+1; kv1<m->n; kv1++){
      for(int v1=alo[kv1]; v1<alo[kv1+1]; v1++){
        double dA = 0.0;
        double dB = 0.0;
        for(int a=alv[ku1]; a<alv[ku1+1]; a++){
          for(int u=alo[kv1]; u<alo[kv1+1]; u++){
            double R = R2(a,u1,v1,u,pmmm,bo,bv,m,qmd);
            int au = a*Mo+u;
            dA -= R * dEdFa[au];
            dB -= R * dEdFb[au];
          }
        }
        for(int a=alv[kv1]; a<alv[kv1+1]; a++){
          for(int u=alo[ku1]; u<alo[ku1+1]; u++){
            double R = R2(a,v1,u1,u,pmmm,bo,bv,m,qmd);
            int au = a*Mo+u;
            dA -= R * dEdFa[au] ;
            dB -= R * dEdFb[au] ;
          }
        }
        int i = mpos(u1,v1);
        FA[i] = Fa[i] + F2a[i] + 0.5*dA;
        FB[i] = Fb[i] + F2b[i] + 0.5*dB;
      }
    }

  }

  free(dEdFa);
  return;
}

void Heff_test(int Na, int Nb,
    double * Ca, double * Cb, double * H, double * Hmp,
    int * alo, int * alv, double * mmmm, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd){
  int Mo = bo->M;
  int Mv = bv->M;
  double * Da   = malloc(sizeof(double)*symsize(Mo));
  double * Db   = malloc(sizeof(double)*symsize(Mo));
  double * Fa   = malloc(sizeof(double)*symsize(Mo));
  double * Fb   = malloc(sizeof(double)*symsize(Mo));
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
  double * FA   = malloc(sizeof(double)*symsize(Mo));
  double * FB   = malloc(sizeof(double)*symsize(Mo));
  double E1, E2, e1, e2, de;
  D_eq9  (Na, Mo, Ca, Da);
  D_eq9  (Nb, Mo, Cb, Db);
  F_eq4  (Da,Db, H, Fa,Fb, alo, mmmm, bo,m,qmd);
  F2_8_7_14_15_6(Da, Db, Hmp, Fmpa, Fmpb, Xa, Xb, FaXa, FbXb, sa, sb, F2a, F2b, alo, alv, pmmm, bo, bv, m, qmd);
  Heff   (Da, Db, Xa, Xb, FaXa, FbXb, sa, sb, Fa, Fb, F2a, F2b, Fmpa, Fmpb, FA, FB, alo, alv, pmmm, bo, bv, m, qmd);
  mx_nosym_print(Mo, FA, stdout);
  for(int u=0; u<Mo; u++){
    for(int v=u; v<Mo; v++){
      int i = mpos(u,v);
      double d = 1e-4;
      double D = Da[i];
      Da[i] = D+d;
      F_eq4 (Da,Db, H, Fa,Fb, alo, mmmm, bo,m,qmd);
      F2_8_7_14_15_6(Da, Db, Hmp, Fmpa, Fmpb, Xa, Xb, FaXa, FbXb, sa, sb, F2a, F2b, alo, alv, pmmm, bo, bv, m, qmd);
      E1 = E1_eq3(Mo, H, Da, Db, Fa, Fb);
      E2 = E2_eq5(Mo, Da, Db, F2a, F2b);
      e1 = E1+E2;
      Da[i] = D-d;
      F_eq4(Da,Db, H, Fa,Fb, alo, mmmm, bo,m,qmd);
      F2_8_7_14_15_6(Da, Db, Hmp, Fmpa, Fmpb, Xa, Xb, FaXa, FbXb, sa, sb, F2a, F2b, alo, alv, pmmm, bo, bv, m, qmd);
      E1 = E1_eq3(Mo, H, Da, Db, Fa, Fb);
      E2 = E2_eq5(Mo, Da, Db, F2a, F2b);
      e2 = E1+E2;
      Da[i] = D;
      de = (e1-e2)/d*0.5;
      if(u!=v){
        de *= 0.5;
      }
      printf("%2d,%2d:   0=  % lf  a=  % lf  n=  % lf  a-n=  % .1e  n-0=  % .1e  a-0=  % .1e\n",
          u, v, Fa[i]+F2a[i], FA[i], de, FA[i]-de, de-Fa[i]-F2a[i], FA[i]-Fa[i]-F2a[i]);
    }
  }
  free(Da);
  free(Db);
  free(Fa);
  free(Fb);
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
  free(FA);
  free(FB);
}

