#include "eq.h"
#include "gradient.h"
#include "lowfunc_grad.h"
#include "matrix.h"
#define EPS 1e-15

fgpair_t F_eq47_grad(int fbi, int lu, int lv, int qu, int qv, double r, qmdata * qmd){
  double au = qmd->afb[qu*qmd->nLo+lu].a;
  double av = qmd->afb[qv*qmd->nLo+lv].a;
  double cv = qmd->afb[qv*qmd->nLo+lv].c;
  double cu = qmd->afb[qu*qmd->nLo+lu].c;
  double cm = qmd->fb[fbi].c0;
  double am = qmd->fb[fbi].a;
  double c = cu*cv*cm;
  double a = sqrt(au*av)*am;
  return f_eq42_n0_grad(r, a, c);
}

fgpair_t F_eq48_grad(int f1bi, int la, int lv, int qa, int qv, double r, qmdata * qmd){
  double av = qmd->af0b[qv*qmd->nLo+lv].a;
  double cv = qmd->af0b[qv*qmd->nLo+lv].c;
  double aa = qmd->af1b[qa*qmd->nLv+la].a;
  double ca = qmd->af1b[qa*qmd->nLv+la].c;
  double am = qmd->f1b[f1bi].a;
  double cm = qmd->f1b[f1bi].c0;
  double a  = 2.0*aa*av/(aa+av)*am;
  double b  = qmd->f1b[f1bi].b;
  double c  = ca*cv*cm;
  return fo_eq46_grad(r,a,b,c);
}

fgpair_t V_eq49_mm_grad(int ubi, int lu, int lv, int qu, int qk, double r, qmdata * qmd){
  double au  = qmd->aub[qu*qmd->nLo+lu].a;
  double cu  = qmd->aub[qu*qmd->nLo+lu].c;
  double av  = qmd->aub[qu*qmd->nLo+lv].a;
  double cv  = qmd->aub[qu*qmd->nLo+lv].c;
  double ak  = qmd->bub[qk].a;
  double ck  = qmd->bub[qk].c;
  double cl  = qmd->ub[ubi].c0;
  double al  = qmd->ub[ubi].a;
  double aa  = au+av;
  double akk = 2.0*ak;
  double a   = aa*akk/(aa+akk)*al;
  double b   = qmd->ub[ubi].b;
  double c   = 0.5*(cu+cv)*ck*cl;
  return fo_eq46_grad(r,a,b,c);
}

fgpair_t V_eq49_mp_grad(int u1bi, int lu, int lv, int qu, int qk, double r, qmdata * qmd){
  double au  = qmd->au1b[qu*qmd->nLv+lu].a;
  double cu  = qmd->au1b[qu*qmd->nLv+lu].c;
  double av  = qmd->au0b[qu*qmd->nLo+lv].a;
  double cv  = qmd->au0b[qu*qmd->nLo+lv].c;
  double ak  = qmd->bu1b[qk].a;
  double ck  = qmd->bu1b[qk].c;
  double cl  = qmd->u1b[u1bi].c0;
  double al  = qmd->u1b[u1bi].a;
  double aa  = au+av;
  double akk = 2.0*ak;
  double a   = aa*akk/(aa+akk)*al;
  double b   = qmd->u1b[u1bi].b;
  double c   = 0.5*(cu+cv)*ck*cl;
  return fo_eq46_grad(r,a,b,c);
}

fgpair_t S_eq50_grad(int foi, int lu, int ld, int qu, int qd, double r, qmdata * qmd){
  double au = qmd->afo[qu*qmd->nLo+lu].a;
  double cu = qmd->afo[qu*qmd->nLo+lu].c;
  double ad = qmd->aof[qd*qmd->nLp+ld].a;
  double cd = qmd->aof[qd*qmd->nLp+ld].c;
  double cm = qmd->fo[foi].c0;
  double a  = 2.0*(au*ad)/(au+ad);
  double b  = qmd->fo[foi].b;
  double c  = cu*cd*cm;
  return fo_eq46_grad(r,a,b,c);
}

double V_eq51_grad(int lu, int lu1, int qu, int qk, double r, qmdata * qmd){
  if(lu>lu1){
    int t;
    SWAP(lu,lu1,t);
  }
  int qq  = qu<=qk?mpos(qu,qk):mpos(qk,qu);
  int bra = qmd->qq_list[qq-1].vb;
  int ket = qmd->qq_list[qq  ].vb;
  for(int i=bra; i<ket; i++){
    if( (qmd->vb[i].qu  != qu ) ||
        (qmd->vb[i].qv  != qk ) ||
        (qmd->vb[i].lu  != lu ) ||
        (qmd->vb[i].lu1 != lu1) ){
      continue;
    }
    if( fabs(qmd->vb[i].c0) > EPS){
      double au  = qmd->avb[qu*qmd->nLo+lu].a;
      double cu  = qmd->avb[qu*qmd->nLo+lu].c;
      double au1 = qmd->avb[qu*qmd->nLo+lu1].a;
      double cu1 = qmd->avb[qu*qmd->nLo+lu1].c;
      double ck  = qmd->bvb[qk].c;
      double ak  = qmd->bvb[qk].a;
      double cuk = qmd->vb[i].c0;
      double aa  = au+au1;
      double a   = 2.0*aa*ak/(aa+2.0*ak);
      double b   = qmd->vb[i].b;
      double c   = 0.5*(cu+cu1)*ck*cuk;
      return (fo_eq46_grad(r,a,b,c)).g;
    }
  }
  return 0.0;
}

fgpair_t G_eq52_mmmm_grad(int m, int l, int l1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double r, qmdata * qmd){
  double q1,q2,q3;
  if( ! qlll_mm(qu,  lu, lv, l,  &q1, qmd)) GOTOHELL;
  if( ! qlll_mm(qu1, lu1,lv1,l1, &q2, qmd)) GOTOHELL;
  q3 = fundconst(l, l1, m);
  double au  = qmd->agb[qu*qmd->nLo+lu];
  double av  = qmd->agb[qu*qmd->nLo+lv];
  double au1 = qmd->agb[qu1*qmd->nLo+lu1];
  double av1 = qmd->agb[qu1*qmd->nLo+lv1];
  double q   = q1*q2*q3;
  double aa  = au+av;
  double aa1 = au1+av1;
  double a   = aa*aa1/(aa+aa1);
  if( (l==0) && (l1==0)){
    return g_eq43_l0_grad(r, q, a, -0.5);
  }
  else{
    int aph = (l+m)%2?-1:1;
    fgpair_t G = g_eq43_c0_grad(l+l1, r, q, a);
    G.f *= aph, G.g *= aph;
    return G;
  }
}

fgpair_t G_eq52_mmmp_grad(int m, int l, int l1, int lu, int lv, int lu1, int lv1, int qu, int qu1, double q1, double q2, double r, qmdata * qmd){
  double q3 = fundconst(l, l1, m);
  double au  = qmd->ag1b[qu *qmd->nLv+lu ];
  double av  = qmd->ag0b[qu *qmd->nLo+lv ];
  double au1 = qmd->ag0b[qu1*qmd->nLo+lu1];
  double av1 = qmd->ag0b[qu1*qmd->nLo+lv1];
  double q   = q1*q2*q3;
  double aa  = au+av;
  double aa1 = au1+av1;
  double a   = aa*aa1/(aa+aa1);
  int    aph = (l+m)%2?-1:1;
  fgpair_t G = g_eq43_c0_grad(l+l1, r, q, a);
  G.f *= aph, G.g *= aph;
  return G;
}

double G6_eq53_grad(int lu, int lu1, int qu, int qu1, double r, qmdata * qmd){
  double au  = qmd->ag6[qu *qmd->nLo+lu ].a;
  double cu  = qmd->ag6[qu *qmd->nLo+lu ].c;
  double au1 = qmd->ag6[qu1*qmd->nLo+lu1].a;
  double cu1 = qmd->ag6[qu1*qmd->nLo+lu1].c;
  double a = 2.0*au*au1/(au+au1);
  double c = cu*cu1;
  return g6_eq45_grad(r, a, c);
}

