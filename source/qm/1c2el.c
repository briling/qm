#include "eq.h"
#include "common.h"
#include "vec3.h"
#include "matrix.h"

double R0_eq34_old(int mu, int mv, int mu1, int mv1, int lu, int lv, int lu1, int lv1, int ku, mol * m, qmdata * qmd){
  double R00 = R00_eq35(mu,mv,mu1,mv1,lu,lv,lu1,lv1,m->q[ku],qmd);
  double R0k = 0.0;
  if( (mu==mv) && (mu1==mv1) && (lu==lv) && (lu1==lv1) ){
    for(int k=0; k<m->n; k++){
      if(k==ku){
        continue;
      }
      // eq36:
      double r = sqrt(r3d2(m->r+3*ku, m->r+3*k));
      R0k += V_eq51(lu, lu1, m->q[ku], m->q[k], r, qmd);
    }
  }
  return R00+R0k;
}

double R00_eq35(int mu, int mv, int mu1, int mv1, int lu, int lv, int lu1, int lv1, int q, qmdata * qmd){
  int t;
  if(lu  > lv ){
    SWAP(lu,  lv,  t);
    SWAP(mu,  mv,  t);
  }
  if(lu1 > lv1){
    SWAP(lu1, lv1, t);
    SWAP(mu1, mv1, t);
  }
  if(lu  > lu1){
    SWAP(lu,  lu1, t);
    SWAP(lv,  lv1, t);
    SWAP(mu,  mu1, t);
    SWAP(mv,  mv1, t);
  }
  int bra = qmd->q_list[q-1].ga;
  int ket = qmd->q_list[q  ].ga;
  double r = 0.0;
  for(int i=bra; i<ket; i++){
    if( (lu  != qmd->ga[i].lu ) ||
        (lv  != qmd->ga[i].lv ) ||
        (lu1 != qmd->ga[i].lu1) ||
        (lv1 != qmd->ga[i].lv1) ){
      continue;
    }
    int l = qmd->ga[i].l;
    double BB = 0.0;
    double G  = qmd->ga[i].v;
    for(int m=-l; m<=l; m++){
      double B0 = B(l,lu,lv, m,mu,mv);
      double B1 = B(l,lu1,lv1,m,mu1,mv1);
      BB += B0*B1;
    }
    r += BB*G;
  }
  return r;
}

double R0_eq34(int u, int v, int u1, int v1, int ku, double * rij, basis * bo, mol * m, qmdata * qmd){
  int mu  = bo->m[u ];
  int lu  = bo->l[u ];
  int mv  = bo->m[v ];
  int lv  = bo->l[v ];
  int mu1 = bo->m[u1];
  int lu1 = bo->l[u1];
  int mv1 = bo->m[v1];
  int lv1 = bo->l[v1];
  int qu  = bo->Q[u ];
  double R00 = R00_eq35(mu,mv,mu1,mv1,lu,lv,lu1,lv1,qu,qmd);
  if(u!=v || u1!=v1){
    return R00;
  }
  double R0k = 0.0;
  for(int k=0; k<m->n; k++){
    if(k==ku){
      continue;
    }
    // eq36:
    double r = rij[MPOSIF(ku,k)];
    R0k += V_eq51(lu, lu1, qu, m->q[k], r, qmd);
  }
  return R00+R0k;
}

