#include "2el.h"
#include "eq.h"
#include "vec3.h"

double R_old(int u, int v, int u1, int v1, basis * bo, mol * m, qmdata * qmd){
  int ku  = bo->k[u ];
  int kv  = bo->k[v ];
  int ku1 = bo->k[u1];
  int kv1 = bo->k[v1];
  if( (ku!=kv) || (ku1!=kv1) ){
    return 0.0;
  }
  int mu  = bo->m[u ];
  int lu  = bo->l[u ];
  int mv  = bo->m[v ];
  int lv  = bo->l[v ];
  int mu1 = bo->m[u1];
  int lu1 = bo->l[u1];
  int mv1 = bo->m[v1];
  int lv1 = bo->l[v1];
  // one-center
  if(ku==ku1){
    return R0_eq34_old(mu,mv,mu1,mv1,lu,lv,lu1,lv1,ku,m,qmd);
  }
  // two-center
  else{
    int qu  = m->q[ku ];
    int qu1 = m->q[ku1];
    double ruu1[3];
    r3diff(ruu1, m->r+3*ku, m->r+3*ku1);
    double R_0 = R0_eq39_mmmm_old(mu, mv, mu1, mv1, lu, lv, lu1, lv1, qu, qu1, ruu1, qmd);
    // eq20, eq40:
    if( (u == v) && (u1 == v1) ){
      double r   = sqrt(r3dot(ruu1,ruu1));
      double R_6 = G6_eq53(lu, lu1, qu, qu1, r, qmd);
      return R_0-R_6;
    }
    else{
      return R_0;
    }
  }
}

double R2_old(int a, int v, int u1, int v1, basis * bo, basis * bv, mol * m, qmdata * qmd){
  int ka  = bv->k[a ];
  int kv  = bo->k[v ];
  int ku1 = bo->k[u1];
  int kv1 = bo->k[v1];
  if( (ka!=kv) || (ku1!=kv1) || (ka==ku1)){
    return 0.0;
  }
  int ma  = bv->m[a ];
  int la  = bv->l[a ];
  int mv  = bo->m[v ];
  int lv  = bo->l[v ];
  int mu1 = bo->m[u1];
  int lu1 = bo->l[u1];
  int mv1 = bo->m[v1];
  int lv1 = bo->l[v1];
  int qa  = m->q[ka];
  int qu1 = m->q[ku1];
  double r[3];
  r3diff(r, m->r+3*ka, m->r+3*ku1);
  return R0_eq39_mmmp_old(ma,mv,mu1,mv1,la,lv,lu1,lv1,qa,qu1,r,qmd);
}

