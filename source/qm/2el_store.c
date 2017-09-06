#include "2el.h"
#include "eq.h"
#include "matrix.h"

/* mo(mv)  -- max number of functions in min.(pol.) basis for one atom
 * b4      -- block size for a pair of atoms
 * bo->n[u] == u-alo[bo->k[u]] , etc
 */

double R(int u, int v, int u1, int v1, double * mmmm, basis * bo, mol * m, qmdata * qmd){
  int ku = bo->k[u];
  int kv = bo->k[v];
  if(ku!=kv){
    return 0.0;
  }
  int ku1 = bo->k[u1];
  int kv1 = bo->k[v1];
  if(ku1!=kv1){
    return 0.0;
  }
  int mo = qmd->nLo*qmd->nLo;
  int b4 = mo*mo*mo*mo;
  int paind = (ku*m->n+ku1)*b4;
  int nu  = bo->n[u ];
  int nv  = bo->n[v ];
  int nu1 = bo->n[u1];
  int nv1 = bo->n[v1];
  int ind = mpos4(nu,nv,nu1,nv1,mo);
  return mmmm[paind+ind];
}

double R2(int a, int v, int u1, int v1, double * pmmm, basis * bo, basis * bv, mol * m, qmdata * qmd){
  int ka = bv->k[a];
  int kv = bo->k[v];
  if(ka!=kv){
    return 0.0;
  }
  int ku1 = bo->k[u1];
  int kv1 = bo->k[v1];
  if( (ku1!=kv1) || (kv==kv1) ){
    return 0.0;
  }
  int mo   = qmd->nLo*qmd->nLo; // максимальное число орбиталей минимального базиса у атома
  int mv   = qmd->nLv*qmd->nLv; // максимальное число орбиталей поляризационного базиса у атома
  int b4   = mv*mo*mo*mo;       // размер блока интегралов для пары атомов
  int paind = (ka*m->n+ku1)*b4;
  int na  = bv->n[a ];
  int nv  = bo->n[v ];
  int nu1 = bo->n[u1];
  int nv1 = bo->n[v1];
  int ind = mpos4(na,nv,nu1,nv1,mo);
  return pmmm[paind+ind];
}

void mmmm6_add(int * alo, double * mmmm, double * rij, basis * bo, mol * m, qmdata * qmd){
  int mo = qmd->nLo*qmd->nLo;
  int b4 = mo*mo*mo*mo;
  for(int ku=0; ku<m->n; ku++){
    int qu = m->q[ku];
    for(int ku1=ku+1; ku1<m->n; ku1++){
      int qu1    = m->q[ku1];
      int paind0 = (ku *m->n+ku1)*b4;
      int paind1 = (ku1*m->n+ku )*b4;
      for(int u=alo[ku]; u<alo[ku+1]; u++){
        for(int u1=alo[ku1]; u1<alo[ku1+1]; u1++){
          int lu    = bo->l[u ];
          int lu1   = bo->l[u1];
          double r  = rij[MPOSIF(ku, ku1)];
          double r6 = G6_eq53(lu, lu1, qu, qu1, r, qmd);
          int ind0  = mpos4( bo->n[u ], bo->n[u ], bo->n[u1], bo->n[u1], mo);
          int ind1  = mpos4( bo->n[u1], bo->n[u1], bo->n[u ], bo->n[u ], mo);
          mmmm[paind0+ind0] -= r6;
          mmmm[paind1+ind1] -= r6;
        }
      }
    }
  }
  return;
}

double * mmmm0_fill(int * alo, double * rij, euler * z, basis * bo, mol * m, qmdata * qmd){
  int mo   = qmd->nLo*qmd->nLo;
  int b4   = mo*mo*mo*mo;
  int nnb4 = m->n*m->n*b4;
  double * mmmm = malloc(nnb4*sizeof(double));
  for(int k=0; k<m->n; k++){
    int paind = (k*m->n+k)*b4;
    for(int u=alo[k]; u<alo[k+1]; u++){
      for(int v=u; v<alo[k+1]; v++){
        for(int u1=alo[k]; u1<alo[k+1]; u1++){
          for(int v1=u1; v1<alo[k+1]; v1++){
            int nu  = bo->n[u ];
            int nv  = bo->n[v ];
            int nu1 = bo->n[u1];
            int nv1 = bo->n[v1];
            int ind1 = mpos4(nu, nv, nu1, nv1, mo);
            int ind2 = mpos4(nv, nu, nu1, nv1, mo);
            int ind3 = mpos4(nu, nv, nv1, nu1, mo);
            int ind4 = mpos4(nv, nu, nv1, nu1, mo);
            double r = R0_eq34(u,v,u1,v1,k,rij,bo,m,qmd);
            mmmm[paind+ind1] = r;
            mmmm[paind+ind2] = r;
            mmmm[paind+ind3] = r;
            mmmm[paind+ind4] = r;
          }
        }
      }
    }
    for(int k1=k+1; k1<m->n; k1++){
      int paind0 = (k *m->n+k1)*b4;
      int paind1 = (k1*m->n+k )*b4;
      euler * zuu1 = z+(k*m->n+k1);
      for(int u=alo[k]; u<alo[k+1]; u++){
        for(int v=u; v<alo[k+1]; v++){
          for(int u1=alo[k1]; u1<alo[k1+1]; u1++){
            for(int v1=u1; v1<alo[k1+1]; v1++){
              int nu  = bo->n[u ];
              int nv  = bo->n[v ];
              int nu1 = bo->n[u1];
              int nv1 = bo->n[v1];
              int ind000 = mpos4(nu , nv , nu1, nv1, mo);
              int ind010 = mpos4(nu , nv , nv1, nu1, mo);
              int ind001 = mpos4(nv , nu , nu1, nv1, mo);
              int ind011 = mpos4(nv , nu , nv1, nu1, mo);
              int ind100 = mpos4(nu1, nv1, nu , nv , mo);
              int ind110 = mpos4(nv1, nu1, nu , nv , mo);
              int ind101 = mpos4(nu1, nv1, nv , nu , mo);
              int ind111 = mpos4(nv1, nu1, nv , nu , mo);
              double r = R0_eq39_mmmm(u,v,u1,v1, m->q[k],m->q[k1], zuu1, bo,qmd);
              mmmm[paind0+ind000] = r;
              mmmm[paind0+ind010] = r;
              mmmm[paind0+ind001] = r;
              mmmm[paind0+ind011] = r;
              mmmm[paind1+ind100] = r;
              mmmm[paind1+ind110] = r;
              mmmm[paind1+ind101] = r;
              mmmm[paind1+ind111] = r;
            }
          }
        }
      }
    }
  }
  return mmmm;
}

double * pmmm_fill(int * alo, int * alv, euler * z, basis * bo, basis * bv, mol * m, qmdata * qmd){
  int mo   = qmd->nLo*qmd->nLo;
  int mv   = qmd->nLv*qmd->nLv;
  int b4   = mv*mo*mo*mo;
  int nnb4 = m->n*m->n*b4;
  double * pmmm = malloc(nnb4*sizeof(double));
  for(int k=0; k<m->n; k++){
    for(int k1=0; k1<m->n; k1++){
      if(k==k1){
        continue;
      }
      int paind = (k*m->n+k1)*b4;
      euler * zuu1 = z+(k*m->n+k1);
      for(int a=alv[k]; a<alv[k+1]; a++){
        for(int v=alo[k]; v<alo[k+1]; v++){
          for(int u1=alo[k1]; u1<alo[k1+1]; u1++){
            for(int v1=u1; v1<alo[k1+1]; v1++){
              int na   = bv->n[a ];
              int nv   = bo->n[v ];
              int nu1  = bo->n[u1];
              int nv1  = bo->n[v1];
              int ind0 = mpos4(na, nv, nu1, nv1, mo);
              int ind1 = mpos4(na, nv, nv1, nu1, mo);
              double r = R0_eq39_mmmp(a, v, u1, v1, m->q[k], m->q[k1], zuu1, bo, bv, qmd);
              pmmm[paind+ind0] = r;
              pmmm[paind+ind1] = r;
            }
          }
        }
      }
    }
  }
  return pmmm;
}

void mmmm_check(double * mmmm, basis * bo, mol * m, qmdata * qmd){
  int Mo = bo->M;
  for(int u=0; u<Mo; u++){
    for(int v=0; v<Mo; v++){
      for(int u1=0; u1<Mo; u1++){
        for(int v1=0; v1<Mo; v1++){
          double r1 = R_old(u,v,u1,v1, bo,m,qmd);
          double r2 = R    (u,v,u1,v1, mmmm, bo,m,qmd);
          printf("%3d %3d %3d %3d    % e   % e\n", u,v,u1,v1, r1,r2);
          if(fabs(r1-r2)>1e-15){
            abort();
          }
        }
      }
    }
  }
  return;
}

void pmmm_check(double * pmmm, basis * bo, basis * bv, mol * m, qmdata * qmd){
  int Mo = bo->M;
  int Mv = bv->M;
  for(int a=0; a<Mv; a++){
    for(int v=0; v<Mo; v++){
      for(int u1=0; u1<Mo; u1++){
        for(int v1=0; v1<Mo; v1++){
          double r1 = R2_old(a,v,u1,v1, bo,bv,m,qmd);
          double r2 = R2    (a,v,u1,v1, pmmm, bo,bv,m,qmd);
          printf("%3d %3d %3d %3d    % e   % e\n", a,v,u1,v1, r1,r2);
          if(fabs(r1-r2)>1e-15){
            abort();
          }
        }
      }
    }
  }
  return;
}

