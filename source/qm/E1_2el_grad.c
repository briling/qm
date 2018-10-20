#include "qm.h"
#include "matrix.h"
#include "eq.h"
#include "2el.h"
#include "gradient.h"
#include "vec3.h"
#define EPS 1e-15

void mmmm6_grad(double * g, double * Da, double * Db, int * alo, basis * bo, mol * m, qmdata * qmd){
  for(int ku=0; ku<m->n; ku++){
    int qu = m->q[ku];
    for(int ku1=ku+1; ku1<m->n; ku1++){
      int qu1 = m->q[ku1];
      for(int u=alo[ku]; u<alo[ku+1]; u++){
        int lu = bo->l[u ];
        int uu = mpos(u,u);
        for(int u1=alo[ku1]; u1<alo[ku1+1]; u1++){
          int lu1  = bo->l[u1];
          int u1u1 = mpos(u1,u1);
          int uu1  = MPOSIF(u,u1);
          double r  = sqrt(r3d2(m->r+3*ku,m->r+3*ku1));
          double r6 = G6_eq53_grad(lu, lu1, qu, qu1, r, qmd);
          double dd = (
              +Da[uu]*Da[u1u1]
              +Db[uu]*Da[u1u1]
              +Da[uu]*Db[u1u1]
              +Db[uu]*Db[u1u1]
              -Da[uu1]*Da[uu1]
              -Db[uu1]*Db[uu1]);
          double dr[3];
          r3diff(dr, m->r+ku*3, m->r+ku1*3);
          r3scal(dr, -dd*r6/r);
          r3add(g+3*ku,  dr);
          r3min(g+3*ku1, dr);
        }
      }
    }
  }
  return;
}

void mmmm6_grad_r(double * g, double * D, int * alo, basis * bo, mol * m, qmdata * qmd){
  for(int ku=0; ku<m->n; ku++){
    int qu = m->q[ku];
    for(int ku1=ku+1; ku1<m->n; ku1++){
      int qu1 = m->q[ku1];
      for(int u=alo[ku]; u<alo[ku+1]; u++){
        int lu = bo->l[u ];
        int uu = mpos(u,u);
        for(int u1=alo[ku1]; u1<alo[ku1+1]; u1++){
          int lu1  = bo->l[u1];
          int u1u1 = mpos(u1,u1);
          int uu1  = MPOSIF(u,u1);
          double r  = sqrt(r3d2(m->r+3*ku,m->r+3*ku1));
          double r6 = G6_eq53_grad(lu, lu1, qu, qu1, r, qmd);
          double dd = 4.0*D[uu]*D[u1u1] - 2.0*D[uu1]*D[uu1];
          double dr[3];
          r3diff(dr, m->r+ku*3, m->r+ku1*3);
          r3scal(dr, -dd*r6/r);
          r3add(g+3*ku,  dr);
          r3min(g+3*ku1, dr);
        }
      }
    }
  }
  return;
}

static void R0_eq39_mmmm_grad(int u, int v, int u1, int v1,
    double g[3], double ru[3], double ru1[3], basis * bo, qmdata * qmd){

  double z[3];
  r3diff(z, ru, ru1);
  double r  = sqrt(r3dot(z,z));
  double r1 = 1.0/r;
  r3scal(z, r1);

  int qu  = bo->Q[u ];
  int qu1 = bo->Q[u1];
  int mu  = bo->m[u ];
  int lu  = bo->l[u ];
  int mv  = bo->m[v ];
  int lv  = bo->l[v ];
  int mu1 = bo->m[u1];
  int lu1 = bo->l[u1];
  int mv1 = bo->m[v1];
  int lv1 = bo->l[v1];

  double A1, A2, A3, A4;
  double gA1[3], gA2[3], gA3[3], gA4[3];
  double g1[3]={0}, g2[3]={0};

  for(int m_=-lu;m_<=lu;m_++){
    A1 = A_grad_z(lu,m_,mu,gA1,z);
    for(int m__=-lv;m__<=lv;m__++){
      A2 = A_grad_z(lv,m__,mv,gA2,z);
      for(int m1_=-lu1;m1_<=lu1;m1_++){
        A3 = A_grad_z(lu1,m1_,mu1,gA3,z);
        for(int m1__=-lv1;m1__<=lv1;m1__++){
          A4 = A_grad_z(lv1,m1__,mv1,gA4,z);

          double AAAA = A1*A2*A3*A4;
          double BBG  = 0.0;
          double BBgG = 0.0;
          double tg[3]={0};
          r3adds(tg, gA1, A2*A3*A4);
          r3adds(tg, gA2, A1*A3*A4);
          r3adds(tg, gA3, A1*A2*A4);
          r3adds(tg, gA4, A1*A2*A3);

          for(int l=abs(lu-lv); l<=lu+lv; l++){
            for(int l1=abs(lu1-lv1); l1<=lu1+lv1; l1++){
              int lm = MIN(l,l1);
              for(int m=-lm; m<=lm; m++){
                double B1 = B(l,  lu,  lv,  m, m_,  m__);
                if(fabs(B1)<EPS){
                  continue;
                }
                double B2 = B(l1, lu1, lv1, m, m1_, m1__);
                if(fabs(B2)<EPS){
                  continue;
                }
                fgpair_t G = G_eq52_mmmm_grad(m,l,l1,lu,lv,lu1,lv1,qu,qu1,r,qmd);
                BBgG += B1*B2*G.g;
                BBG  += B1*B2*G.f;
              }
            }
          }
          r3adds(g1, z, AAAA*BBgG);
          r3adds(g2, tg, BBG);
        }
      }
    }
  }
  A_grad_z2r(g2, z, r1);
  r3sum(g, g1, g2);
  return;
}

void R_mmmm_grad(double * g, double * Da, double * Db, int * alo, basis * bo, mol * m, qmdata * qmd){
  // two-center corrections to one-center
  for(int ku=0; ku<m->n; ku++){
    int qu = m->q[ku];
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      int uu = mpos(u,u);
      int lu = bo->l[u];
      double pu = qmd->p[qu*(qmd->nLo)+lu ]/(2.0*lu +1.0);
      double Du = Da[uu] + Db[uu];
      for(int u1=u; u1<alo[ku+1]; u1++){
        int lu1 = bo->l[u1];
        double pu1 = qmd->p[qu*(qmd->nLo)+lu1]/(2.0*lu1+1.0);
        int u1u1 = mpos(u1,u1);
        int uu1  = MPOSIF(u,u1);
        double D1  = Da[uu] * Da[u1u1] + Db[uu] * Db[u1u1];
        double D2  = Da[uu] * Db[u1u1] + Db[uu] * Da[u1u1];
        double D3  = Da[uu1] * Da[uu1] + Db[uu1] * Db[uu1];
        double Du1 = Da[u1u1] + Db[u1u1];
        for(int k1=0; k1<m->n; k1++){
          if(k1==ku) continue;
          double dr[3];
          r3diff(dr, m->r+3*ku, m->r+3*k1);
          double r  = sqrt(r3dot(dr,dr));
          double gV = V_eq51_grad(lu, lu1, qu, m->q[k1], r, qmd);
          double fact1 = (u==u1?0.5:1.0) * (D1+D2-D3);
          double fact2 = (u==u1 ? 0.5*Du*pu1 : Du*pu1+Du1*pu);  // eq22 (R -> H -> E)
          r3scal(dr, (fact1-fact2) * gV/r);
          r3add(g+3*ku, dr);
          r3min(g+3*k1, dr);
        }
      }
    }
  }
  // two-center
  for(int ku=0; ku<m->n; ku++){
    int qu = m->q[ku];
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      int lu = bo->l[u];
      double pu  = qmd->p[qu *(qmd->nLo)+lu ]/(2.0*lu +1.0);
      for(int ku1=ku+1; ku1<m->n; ku1++){
        for(int u1=alo[ku1]; u1<alo[ku1+1]; u1++){
          int qu1 = m->q[ku1];
          int lu1 = bo->l[u1];
          double pu1 = qmd->p[qu1*(qmd->nLo)+lu1]/(2.0*lu1+1.0);
          for(int v1=u1; v1<alo[ku1+1]; v1++){
            for(int v=u; v<alo[ku+1]; v++){
              int uv   = mpos(u,v);
              int u1v1 = mpos(u1,v1);
              int vu1  = mpos(v,u1);
              int uv1  = mpos(u,v1);
              int uu1 = mpos(u,u1);
              int vv1 = mpos(v,v1);
              double D1 = Da[uv]*Db[u1v1];
              double D2 = Da[u1v1]*Db[uv];
              double D3 = Da[uv]*Da[u1v1] + Db[uv]*Db[u1v1];
              double D4 = Da[uv1]*Da[vu1] + Db[uv1]*Db[vu1];
              double D5 = Da[uu1]*Da[vv1] + Db[uu1]*Db[vv1];
              double D  = (D1+D2+D3) - 0.5*(D4+D5);
              D *= (u==v?1:2) * (u1==v1?1:2);
              // eq22 (R -> H -> E)
              if(u==v){
                D -= pu  * (u1==v1?1:2) * (Da[u1v1]+Db[u1v1]);
              }
              if(u1==v1){
                D -= pu1 * (u==v?1:2)   * (Da[uv]+Db[uv]);
              }
              double tg[3];
              R0_eq39_mmmm_grad(u,v,u1,v1, tg, m->r+3*ku, m->r+3*ku1, bo,qmd);
              r3scal(tg, D);
              r3add(g+3*ku , tg);
              r3min(g+3*ku1, tg);
            }
          }
        }
      }
    }
  }
  return;
}

void R_mmmm_grad_r(double * g, double * D, int * alo, basis * bo, mol * m, qmdata * qmd){
  // two-center corrections to one-center
  for(int ku=0; ku<m->n; ku++){
    int qu = m->q[ku];
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      int uu = mpos(u,u);
      int lu = bo->l[u];
      double pu = qmd->p[qu*(qmd->nLo)+lu ]/(2.0*lu +1.0);
      double Du = 2.0*D[uu];
      for(int u1=u; u1<alo[ku+1]; u1++){
        int lu1 = bo->l[u1];
        double pu1 = qmd->p[qu*(qmd->nLo)+lu1]/(2.0*lu1+1.0);
        int u1u1 = mpos(u1,u1);
        int uu1  = MPOSIF(u,u1);
        double D1  = D[uu] * D[u1u1] + D[uu] * D[u1u1];
        double D3  = D[uu1] * D[uu1] + D[uu1] * D[uu1];
        double Du1 = 2.0*D[u1u1];
        for(int k1=0; k1<m->n; k1++){
          if(k1==ku) continue;
          double dr[3];
          r3diff(dr, m->r+3*ku, m->r+3*k1);
          double r  = sqrt(r3dot(dr,dr));
          double gV = V_eq51_grad(lu, lu1, qu, m->q[k1], r, qmd);
          double fact1 = (u==u1?0.5:1.0) * (2.0*D1-D3);
          double fact2 = (u==u1 ? 0.5*Du*pu1 : Du*pu1+Du1*pu);  // eq22 (R -> H -> E)
          r3scal(dr, (fact1-fact2) * gV/r);
          r3add(g+3*ku, dr);
          r3min(g+3*k1, dr);
        }
      }
    }
  }
  // two-center
  for(int ku=0; ku<m->n; ku++){
    int qu = m->q[ku];
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      int lu = bo->l[u];
      double pu  = qmd->p[qu *(qmd->nLo)+lu ]/(2.0*lu +1.0);
      for(int ku1=ku+1; ku1<m->n; ku1++){
        for(int u1=alo[ku1]; u1<alo[ku1+1]; u1++){
          int qu1 = m->q[ku1];
          int lu1 = bo->l[u1];
          double pu1 = qmd->p[qu1*(qmd->nLo)+lu1]/(2.0*lu1+1.0);
          for(int v1=u1; v1<alo[ku1+1]; v1++){
            for(int v=u; v<alo[ku+1]; v++){
              int uv   = mpos(u,v);
              int u1v1 = mpos(u1,v1);
              int vu1  = mpos(v,u1);
              int uv1  = mpos(u,v1);
              int uu1 = mpos(u,u1);
              int vv1 = mpos(v,v1);
              double D1 = D[uv]*D[u1v1];
              double D4 = D[uv1]*D[vu1];
              double D5 = D[uu1]*D[vv1];
              double D0 = 4.0*D1 - D4 - D5;
              D0 *= (u==v?1:2) * (u1==v1?1:2);
              // eq22 (R -> H -> E)
              if(u==v){
                D0 -= pu  * (u1==v1?2.0:4.0) * D[u1v1];
              }
              if(u1==v1){
                D0 -= pu1 * (u==v?2.0:4.0)   * D[uv];
              }
              double tg[3];
              R0_eq39_mmmm_grad(u,v,u1,v1, tg, m->r+3*ku, m->r+3*ku1, bo,qmd);
              r3scal(tg, D0);
              r3add(g+3*ku , tg);
              r3min(g+3*ku1, tg);
            }
          }
        }
      }
    }
  }
  return;
}

