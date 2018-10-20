#include "qm.h"
#include "matrix.h"
#include "eq.h"
#include "2el.h"
#include "tools.h"
#include "gradient.h"
#include "vec3.h"
#define EPS 1e-15

static void f_eq30_mp_grad(int ma, int mv, int la, int lv, int qa, int qk, double * g, double ra[3], double rk[3], qmdata * qmd){

  double z[3];
  r3diff(z, ra, rk);
  double r  = sqrt(r3dot(z,z));
  double r1 = 1.0/r;
  r3scal(z, r1);

  int qq   = MPOSIF(qk,qa);
  int bra  = qmd->qq_list[qq-1].u1b;
  int ket  = qmd->qq_list[qq  ].u1b;
  int lb   = MIN(la,lv);
  double g1[3] = {0};
  double g2[3] = {0};
  for(int i=bra; i<ket; i++){
    if( (qmd->u1b[i].qu == qa) &&
        (qmd->u1b[i].qv == qk) &&
        (qmd->u1b[i].lu == la) &&
        (qmd->u1b[i].lv == lv) ){
      int l  = qmd->u1b[i].l;
      fgpair_t V = V_eq49_mp_grad(i,la,lv,qa,qk,r,qmd);
      double AAB = 0.0;
      double tg[3] = {0};
      for(int m=-lb; m<=lb; m++){
        double gA1[3], A1 = A_grad_z(la, m, ma, gA1, z);
        double gA2[3], A2 = A_grad_z(lv, m, mv, gA2, z);
        double B0 = B(l,la,lv,0,m,m);
        AAB += A1*A2*B0;
        r3adds(tg, gA1, A2*B0);
        r3adds(tg, gA2, A1*B0);
      }
      r3adds(g1, z, AAB*V.g);
      r3adds(g2, tg,    V.f);
    }
  }
  A_grad_z2r(g2, z, r1);
  r3add(g, g1);
  r3add(g, g2);
  return;
}

static void f_eq28_mp_grad(int ma, int mv, int la, int lv, int qa, int qv, double g[3], double ra[3], double rv[3], qmdata * qmd){

  double z[3];
  r3diff(z, ra, rv);
  double r  = sqrt(r3dot(z,z));;
  double r1 = 1.0/r;
  r3scal(z, r1);

  int qq  = MPOSIF(qa,qv);
  int bra = qmd->qq_list[qq-1].f1b;
  int ket = qmd->qq_list[qq  ].f1b;
  double g1[3] = {0};
  double g2[3] = {0};
  for(int i=bra; i<ket; i++){
    if( (qmd->f1b[i].qu == qa) &&
        (qmd->f1b[i].qv == qv) &&
        (qmd->f1b[i].lu == la) &&
        (qmd->f1b[i].lv == lv) ){
      int m   = qmd->f1b[i].m;
      int aph = (((la+m)%2)?(-1):(1));
      double gA1[3], A1 = A_grad_z(la, m, ma, gA1, z);
      double gA2[3], A2 = A_grad_z(lv, m, mv, gA2, z);
      double AA = A1*A2;
      fgpair_t F = F_eq48_grad(i, la,lv, qa,qv, r, qmd);
      double tg[3];
      r3sums(tg, gA1, A2, gA2, A1);
      if(m){
        double gA1_[3], A1_ = A_grad_z(la, -m, ma, gA1_, z);
        double gA2_[3], A2_ = A_grad_z(lv, -m, mv, gA2_, z);
        AA += A1_*A2_;
        r3adds(tg, gA1_, A2_);
        r3adds(tg, gA2_, A1_);
      }
      r3adds(g1, z, AA*aph*F.g);
      r3adds(g2, tg,   aph*F.f);
    }
  }
  A_grad_z2r(g2, z, r1);
  r3sum(g, g1, g2);
  return;
}

static void R0_eq39_mmmp_grad(int a, int v, int u1, int v1,
    double g[3], double ra[3], double ru1[3],
    basis * bo,  basis * bv,   qmdata * qmd){

  double z[3];
  r3diff(z, ra, ru1);
  double r  = sqrt(r3dot(z,z));
  double r1 = 1.0/r;
  r3scal(z, r1);

  int qa  = bv->Q[a ];
  int qu1 = bo->Q[u1];
  int ma  = bv->m[a ];
  int la  = bv->l[a ];
  int mv  = bo->m[v ];
  int lv  = bo->l[v ];
  int mu1 = bo->m[u1];
  int lu1 = bo->l[u1];
  int mv1 = bo->m[v1];
  int lv1 = bo->l[v1];

  double A1, A2, A3, A4;
  double gA1[3], gA2[3], gA3[3], gA4[3];
  double g1[3]={0}, g2[3]={0};

  for(int m_=-la;m_<=la;m_++){
    A1 = A_grad_z(la,m_,ma,gA1,z);
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

          for(int l=abs(la-lv); l<=la+lv; l++){
            double q1;
            if( ! qlll_pm(qa,  la,  lv,  l,  &q1, qmd)){
              continue;
            }

            for(int l1=abs(lu1-lv1); l1<=lu1+lv1; l1++){
              double q2;
              if( ! qlll_mm(qu1, lu1, lv1, l1, &q2, qmd)){
                continue;
              }
              int lm = MIN(l,l1);
              for(int m=-lm; m<=lm; m++){
                double B1 = B(l,  la,  lv,  m, m_,  m__);
                if(fabs(B1)<EPS){
                  continue;
                }
                double B2 = B(l1, lu1, lv1, m, m1_, m1__);
                if(fabs(B2)<EPS){
                  continue;
                }
                fgpair_t G = G_eq52_mmmp_grad(m,l,l1,la,lv,lu1,lv1,qa,qu1,q1,q2,r,qmd);
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

static void dFav_dRk_1c(int a, int v, int k,
    double g[3], double * Da, double * Db,
    int * alo, basis * bo, basis * bv, mol * m, qmdata * qmd){
  int ka = bv->k[a];
  int qa = bv->Q[a];
  int la = bv->l[a];
  int ma = bv->m[a];
  int lv = bo->l[v];
  int mv = bo->m[v];
  f_eq30_mp_grad(ma, mv, la, lv, qa, m->q[k], g, m->r+3*ka, m->r+3*k, qmd);
  for(int u1=alo[k]; u1<alo[k+1]; u1++){
    for(int v1=u1; v1<alo[k+1]; v1++){
      int u1v1 = mpos(u1,v1);
      double tg[3];
      R0_eq39_mmmp_grad(a, v, u1, v1, tg, m->r+3*ka, m->r+3*k, bo,  bv, qmd);
      double D = Da[u1v1]+Db[u1v1];
      if(u1!=v1){
        D *= 2.0;
      }
      else{
        int q1 = bo->Q[u1];
        int l1 = bo->l[u1];
        D -= qmd->p[q1*(qmd->nLo)+l1]/(2.0*l1+1.0);
      }
      r3adds(g, tg, D);
    }
  }
  return;
}

static void dFav_dRk_2c(int a, int v,
    double ga[3], double gb[3], double * Da, double * Db,
    int * alo, basis * bo, basis * bv, mol * m, qmdata * qmd){

  int ka = bv->k[a];
  int qa = bv->Q[a];
  int la = bv->l[a];
  int ma = bv->m[a];
  int kv = bo->k[v];
  int qv = bo->Q[v];
  int lv = bo->l[v];
  int mv = bo->m[v];

  f_eq28_mp_grad(ma,mv, la,lv, qa,qv, ga, m->r+3*ka, m->r+3*kv, qmd);
  r3cp(gb, ga);

  for(int v1=alo[ka]; v1<alo[ka+1]; v1++){
    for(int u1=alo[kv]; u1<alo[kv+1]; u1++){
      int u1v1 = MPOSIF(u1,v1);
      double tg[3];
      R0_eq39_mmmp_grad(a, v1, u1, v, tg, m->r+3*ka, m->r+3*kv, bo,  bv, qmd);
      r3adds(ga, tg, -Da[u1v1]);
      r3adds(gb, tg, -Db[u1v1]);
    }
  }
  return;
}

void E2_grad(double * g,
    double * Da, double * Db,
    double * Fa, double * Fb,
    double * Xa, double * Xb,
    double * sa, double * sb,
    double * ja, double * jb,
    int * alo, int * alv, basis * bo , basis * bv, mol * m, qmdata * qmd){

  int Mo = bo->M;
  int Mv = bv->M;

  double * ta = malloc(sizeof(double)*Mo*Mv);
  double * tb = malloc(sizeof(double)*Mo*Mv);

  for(int v=0; v<Mo; v++){
    int lv = bo->l[v];
    double clv = 2.0 / (1.0 + 2.0 * lv);
    double ev  = qmd->fa [bo->Q[v ]*qmd->nLo+lv ];
    double multa = sa[v] * sa[v] * ja[v] * clv;
    double multb = sb[v] * sb[v] * jb[v] * clv;
    for(int a=0; a<Mv; a++){
      int av = a*Mo+v;
      double ea  = qmd->f2a[bv->Q[a ]*qmd->nLv+bv->l[a ]];
      double eva = 1.0/(ev-ea);
      double fa = 0.0;
      double fb = 0.0;
      for(int u=0; u<Mo; u++){
        int uv = MPOSIF(u,v);
        int au = a*Mo + u;
        fa += Da[uv] * (sa[v]*eva*Fa[au] + sa[u]*Xa[au]);
        fb += Db[uv] * (sb[v]*eva*Fb[au] + sb[u]*Xb[au]);
      }
      ta[av] = fa - Xa[av] * eva * multa;
      tb[av] = fb - Xb[av] * eva * multb;
    }
  }

  for(int k=0; k<m->n; k++){
    for(int ka=0; ka<m->n; ka++){
      if(k==ka) continue;
      for(int a=alv[ka]; a<alv[ka+1]; a++){
        for(int v=alo[ka]; v<alo[ka+1]; v++){
          int av = a*Mo+v;
          double tg[3]={0};
          dFav_dRk_1c(a, v, k, tg, Da, Db, alo, bo, bv, m, qmd);
          r3scal(tg, ta[av]+tb[av]);
          r3add(g+3*ka, tg);
          r3min(g+3*k,  tg);
        }
      }
    }
  }

  for(int ka=0; ka<m->n; ka++){
    for(int kv=0; kv<m->n; kv++){
      if(kv==ka) continue;
      for(int a=alv[ka]; a<alv[ka+1]; a++){
        for(int v=alo[kv]; v<alo[kv+1]; v++){
          int av = a*Mo+v;
          double ga[3], gb[3];
          dFav_dRk_2c(a, v, ga, gb, Da, Db, alo, bo, bv, m, qmd);
          double tg[3];
          r3sums(tg, ga, ta[av], gb, tb[av]);
          r3add(g+3*ka, tg);
          r3min(g+3*kv, tg);
        }
      }
    }
  }

  free(ta);
  free(tb);
  return;
}

void E2_grad_r(double * g,
    double * D, double * F, double * X, double * s, double * j,
    int * alo, int * alv, basis * bo , basis * bv, mol * m, qmdata * qmd){

  int Mo = bo->M;
  int Mv = bv->M;

  double * t = malloc(sizeof(double)*Mo*Mv);

  for(int v=0; v<Mo; v++){
    int lv = bo->l[v];
    double clv  = 2.0 / (1.0 + 2.0 * lv);
    double ev   = qmd->fa [bo->Q[v ]*qmd->nLo+lv ];
    double mult = s[v] * s[v] * j[v] * clv;
    for(int a=0; a<Mv; a++){
      int av = a*Mo+v;
      double ea  = qmd->f2a[bv->Q[a ]*qmd->nLv+bv->l[a ]];
      double eva = 1.0/(ev-ea);
      double f = 0.0;
      for(int u=0; u<Mo; u++){
        int uv = MPOSIF(u,v);
        int au = a*Mo + u;
        f += D[uv] * (s[v]*eva*F[au] + s[u]*X[au]);
      }
      t[av] = f - X[av] * eva * mult;
    }
  }

  for(int k=0; k<m->n; k++){
    for(int ka=0; ka<m->n; ka++){
      if(k==ka) continue;
      for(int a=alv[ka]; a<alv[ka+1]; a++){
        for(int v=alo[ka]; v<alo[ka+1]; v++){
          int av = a*Mo+v;
          double tg[3]={0};
          dFav_dRk_1c(a, v, k, tg, D, D, alo, bo, bv, m, qmd);
          r3scal(tg, 2.0*t[av]);
          r3add(g+3*ka, tg);
          r3min(g+3*k,  tg);
        }
      }
    }
  }

  for(int ka=0; ka<m->n; ka++){
    for(int kv=0; kv<m->n; kv++){
      if(kv==ka) continue;
      for(int a=alv[ka]; a<alv[ka+1]; a++){
        for(int v=alo[kv]; v<alo[kv+1]; v++){
          int av = a*Mo+v;
          double ga[3], gb[3];
          dFav_dRk_2c(a, v, ga, gb, D, D, alo, bo, bv, m, qmd);
          double tg[3];
          r3sums(tg, ga, t[av], gb, t[av]);
          r3add(g+3*ka, tg);
          r3min(g+3*kv, tg);
        }
      }
    }
  }

  free(t);
  return;
}

