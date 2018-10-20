#include "qm.h"
#include "matrix.h"
#include "eq.h"
#include "gradient.h"
#include "vec3.h"

static void f_eq28_mm_grad(int mu, int mv, int lu, int lv, int qu, int qv, double g[3], double ru[3], double rv[3], qmdata * qmd){

  double z[3];
  r3diff(z, ru, rv);
  double r  = sqrt(r3dot(z,z));;
  double r1 = 1.0/r;
  r3scal(z, r1);

  int qq  = mpos(qu,qv); // qu<=qv
  int bra = qmd->qq_list[qq-1].fb;
  int ket = qmd->qq_list[qq  ].fb;
  double g1[3] = {0};
  double g2[3] = {0};
  for(int i=bra; i<ket; i++){
    if( (qmd->fb[i].lu == lu) && (qmd->fb[i].lv == lv) ){
      int m     = qmd->fb[i].m;
      int aph   = (((lu+m)%2)?(-1):(1));
      double gA1[3], A1 = A_grad_z(lu, m, mu, gA1, z);
      double gA2[3], A2 = A_grad_z(lv, m, mv, gA2, z);
      double AA = A1*A2;
      fgpair_t F = F_eq47_grad(i, lu, lv, qu, qv, r, qmd);
      double tg[3];
      r3sums(tg, gA1, A2, gA2, A1);
      if(m){
        double gA1_[3], A1_ = A_grad_z(lu, -m, mu, gA1_, z);
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

static void f_eq30_mm_grad(int mu, int mv, int lu, int lv, int qu, int qk, double * g, double ru[3], double rk[3], qmdata * qmd){

  double z[3];
  r3diff(z, ru, rk);
  double r  = sqrt(r3dot(z,z));
  double r1 = 1.0/r;
  r3scal(z, r1);

  if(lu>lv){
    int t;
    SWAP(lu,lv,t);
    SWAP(mu,mv,t);
  }

  int qq  = MPOSIF(qk,qu);
  int bra = qmd->qq_list[qq-1].ub;
  int ket = qmd->qq_list[qq  ].ub;
  int lb  = MIN(lu,lv);
  double g1[3] = {0};
  double g2[3] = {0};
  for(int i=bra; i<ket; i++){
    if( (qmd->ub[i].qu == qu) &&
        (qmd->ub[i].qv == qk) &&
        (qmd->ub[i].lu == lu) &&
        (qmd->ub[i].lv == lv) ){
      int l = qmd->ub[i].l;
      fgpair_t V = V_eq49_mm_grad(i,lu,lv,qu,qk,r,qmd);
      double AAB = 0.0;
      double tg[3] = {0};
      for(int m=-lb; m<=lb; m++){
        double gA1[3], A1 = A_grad_z(lu, m, mu, gA1, z);
        double gA2[3], A2 = A_grad_z(lv, m, mv, gA2, z);
        double B0 = B(l,lu,lv,0,m,m);
        AAB += A1*A2*B0;
        r3adds(tg, gA1, A2*B0);
        r3adds(tg, gA2, A1*B0);
      }
      r3adds(g1, z, AAB*V.g);
      r3adds(g2, tg,    V.f);
    }
  }
  A_grad_z2r(g2, z, r1);
  r3sum(g, g1, g2);
  return;
}

static double S_eq32_grad(int mu, int md, int lu, int ld, int qu, int qd, double g[3], double ru[3], double rd[3], qmdata * qmd){

  double z[3];
  r3diff(z, ru, rd);
  double r  = sqrt(r3dot(z,z));;
  double r1 = 1.0/r;
  r3scal(z, r1);

  int qq  = MPOSIF(qu,qd);
  int bra = qmd->qq_list[qq-1].fo;
  int ket = qmd->qq_list[qq  ].fo;
  double g1[3] = {0};
  double g2[3] = {0};
  double s = 0;
  for(int i=bra; i<ket; i++){
    if( (qmd->fo[i].qu != qu) ||
        (qmd->fo[i].qv != qd) ||
        (qmd->fo[i].lu != lu) ||
        (qmd->fo[i].lv != ld) ){
      continue;
    }
    int m   = qmd->fo[i].m;
    int aph = (lu+m)%2?-1:1;
    fgpair_t S = S_eq50_grad(i, lu, ld, qu, qd, r, qmd);
    double gA1[3], A1 = A_grad_z(lu, m, mu, gA1, z);
    double gA2[3], A2 = A_grad_z(ld, m, md, gA2, z);
    double AA = A1*A2;
    double tg[3];
    r3sums(tg, gA1, A2, gA2, A1);
    if(m){
      double gA1_[3], A1_ = A_grad_z(lu, -m, mu, gA1_, z);
      double gA2_[3], A2_ = A_grad_z(ld, -m, md, gA2_, z);
      AA += A1_*A2_;
      r3adds(tg, gA1_, A2_);
      r3adds(tg, gA2_, A1_);
    }
    r3adds(g1, z, AA*aph*S.g);
    r3adds(g2, tg,   aph*S.f);
    s += aph*AA*S.f;
  }
  A_grad_z2r(g2, z, r1);
  r3sum(g, g1, g2);
  return s;
}

static void f_eq31_grad(int mu, int mv, int lu, int lv, int qu, int qv, int qd,
    double gu[3], double gv[3], double ru[3], double rv[3], double rd[3], qmdata * qmd){
  r3set(gu, 0.0);
  r3set(gv, 0.0);
  int Lp = qmd->Lp[qd];
  for(int ld=0; ld<=Lp; ld++){
    double fd = qmd->cf[qd*qmd->nLp+ld];
    for(int md=-ld; md<=ld; md++){
      double gSud[3], Sud = S_eq32_grad(mu, md, lu, ld, qu, qd, gSud, ru, rd, qmd);
      double gSvd[3], Svd = S_eq32_grad(mv, md, lv, ld, qv, qd, gSvd, rv, rd, qmd);
      r3adds(gu, gSud, Svd*fd);
      r3adds(gv, gSvd, Sud*fd);
    }
  }
  return;
}

void f_eq25_mm_grad(double * g, double * Da, double * Db, int * alo, basis * bo, mol * m, qmdata * qmd){
  for(int ku=0; ku<m->n; ku++){
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      for(int v=u; v<alo[ku+1]; v++){
        double duv = (Da[mpos(u,v)]+Db[mpos(u,v)])*(u==v?1.0:2.0);
        for(int k=0; k<m->n; k++){
          if(k==ku) continue;
          int qu = m->q[ku];
          int lu = bo->l[u];
          int mu = bo->m[u];
          int lv = bo->l[v];
          int mv = bo->m[v];
          int qk = m->q[k];
          double tg[3];
          f_eq30_mm_grad(mu, mv, lu, lv, qu, qk, tg, m->r+3*ku, m->r+3*k, qmd);
          r3scal(tg, duv);
          r3add(g+3*ku, tg);
          r3min(g+3*k,  tg);
        }
      }
    }
  }

  for(int ku=0; ku<m->n; ku++){
    for(int kv=ku+1; kv<m->n; kv++){
      int qu = m->q[ku];
      int qv = m->q[kv];
      double * ru = m->r+3*ku;
      double * rv = m->r+3*kv;
      for(int u=alo[ku]; u<alo[ku+1]; u++){
        for(int v=alo[kv]; v<alo[kv+1]; v++){
          int lu = bo->l[u];
          int mu = bo->m[u];
          int lv = bo->l[v];
          int mv = bo->m[v];
          double duv = 2.0*(Da[mpos(u,v)]+Db[mpos(u,v)]);
          double tg[3];
          if(qu<=qv){
            f_eq28_mm_grad(mu, mv, lu, lv, qu, qv, tg, ru, rv, qmd);
            r3scal(tg, duv);
          }
          else{
            f_eq28_mm_grad(mv, mu, lv, lu, qv, qu, tg, rv, ru, qmd);
            r3scal(tg, -duv);
          }
          r3add(g+3*ku, tg);
          r3min(g+3*kv, tg);
          for(int k1=0; k1<m->n; k1++){
            if( (k1==ku) || (k1==kv) ){
              continue;
            }
            double gu[3], gv[3];
            f_eq31_grad(mu, mv, lu, lv, qu, qv, m->q[k1], gu, gv, ru, rv, m->r+3*k1, qmd);
            r3scal(gu, duv);
            r3scal(gv, duv);
            r3add(g+3*ku, gu);
            r3add(g+3*kv, gv);
            r3min(g+3*k1, gu);
            r3min(g+3*k1, gv);
          }
        }
      }
    }
  }
  return;
}

void f_eq25_mm_grad_r(double * g, double * D, int * alo, basis * bo, mol * m, qmdata * qmd){
  for(int ku=0; ku<m->n; ku++){
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      for(int v=u; v<alo[ku+1]; v++){
        double duv = D[mpos(u,v)] * (u==v?2.0:4.0);
        for(int k=0; k<m->n; k++){
          if(k==ku) continue;
          int qu = m->q[ku];
          int lu = bo->l[u];
          int mu = bo->m[u];
          int lv = bo->l[v];
          int mv = bo->m[v];
          int qk = m->q[k];
          double tg[3];
          f_eq30_mm_grad(mu, mv, lu, lv, qu, qk, tg, m->r+3*ku, m->r+3*k, qmd);
          r3scal(tg, duv);
          r3add(g+3*ku, tg);
          r3min(g+3*k,  tg);
        }
      }
    }
  }

  for(int ku=0; ku<m->n; ku++){
    for(int kv=ku+1; kv<m->n; kv++){
      int qu = m->q[ku];
      int qv = m->q[kv];
      double * ru = m->r+3*ku;
      double * rv = m->r+3*kv;
      for(int u=alo[ku]; u<alo[ku+1]; u++){
        for(int v=alo[kv]; v<alo[kv+1]; v++){
          int lu = bo->l[u];
          int mu = bo->m[u];
          int lv = bo->l[v];
          int mv = bo->m[v];
          double duv = 4.0 * D[mpos(u,v)];
          double tg[3];
          if(qu<=qv){
            f_eq28_mm_grad(mu, mv, lu, lv, qu, qv, tg, ru, rv, qmd);
            r3scal(tg, duv);
          }
          else{
            f_eq28_mm_grad(mv, mu, lv, lu, qv, qu, tg, rv, ru, qmd);
            r3scal(tg, -duv);
          }
          r3add(g+3*ku, tg);
          r3min(g+3*kv, tg);
          for(int k1=0; k1<m->n; k1++){
            if( (k1==ku) || (k1==kv) ){
              continue;
            }
            double gu[3], gv[3];
            f_eq31_grad(mu, mv, lu, lv, qu, qv, m->q[k1], gu, gv, ru, rv, m->r+3*k1, qmd);
            r3scal(gu, duv);
            r3scal(gv, duv);
            r3add(g+3*ku, gu);
            r3add(g+3*kv, gv);
            r3min(g+3*k1, gu);
            r3min(g+3*k1, gv);
          }
        }
      }
    }
  }
  return;
}

