#include "eq.h"
#include "common.h"
#include "vec3.h"
#include "matrix.h"
#include "2el.h"

static double f_eq28_mm(int mu, int mv, int lu, int lv, int qu, int qv, euler * z, qmdata * qmd);
static double f_eq28_mp(int ma, int mv, int la, int lv, int qa, int qv, euler * z, qmdata * qmd);
static double f_eq30_mm(int mu, int mv, int lu, int lv, int qu, int qk, euler * z, qmdata * qmd);
static double f_eq30_mp(int ma, int mv, int la, int lv, int qa, int qk, euler * z, qmdata * qmd);
static double f_eq31(int mu, int mv, int lu, int lv, int qu, int qv, int qd, euler * zud, euler * zvd, qmdata * qmd);
static double S_eq32(int mu, int md, int lu, int ld, int qu, int qd, euler * zud, qmdata * qmd);

void H_eq22_mm(double * f, double * H, int * alo, double * mmmm, basis * bo, mol * m, qmdata * qmd){
  int Mo = bo->M;
  for(int ku=0; ku<m->n; ku++){
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      for(int v=u; v<alo[ku+1]; v++){
        int uv = mpos(u,v);
        double t = f[uv];
        for(int u1=0; u1<Mo; u1++){
          int ku1 = bo->k[u1];
          int qu1 = bo->Q[u1];
          int lu1 = bo->l[u1];
          double p1 = qmd->p[qu1*(qmd->nLo)+lu1]/(2.0*lu1+1.0);
          double r  = R(u,v,u1,u1, mmmm, bo, m, qmd);
          if(ku==ku1){
            r -= 0.5 * R(u,u1,u1,v, mmmm, bo, m, qmd);
          }
          t -= r * p1;
        }
        H[uv] = t;
      }
      for(int kv=ku+1; kv<m->n; kv++){
        for(int v=alo[kv]; v<alo[kv+1]; v++){
          int uv = mpos(u,v);
          H[uv] = f[uv];
        }
      }
    }
  }
  return;
}

void H_eq22_mp(double * f, double * H, int * alo, int * alv, double * pmmm, basis * bo, basis * bv, mol * m, qmdata * qmd){
  int Mo = bo->M;
  for(int ka=0; ka<m->n; ka++){
    for(int a=alv[ka]; a<alv[ka+1]; a++){
      for(int v=alo[ka]; v<alo[ka+1]; v++){
        int av = a*Mo+v;
        double t = f[av];
        for(int u1=0; u1<Mo; u1++){
          int ku1 = bo->k[u1];
          if(ku1==ka){
            continue;
          }
          int q1 = bo->Q[u1];
          int l1 = bo->l[u1];
          double p1 = qmd->p[q1*(qmd->nLo)+l1]/(2.0*l1+1.0);
          t -= p1 * R2(a,v,u1,u1, pmmm, bo,bv,m,qmd);
        }
        H[av] = t;
      }
      for(int kv=0; kv<m->n; kv++){
        if(kv==ka){
          continue;
        }
        for(int v=alo[kv]; v<alo[kv+1]; v++){
          int av = a*Mo+v;
          H[av] = f[av];
        }
      }
    }
  }
  return;
}

void f_eq25_mm(double * f, euler * z, int * alo, basis * bo, mol * m, qmdata * qmd){
  for(int ku=0; ku<m->n; ku++){
    int qu = m->q[ku];
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      int lu = bo->l[u];
      int mu = bo->m[u];
      // one-center
      for(int v=u; v<alo[ku+1]; v++){
        int lv = bo->l[v];
        int mv = bo->m[v];
        double f0 = 0.0;
        if( (lu == lv) && (mu == mv) ){
          f0 = qmd->fa[qu*(qmd->nLo)+lu]; // eq26
        }
        double fk = 0.0;
        for(int k1=0; k1<m->n; k1++){
          if(k1==ku){
            continue;
          }
          euler * zuk = z+(ku*m->n+k1);
          fk += f_eq30_mm(mu,mv, lu,lv, qu, m->q[k1], zuk, qmd);
        }
        f[mpos(u,v)] = f0 + fk;
      }
      // two-center
      for(int kv=ku+1; kv<m->n; kv++){
        int qv = m->q[kv];
        euler * zuv = z+(ku*m->n+kv);
        euler * zvu = z+(kv*m->n+ku);
        for(int v=alo[kv]; v<alo[kv+1]; v++){
          int lv = bo->l[v];
          int mv = bo->m[v];
          double f0;
          if(qu<=qv){
            f0 = f_eq28_mm(mu,mv, lu,lv, qu,qv, zuv, qmd);
          }
          else{
            f0 = f_eq28_mm(mv,mu, lv,lu, qv,qu, zvu, qmd);
          }
          double fk = 0.0;
          for(int k1=0; k1<m->n; k1++){
            if( (k1==ku) || (k1==kv) ){
              continue;
            }
            euler * zud = z+(ku*m->n+k1);
            euler * zvd = z+(kv*m->n+k1);
            fk += f_eq31(mu,mv, lu,lv, qu,qv, m->q[k1], zud,zvd, qmd);
          }
          f[mpos(u,v)] = f0+fk;
        }
      }
    }
  }
  return;
}

void f_eq25_mp(double * f, euler * z, int * alo, int * alv, basis * bo, basis * bv, mol * m, qmdata * qmd){
  int Mo = bo->M;
  for(int ka=0; ka<m->n; ka++){
    int qa = m->q[ka];
    for(int a=alv[ka]; a<alv[ka+1]; a++){
      int la = bv->l[a];
      int ma = bv->m[a];
      // one-center
      for(int v=alo[ka]; v<alo[ka+1]; v++){
        int lv = bo->l[v];
        int mv = bo->m[v];
        double fk = 0.0;
        for(int k1=0; k1<m->n; k1++){
          if(k1==ka){
            continue;
          }
          euler * zak = z+(ka*m->n+k1);
          fk += f_eq30_mp(ma,mv, la,lv, qa, m->q[k1], zak, qmd);
        }
        f[a*Mo+v] = fk;
      }
      // two-center
      for(int kv=0; kv<m->n; kv++){
        if(ka==kv){
          continue;
        }
        int qv = m->q[kv];
        euler * zav = z+(ka*m->n+kv);
        for(int v=alo[kv]; v<alo[kv+1]; v++){
          int lv = bo->l[v];
          int mv = bo->m[v];
          f[a*Mo+v] = f_eq28_mp(ma,mv, la,lv, qa,qv, zav, qmd);
        }
      }
    }
  }
  return;
}

static double f_eq28_mm(int mu, int mv, int lu, int lv, int qu, int qv, euler * z, qmdata * qmd){
  double r = z->r;
  int qq  = mpos(qu,qv); // qu<=qv
  int bra = qmd->qq_list[qq-1].fb;
  int ket = qmd->qq_list[qq  ].fb;
  double f0 = 0;
  for(int i=bra; i<ket; i++){
    if( (qmd->fb[i].lu == lu) && (qmd->fb[i].lv == lv) ){
      int m     = qmd->fb[i].m;
      int aph   = (((lu+m)%2)?(-1):(1));
      double F  = F_eq47(i, lu, lv, qu, qv, r, qmd);
      double AA = A_new(lu, m,mu,z) * A_new(lv, m,mv,z);
      if(m){
        AA     += A_new(lu,-m,mu,z) * A_new(lv,-m,mv,z);
      }
      f0 += aph * F * AA;
    }
  }
  return f0;
}

static double f_eq28_mp(int ma, int mv, int la, int lv, int qa, int qv, euler * z, qmdata * qmd){
  double r = z->r;
  int qq  = MPOSIF(qa,qv);
  int bra = qmd->qq_list[qq-1].f1b;
  int ket = qmd->qq_list[qq  ].f1b;
  double f0 = 0;
  for(int i=bra; i<ket; i++){
    if( (qmd->f1b[i].qu == qa) &&
        (qmd->f1b[i].qv == qv) &&
        (qmd->f1b[i].lu == la) &&
        (qmd->f1b[i].lv == lv) ){
      int m   = qmd->f1b[i].m;
      int aph = (((la+m)%2)?(-1):(1));
      double F = F_eq48(i, la,lv, qa,qv, r, qmd);
      double AA = A_new(la, m,ma,z) * A_new(lv, m,mv,z);
      if(m){
        AA     += A_new(la,-m,ma,z) * A_new(lv,-m,mv,z);
      }
      f0 += AA * F * aph;
    }
  }
  return f0;
}

static double S_eq32(int mu, int md, int lu, int ld, int qu, int qd, euler * zud, qmdata * qmd){
  double r = zud->r;
  int qq  = MPOSIF(qu,qd);
  int bra = qmd->qq_list[qq-1].fo;
  int ket = qmd->qq_list[qq  ].fo;
  double s = 0.0;
  for(int i=bra; i<ket; i++){
    if( (qmd->fo[i].qu != qu) ||
        (qmd->fo[i].qv != qd) ||
        (qmd->fo[i].lu != lu) ||
        (qmd->fo[i].lv != ld) ){
      continue;
    }
    int m   = qmd->fo[i].m;
    int aph = (lu+m)%2?-1:1;
    double S  = S_eq50(i, lu, ld, qu, qd, r, qmd);
    double AA = A_new(lu, m,mu,zud) * A_new(ld, m,md,zud);
    if(m){
      AA     += A_new(lu,-m,mu,zud) * A_new(ld,-m,md,zud);
    }
    s += aph*AA*S;
  }
  return s;
}

static double f_eq30_mm(int mu, int mv, int lu, int lv, int qu, int qk, euler * z, qmdata * qmd){

  if(lu>lv){
    int t;
    SWAP(lu,lv,t);
    SWAP(mu,mv,t);
  }
  double r = z->r;

  int qq   = MPOSIF(qk,qu);
  int bra  = qmd->qq_list[qq-1].ub;
  int ket  = qmd->qq_list[qq  ].ub;
  int lb   = MIN(lu,lv);
  double s = 0.0;
  for(int i=bra; i<ket; i++){
    if( (qmd->ub[i].qu == qu) &&
        (qmd->ub[i].qv == qk) &&
        (qmd->ub[i].lu == lu) &&
        (qmd->ub[i].lv == lv) ){
      int l = qmd->ub[i].l;
      double V = V_eq49_mm(i,lu,lv,qu,qk,r,qmd);
      double AAB = 0.0;
      for(int m=-lb; m<=lb; m++){
        double A1 = A_new(lu, m, mu, z);
        double A2 = A_new(lv, m, mv, z);
        double B0 = B(l,lu,lv,0,m,m);
        AAB += A1*A2*B0;
      }
      s += AAB*V;
    }
  }
  return s;
}

static double f_eq30_mp(int ma, int mv, int la, int lv, int qa, int qk, euler * z, qmdata * qmd){
  double r = z->r;
  int qq   = MPOSIF(qk,qa);
  int bra  = qmd->qq_list[qq-1].u1b;
  int ket  = qmd->qq_list[qq  ].u1b;
  int lb   = MIN(la,lv);
  double s = 0.0;
  for(int i=bra; i<ket; i++){
    if( (qmd->u1b[i].qu == qa) &&
        (qmd->u1b[i].qv == qk) &&
        (qmd->u1b[i].lu == la) &&
        (qmd->u1b[i].lv == lv) ){
      int l  = qmd->u1b[i].l;
      double V = V_eq49_mp(i,la,lv,qa,qk,r,qmd);
      double AAB = 0.0;
      for(int m=-lb; m<=lb; m++){
        double A1 = A_new(la, m, ma, z);
        double A2 = A_new(lv, m, mv, z);
        double B0 = B(l,la,lv,0,m,m);
        AAB += A1*A2*B0;
      }
      s += AAB*V;
    }
  }
  return s;
}

static double f_eq31(int mu, int mv, int lu, int lv, int qu, int qv, int qd, euler * zud, euler * zvd, qmdata * qmd){
  int Lp = qmd->Lp[qd];
  double s = 0.0;
  for(int ld=0; ld<=Lp; ld++){
    double fd = qmd->cf[qd*qmd->nLp+ld];
    for(int md=-ld; md<=ld; md++){
      double Sud = S_eq32(mu, md, lu, ld, qu, qd, zud, qmd);
      double Svd = S_eq32(mv, md, lv, ld, qv, qd, zvd, qmd);
      s += Sud*Svd*fd;
    }
  }
  return s;
}

