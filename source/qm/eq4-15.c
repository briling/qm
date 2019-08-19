#include "eq.h"
#include "matrix.h"
#include "2el.h"

void D_eq9(int N, int M, double * C, double * D){
  for(int i=0; i<M; i++){
    for(int j=i; j<M; j++){
      double s = 0.0;
      for(int k=0; k<N; k++){
        s += C[k*M+i]*C[k*M+j];
      }
      D[mpos(i,j)] = s;
    }
  }
  return;
}

void F_eq4(double * Da, double * Db,
           double * H, double * Fa, double * Fb,
           int * alo, double * mmmm,
           basis * bo, mol * m, qmdata * qmd){
  int Mo = bo->M;
  for(int i=0; i<symsize(Mo); i++){
    Fa[i] = H[i];
    Fb[i] = H[i];
  }
  for(int ku=0; ku<m->n; ku++){
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      for(int v=u; v<alo[ku+1]; v++){
        int uv = mpos(u,v);
        double df = 0.0;
        for(int ku1=0; ku1<m->n; ku1++){
          for(int u1=alo[ku1]; u1<alo[ku1+1]; u1++){
            int v1 = u1;
            int i = mpos(u1,v1);
            double R1234 = R(u,v,u1,v1, mmmm, bo,m,qmd);
            df += R1234 * (Da[i] + Db[i]);
            for(int v1=u1+1; v1<alo[ku1+1]; v1++){
              int i = mpos(u1,v1);
              double R1234 = R(u,v,u1,v1, mmmm, bo,m,qmd);
              df += 2.0 * R1234 * (Da[i] + Db[i]);
            }
          }
        }
        Fa[uv] += df;
        Fb[uv] += df;
      }
      for(int v=u; v<Mo; v++){
        int kv = bo->k[v];
        int uv = mpos(u,v);
        double ta = 0.0;
        double tb = 0.0;
        for(int u1=alo[kv]; u1<alo[kv+1]; u1++){
          for(int v1=alo[ku]; v1<alo[ku+1]; v1++){
            int i = MPOSIF(u1,v1);
            double R1432 = R(u,v1,u1,v, mmmm, bo,m,qmd);
            ta -= R1432 * Da[i];
            tb -= R1432 * Db[i];
          }
        }
        Fa[uv] += ta;
        Fb[uv] += tb;
      }
    }
  }
  return;
}

void FX(int Mo, int Mv, double * Fa, double * Fb, double * Xa, double * Xb, double * FaXa, double * FbXb){
  for(int u=0; u<Mo; u++){
    for(int v=0; v<Mo; v++){
      double da = 0.0;
      double db = 0.0;
      for(int a=0; a<Mv; a++){
        int au = a*Mo+u;
        int av = a*Mo+v;
        da += Fa[au] * Xa[av];
        db += Fb[au] * Xb[av];
      }
      FaXa[u*Mo+v] = da;
      FbXb[u*Mo+v] = db;
    }
  }
  return;
}

void F2_eq6(int Mo, double * FaXa, double * FbXb, double * sa, double * sb, double * F2a, double * F2b){
  for(int v=0; v<Mo; v++){
    for(int u=0; u<=v; u++){
      int uv_sym = mpos(u,v);
      int uv = Mo*u+v;
      int vu = Mo*v+u;
      F2a[uv_sym] = ( sa[u]*FaXa[vu] + sa[v]*FaXa[uv] ) * 0.5;
      F2b[uv_sym] = ( sb[u]*FbXb[vu] + sb[v]*FbXb[uv] ) * 0.5;
    }
  }
  return;
}

void X_eq7(double * Fa, double * Fb, double * Xa, double * Xb, basis * bo, basis * bv, qmdata * qmd){
  int Mo = bo->M;
  int Mv = bv->M;
  for(int a=0; a<Mv; a++){
    int qa = bv->Q[a];
    int la = bv->l[a];
    double ea = qmd->f2a[qa*qmd->nLv+la]; //eq41
    for(int v=0; v<Mo; v++){
      int    av = a*Mo+v;
      int    qv = bo->Q[v];
      int    lv = bo->l[v];
      double ev = qmd->fa[qv*qmd->nLo+lv]; //eq41
      double denom = 1.0/(ev-ea);
      Xa[av] = Fa[av]*denom;
      Xb[av] = Fb[av]*denom;
    }
  }
  return;
}

void F_eq8(double * Da, double * Db,
    double * H, double * Fa, double * Fb,
    int * alo, int * alv, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd){
  int Mo = bo->M;
  int Mv = bv->M;
  for(int i=0; i<Mo*Mv; i++){
    Fa[i] = H[i];
    Fb[i] = H[i];
  }
  for(int ka=0; ka<m->n; ka++){
    for(int a=alv[ka]; a<alv[ka+1]; a++){
      for(int v=alo[ka]; v<alo[ka+1]; v++){
        int av = a*Mo+v;
        double df = 0.0;
        for(int ku1=0; ku1<m->n; ku1++){
          if(ka==ku1){
            continue;
          }
          for(int u1=alo[ku1]; u1<alo[ku1+1]; u1++){
            int v1 = u1;
            int i = mpos(u1,v1);
            double R1234 = R2(a, v, u1, v1, pmmm, bo, bv, m, qmd);
            df += R1234 * (Da[i] + Db[i]);
            for(int v1=u1+1; v1<alo[ku1+1]; v1++){
              int i = mpos(u1,v1);
              double R1234 = R2(a, v, u1, v1, pmmm, bo, bv, m, qmd);
              df += 2.0 * R1234 * (Da[i] + Db[i]);
            }
          }
        }
        Fa[av] += df;
        Fb[av] += df;
      }
      for(int v=0; v<Mo; v++){
        int kv = bo->k[v];
        if(kv==ka){
          continue;
        }
        int av = a*Mo+v;
        double ta = 0.0;
        double tb = 0.0;
        for(int u1=alo[kv]; u1<alo[kv+1]; u1++){
          for(int v1=alo[ka]; v1<alo[ka+1]; v1++){
            int i = MPOSIF(u1,v1);
            double R1432 = R2(a, v1, u1, v,  pmmm, bo, bv, m, qmd);
            ta -= R1432 * Da[i];
            tb -= R1432 * Db[i];
          }
        }
        Fa[av] += ta;
        Fb[av] += tb;
      }
    }
  }
  return;
}

void s_eq15(int Mv, double * X, double * s, int * alo, basis * bo, mol * m, qmdata * qmd){
  int Mo = bo->M;
  for(int ku=0; ku<m->n; ku++){
    int q = m->q[ku];
    int l = qmd->Lo[q];
    for(int lu=0; lu<=l; lu++){
      double t = 0.0;
      int bra = alo[ku] + lu * lu;
      int ket = alo[ku] + (lu+1)*(lu+1);
      for(int v=bra; v<ket; v++){ /* functions on atom ku with L==lu */
        for(int a=0; a<Mv; a++){
          double x = X[a*Mo+v];
          t += x*x;
        }
      }
      t = (1.0+2.0*lu) / (1.0+2.0*lu + t); /* ( 1 + t/(1+2*lu) )^{-1} */
      for(int v=bra; v<ket; v++){
        s[v] = t;
      }
    }
  }
  return;
}

void F_eq4_r(double * D, double * H, double * F,
             int * alo, double * mmmm,
             basis * bo, mol * m, qmdata * qmd){
  int Mo = bo->M;
  for(int i=0; i<symsize(Mo); i++){
    F[i] = H[i];
  }
  for(int ku=0; ku<m->n; ku++){
    for(int u=alo[ku]; u<alo[ku+1]; u++){
      for(int v=u; v<alo[ku+1]; v++){
        int uv = mpos(u,v);
        double df = 0.0;
        for(int ku1=0; ku1<m->n; ku1++){
          for(int u1=alo[ku1]; u1<alo[ku1+1]; u1++){
            int v1 = u1;
            int i = mpos(u1,v1);
            double R1234 = R(u,v,u1,v1, mmmm, bo,m,qmd);
            df += R1234 * (D[i] + D[i]);
            for(int v1=u1+1; v1<alo[ku1+1]; v1++){
              int i = mpos(u1,v1);
              double R1234 = R(u,v,u1,v1, mmmm, bo,m,qmd);
              df += 4.0 * R1234 * D[i];
            }
          }
        }
        F[uv] += df;
      }
      for(int v=u; v<Mo; v++){
        int kv = bo->k[v];
        int uv = mpos(u,v);
        double ta = 0.0;
        for(int u1=alo[kv]; u1<alo[kv+1]; u1++){
          for(int v1=alo[ku]; v1<alo[ku+1]; v1++){
            int i = MPOSIF(u1,v1);
            double R1432 = R(u,v1,u1,v, mmmm, bo,m,qmd);
            ta -= R1432 * D[i];
          }
        }
        F[uv] += ta;
      }
    }
  }
  return;
}

void FX_r(int Mo, int Mv, double * F, double * X, double * FX){
  for(int u=0; u<Mo; u++){
    for(int v=0; v<Mo; v++){
      double d = 0.0;
      for(int a=0; a<Mv; a++){
        int au = a*Mo+u;
        int av = a*Mo+v;
        d += F[au] * X[av];
      }
      FX[u*Mo+v] = d;
    }
  }
  return;
}

void F2_eq6_r(int Mo, double * FX, double * s, double * F2){
  for(int v=0; v<Mo; v++){
    for(int u=0; u<=v; u++){
      int uv_sym = mpos(u,v);
      int uv = Mo*u+v;
      int vu = Mo*v+u;
      F2[uv_sym] = ( s[u]*FX[vu] + s[v]*FX[uv] ) * 0.5;
    }
  }
  return;
}

void X_eq7_r(double * F, double * X, basis * bo, basis * bv, qmdata * qmd){
  int Mo = bo->M;
  int Mv = bv->M;
  for(int a=0; a<Mv; a++){
    int qa = bv->Q[a];
    int la = bv->l[a];
    double ea = qmd->f2a[qa*qmd->nLv+la]; //eq41
    for(int v=0; v<Mo; v++){
      int    av = a*Mo+v;
      int    qv = bo->Q[v];
      int    lv = bo->l[v];
      double ev = qmd->fa[qv*qmd->nLo+lv]; //eq41
      double denom = 1.0/(ev-ea);
      X[av] = F[av]*denom;
    }
  }
  return;
}

void F_eq8_r(double * D, double * H, double * F,
    int * alo, int * alv, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd){
  int Mo = bo->M;
  int Mv = bv->M;
  for(int i=0; i<Mo*Mv; i++){
    F[i] = H[i];
  }
  for(int ka=0; ka<m->n; ka++){
    for(int a=alv[ka]; a<alv[ka+1]; a++){
      for(int v=alo[ka]; v<alo[ka+1]; v++){
        int av = a*Mo+v;
        double df = 0.0;
        for(int ku1=0; ku1<m->n; ku1++){
          if(ka==ku1){
            continue;
          }
          for(int u1=alo[ku1]; u1<alo[ku1+1]; u1++){
            int v1 = u1;
            int i = mpos(u1,v1);
            double R1234 = R2(a, v, u1, v1, pmmm, bo, bv, m, qmd);
            df += 2.0 * R1234 * D[i];
            for(int v1=u1+1; v1<alo[ku1+1]; v1++){
              int i = mpos(u1,v1);
              double R1234 = R2(a, v, u1, v1, pmmm, bo, bv, m, qmd);
              df += 4.0 * R1234 * D[i];
            }
          }
        }
        F[av] += df;
      }
      for(int v=0; v<Mo; v++){
        int kv = bo->k[v];
        if(kv==ka){
          continue;
        }
        int av = a*Mo+v;
        double ta = 0.0;
        for(int u1=alo[kv]; u1<alo[kv+1]; u1++){
          for(int v1=alo[ka]; v1<alo[ka+1]; v1++){
            int i = MPOSIF(u1,v1);
            double R1432 = R2(a, v1, u1, v,  pmmm, bo, bv, m, qmd);
            ta -= R1432 * D[i];
          }
        }
        F[av] += ta;
      }
    }
  }
  return;
}

