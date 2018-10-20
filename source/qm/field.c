#include "qm.h"
#include "eq.h"
#include "matrix.h"
#include "vec3.h"
#include "tools.h"
#include "gradient.h"

double E0_ext(double field[3], mol * m, qmdata * qmd){
  double dip[3] = {0.0};
  for(int k=0; k<m->n; k++){
    int qk = m->q[k];
    int z  = nel(qk, qmd);
    r3adds(dip, m->r+k*3, z);
  }
  return r3dot(dip, field);
}

void E0_ext_grad(double field[3], double * g, double * Da, double * Db,
                 int * alo, mol * m, qmdata * qmd){
  for(int k=0; k<m->n; k++){
    int qk = m->q[k];
    double z = nel(qk, qmd);
    for(int v=alo[k]; v<alo[k+1]; v++){
      int vv = mpos(v,v);
      z -= Da[vv] + Db[vv];
    }
    r3adds(g+k*3, field, z);
  }
  return;
}

void H_ext(double field[3], double * F, int * alo, basis * bo, mol * m, qmdata * qmd){
  for(int k=0; k<m->n; k++){
    double rf = -r3dot(field, m->r+3*k);
    int qk = m->q[k];
    for(int v=alo[k]; v<alo[k+1]; v++){
      int lv = bo->l[v];
      int mv = bo->m[v];
      for(int u=alo[k]; u<v; u++){ // if u==v, duv==0
        int lu = bo->l[u];
        int mu = bo->m[u];
        int uv = mpos(u,v);
        double duv[3];
        dip_mm(duv, mu, mv, lu, lv, qk, qmd);
        F[uv] += r3dot(field, duv);
      }
      int vv = mpos(v,v);
      F[vv] += rf;
    }
  }
  return;
}

void Hmp_ext(double field[3], double * F,
             int * alo, int * alv, basis * bo, basis * bv,
             mol * m, qmdata * qmd){
  for(int k=0; k<m->n; k++){
    int Qk = m->q[k];
    for(int v=alo[k]; v<alo[k+1]; v++){
      int lv = bo->l[v];
      int mv = bo->m[v];
      for(int a=alv[k]; a<alv[k+1]; a++){
        int la = bv->l[a];
        int ma = bv->m[a];
        int av = a*bo->M+v;
        double dav[3];
        dip_pm(dav, ma, mv, la, lv, Qk, qmd);
        F[av] += r3dot(field, dav);
      }
    }
  }
  return;
}

