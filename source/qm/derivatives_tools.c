#include "qm.h"
#include "matrix.h"
#include "derivatives.h"

void d2E2dF2_j(double * Da, double * Db, double * FaXa, double * FbXb,
               double * ja, double * jb, int * alo, basis * bo){
  /* j[u] == j[v], if lu==lv && ku==kv */
  /* !!! here i do not multiply by 2   */
  int Mo = bo->M;
  for(int u3=0; u3<Mo; u3++){
    int ku3 = bo->k[u3];
    int lu3 = bo->l[u3];
    int bra = alo[ku3]+lu3*lu3;
    int ket = alo[ku3]+(lu3+1)*(lu3+1);
    double da = 0.0;
    double db = 0.0;
    for(int v=bra; v<ket; v++){
      for(int u=0; u<Mo; u++){
        int uv   = u*Mo+v;
        int uv_s = MPOSIF(u,v);
        da += Da[uv_s] * FaXa[uv];
        db += Db[uv_s] * FbXb[uv];
      }
    }
    ja[u3] = da;
    jb[u3] = db;
  }
  return;
}

void d2E2dF2_j_r(double * D, double * FX,
                 double * j, int * alo, basis * bo){
  int Mo = bo->M;
  for(int u3=0; u3<Mo; u3++){
    int ku3 = bo->k[u3];
    int lu3 = bo->l[u3];
    int bra = alo[ku3]+lu3*lu3;
    int ket = alo[ku3]+(lu3+1)*(lu3+1);
    double d = 0.0;
    for(int v=bra; v<ket; v++){
      for(int u=0; u<Mo; u++){
        int uv   = u*Mo+v;
        int uv_s = MPOSIF(u,v);
        d += D[uv_s] * FX[uv];
      }
    }
    j[u3] = d;
  }
  return;
}

