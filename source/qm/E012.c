#include "eq.h"
#include "tools.h"
#include "matrix.h"
#include "vec3.h"

double E0_eq2(mol * m, qmdata * qmd){
  int i,j,qi,qj,ci,cj;
  double E0 = 0.0;
  for(i=0; i<m->n; i++){
    qi  = m->q[i];
    ci  = nel(qi, qmd);
    E0 += qmd->ea[qi];
    for(j=i+1; j<m->n; j++){
      qj  = m->q[j];
      cj  = nel(qj, qmd);
      E0 += ci*cj/sqrt(r3d2(m->r+i*3, m->r+j*3));
    }
  }
  return E0;
}

double E1_eq3(int Mo, double * H, double * Da, double * Db, double * Fa, double * Fb){
  double E1 = 0.0;
  for(int u=0;u<Mo;u++){
    for(int v=0;v<u;v++){
      int uv = mpos(v,u);
      E1 += (H[uv]+Fa[uv])*Da[uv] + (H[uv]+Fb[uv])*Db[uv] ;
    }
    int uv = mpos(u,u);
    E1 += 0.5*((H[uv]+Fa[uv])*Da[uv] + (H[uv]+Fb[uv])*Db[uv]) ;
  }
  return E1;
}

double E2_eq5(int Mo, double * Da, double * Db, double * F2a, double * F2b){
  double E2 = 0.0;
  for(int u=0;u<Mo;u++){
    for(int v=0;v<Mo;v++){
      int uv = MPOSIF(u,v);
      E2 += F2a[uv]*Da[uv];
      E2 += F2b[uv]*Db[uv];
    }
  }
  return E2;
}

