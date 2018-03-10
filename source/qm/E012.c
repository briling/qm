#include "eq.h"
#include "tools.h"
#include "matrix.h"
#include "vec3.h"

double E0_eq2(mol * m, qmdata * qmd){
  double E0 = 0.0;
  for(int i=0; i<m->n; i++){
    int qi = m->q[i];
    int ci = nel(qi, qmd);
    E0 += qmd->ea[qi];
    for(int j=i+1; j<m->n; j++){
      int qj = m->q[j];
      int cj = nel(qj, qmd);
      E0 += ci*cj/sqrt(r3d2(m->r+i*3, m->r+j*3));
    }
  }
  return E0;
}

double E1_eq3(int Mo, double * H, double * Da, double * Db, double * Fa, double * Fb){
  double E1 = 0.0;
  for(int u=0; u<Mo; u++){
    int uu = mpos(u,u);
    E1 += 0.5*((H[uu]+Fa[uu])*Da[uu] + (H[uu]+Fb[uu])*Db[uu]) ;
    for(int v=u+1; v<Mo; v++){
      int uv = mpos(u,v);
      E1 += (H[uv]+Fa[uv])*Da[uv] + (H[uv]+Fb[uv])*Db[uv] ;
    }
  }
  return E1;
}

double E2_eq5(int Mo, double * Da, double * Db, double * F2a, double * F2b){
  double E2 = 0.0;
  for(int u=0; u<Mo; u++){
    int uu = mpos(u,u);
    E2 += 0.5*(F2a[uu]*Da[uu] + F2b[uu]*Db[uu]);
    for(int v=u+1; v<Mo; v++){
      int uv = mpos(u,v);
      E2 += F2a[uv]*Da[uv] + F2b[uv]*Db[uv];
    }
  }
  return 2.0*E2;
}

