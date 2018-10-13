#include "qm.h"
#include "tools.h"
#include "gradient.h"
#include "vec3.h"

void E0_eq2_grad(double * g, mol * m, qmdata * qmd){
  for(int i=0; i<m->n; i++){
    int qi = m->q[i];
    int ci = nel(qi, qmd);
    for(int j=i+1; j<m->n; j++){
      int qj = m->q[j];
      int cj = nel(qj, qmd);
      double dr[3];
      r3diff(dr, m->r+j*3, m->r+i*3);
      double r21 = 1.0/r3dot(dr,dr);
      r3scal(dr, ci*cj*r21*sqrt(r21));
      r3add(g+i*3, dr);
      r3min(g+j*3, dr);
    }
  }
  return;
}

