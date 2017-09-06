#include "qm.h"

basis * basis_fill(int M, mol * ml, int * Lq){
  size_t size = sizeof(basis) + sizeof(int)*5*M;
  basis   * b = malloc(size);
  b->M = M;
  b->Q = (int *)(b + 1);
  b->m = b->Q + M;
  b->l = b->m + M;
  b->k = b->l + M;
  b->n = b->k + M;
  int k,l,m,Q,L;
  int i = 0;
  for(k=0; k<ml->n; k++){
    Q = ml->q[k];
    L = Lq[Q];
    int n = 0;
    for(l=0; l<=L; l++){
      for(m=-l; m<=l; m++){
        b->Q[i] = Q;
        b->m[i] = m;
        b->l[i] = l;
        b->k[i] = k;
        b->n[i] = n;
        i++;
        n++;
      }
    }
  }
  return b;
}

int * basis_al(mol * m, int * Lq){
// (functions on atom #i) \in [ al[i] ; al[i+1] )
  int * al = malloc((m->n+1)*sizeof(int));
  al[0] = 0;
  for(int k=0; k<m->n; k++){
    int L = Lq[m->q[k]];
    al[k+1] = al[k] + (L+1)*(L+1);
  }
  return al;
}

int norb(mol * m, int * L){
  int M = 0;
  for(int i=0; i<m->n; i++){
    int Q  = m->q[i];
    int l  = L[Q];
    M += (l+1)*(l+1);
  }
  return M;
}

