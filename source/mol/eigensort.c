#include "mol.h"
#include "vecn.h"

typedef struct{
  double   value ;
  double * vector;
} eigstr;

static int cmpev(const void *p1, const void *p2){
  double d;
  d = (*((eigstr *)p1)).value - (*((eigstr *)p2)).value;
  if (d > 0)
    return  1;
  else if (d < 0)
    return -1;
  else
    return  0;
}

void eigensort(int n, double * val, double * vec){
  double * tval = malloc((n*n+n)*sizeof(double));
  double * tvec = tval + n;
  if(!tval){
    abort();
  }
  veccp(n,   tval, val);
  veccp(n*n, tvec, vec);

  eigstr * eig = malloc(n*sizeof(eigstr));
  if(!eig){
    abort();
  }
  for(int i=0; i<n; i++){
    eig[i].value  = tval[i];
    eig[i].vector = tvec + n*i;
  }
  qsort(eig, n, sizeof(eigstr), cmpev);

  for(int i=0; i<n; i++){
    val[i] = eig[i].value;
    veccp(n, vec+i*n, eig[i].vector);
  }
  free(tval);
  free(eig);
  return;
}

