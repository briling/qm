#include "qm.h"
#include "matrix.h"
#include "tools.h"

void mo_table(int N, double * V, double * C, basis * b, FILE * fo){
  int i, j, M = b->M;
  fprintf(fo, "                 e:  ");
  for(i=0; i<M; i++){
    fprintf(fo, "% 15.10lf", V[i]);
  }
  fprintf(fo, "\n");
  fprintf(fo, "    N   atom  l  m         ");
  for(i=0; i<M; i++){
    fprintf(fo, "%5d (%c)      ", i+1, i<N?'o':'v');
  }
  fprintf(fo, "\n");
  for(i=0; i<M; i++){
    fprintf(fo, "%5d  %5d % 2d % 2d   ", i+1, b->k[i]+1, b->l[i], b->m[i]);
    for(j=0; j<M; j++){
      fprintf(fo, "% 15.10lf", C[j*M+i]);
    }
    fprintf(fo, "\n");
  }
  fprintf(fo, "\n");
  return;
}

void population(double * Da, double * Db, int * alo, mol * m, qmdata * qmd, FILE * fo){

  fprintf(fo, "charges:\n");
  for(int k=0; k<m->n; k++){
    double z = nel(m->q[k], qmd);
    for(int u=alo[k]; u<alo[k+1]; u++){
      int i = mpos(u,u);
      z -= Da[i] + Db[i];
    }
    fprintf(fo, "%4d(%2d)   %+8.6lf\n", k+1, m->q[k], z);
  }
  fprintf(fo, "\n");

  fprintf(fo, "bonds:\n");
  for(int ku=0; ku<m->n; ku++){
    for(int kv=ku+1; kv<m->n; kv++){
      double b = 0.0;
      for(int u=alo[ku]; u<alo[ku+1]; u++){
        for(int v=alo[kv]; v<alo[kv+1]; v++){
          int i = mpos(u,v);
          b += 2.0*(Da[i]*Da[i] + Db[i]*Db[i]);
        }
      }
      if(b>0.0625){
        fprintf(fo, "%4d(%2d) %4d(%2d)   %8.6lf\n", ku+1, m->q[ku], kv+1, m->q[kv], b);
      }
    }
  }
  fprintf(fo, "\n");
}

double * Sab_fill(int M, double * Ca, double * Cb){
  double * Sab = malloc(M*M*sizeof(double));
  mx_multtrmx (M, Sab, Ca, Cb);
  return Sab;
}

