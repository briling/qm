#include "qm.h"

/* Sab = Ca*S*Cb^T */
static double s2uhf_cont(int M, int Na, int Nb, double * Sab){
  double ds = Nb;
  for(int i=0; i<Na; i++){
    for(int j=0; j<Nb; j++){
      ds -= Sab[i*M+j]*Sab[i*M+j];
    }
  }
  return ds;
}

static double s2uhf_pure(int Na, int Nb){
  int Nu = Na-Nb;
  return Nu*0.5*(Nu*0.5+1.0);
}

void s2uhf(int M, int Na, int Nb, double * Sab, FILE * fo){
  double s0 = s2uhf_pure(Na, Nb);
  double ds = s2uhf_cont(M, Na, Nb, Sab);
  fprintf(fo, " <S^2> :  %.4lf / %.4lf\n\n", s0+ds, s0);
  return;
}

