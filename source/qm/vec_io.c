#include "matrix.h"
#include "qm.h"

static void vecarrange(unsigned int n, double * Ca, double * Ct, basis * bo){
  mx_transpcp(n, Ct, Ca);
  for(int i=0; i<n; i++){
    switch(bo->l[i]){
      case 0:
        veccp(n, Ca+n*i, Ct+n*i);
        break;
      case 1:
        switch(bo->m[i]){
          case -1:
            veccp(n, Ca+n*(i+1), Ct+n*i);
            break;
          case  0:
            veccp(n, Ca+n*(i+1), Ct+n*i);
            break;
          case  1:
            veccp(n, Ca+n*(i-2), Ct+n*i);
            break;
        }
        break;
      default: GOTOHELL;
    }
  }
  mx_transp(n, Ca);
  return;
}

int pvec_read(double * Va, double * Vb,
    double * Ca, double * Cb, const char s[], basis * bo){

  FILE * f;
  if( !(f = fopen(s, "r"))){
    return 0;
  }

  int32_t n;
  if( !fread(&n, sizeof(n), 1, f) || n!=bo->M ){
    fclose(f);
    return 0;
  }

  size_t vsize = sizeof(double)*n;
  size_t csize = sizeof(double)*n*n;
  if( !fread(Va, vsize, 1, f) || !fread(Ca, csize, 1, f) ||
      !fread(Vb, vsize, 1, f) || !fread(Cb, csize, 1, f) ){
    fclose(f);
    return 0;
  }

  double * Ct = malloc(csize);
  vecarrange(n, Ca, Ct, bo);
  vecarrange(n, Cb, Ct, bo);
  free(Ct);
  fclose(f);
  return 1;
}

int pvec_write(double * Va, double * Vb,
    double * Ca, double * Cb, const char s[], basis * bo){

  FILE * f;
  if( !(f = fopen(s, "w"))){
    return 0;
  }

  int32_t n = bo->M;
  size_t vsize = sizeof(double)*n;
  size_t csize = sizeof(double)*n*n;

  double * Ct = malloc(csize);
  vecarrange(n, Ca, Ct, bo);
  vecarrange(n, Ca, Ct, bo);
  vecarrange(n, Cb, Ct, bo);
  vecarrange(n, Cb, Ct, bo);
  free(Ct);

  if( !fwrite(&n, sizeof(n), 1, f) ||
      !fwrite(Va, vsize, 1, f) || !fwrite(Ca, csize, 1, f) ||
      !fwrite(Vb, vsize, 1, f) || !fwrite(Cb, csize, 1, f) ){
    fclose(f);
    return 0;
  }

  fclose(f);
  return 1;
}

