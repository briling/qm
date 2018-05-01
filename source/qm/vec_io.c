#include "matrix.h"
#include "qm.h"
#include <stdint.h>

/* In sake of compatibility with Priroda,
 * we change the order of basis functions
 * when working with MO coefficient files.
 *
 * order in this program (m) | order in file (m1)
 *                           |
 * -- s shell :) --------------------------------
 *              0            |        0
 * -- p shell -----------------------------------
 *             -1            |        0
 *              0            |       +1
 *             +1            |       -1
 * -- d shell -----------------------------------
 *             -2            |        0
 *             -1            |       +1
 *              0            |       -1
 *             +1            |       +2
 *             +2            |       -2
 * ----------------------------------------------
 *
 * convert from file:  veccp( n, C + n*(i+(m1-m)), Ct + n*i );
 * convert to file:    veccp( n, Ct + n*i, C + n*(i+(m1-m)) );
 *
 * we can create an array of m1's:
 * n     : 0,  1,  2,  3,  4,  5  ...
 * M1[n] : 0, +1, -1, +2, -2, +3, ...
 * and take values m1 = M1[l+m],
 * but they can be calculated:
 * M1[n] = n%2 ? (n+1)/2 : -(n+1)/2
 */

#define M1(N) ((N)%2 ? ((N)+1)/2 : -((N)+1)/2)

static void vec_to_p(unsigned int n, double * Ct, double * C, basis * bo){
  mx_transp(n, C);
  for(unsigned int i=0; i<n; i++){
    int l  = bo->l[i];
    int m  = bo->m[i];
    int m1 = M1(m+l);
    veccp( n, Ct + n*i, C + n*(i+(m1-m)) );
  }
  mx_transp(n, Ct);
  mx_transp(n, C);
  return;
}

static void vec_from_p(unsigned int n, double * C, double * Ct, basis * bo){
  mx_transp(n, Ct);
  for(unsigned int i=0; i<n; i++){
    int l  = bo->l[i];
    int m  = bo->m[i];
    int m1 = M1(m+l);
    veccp( n, C + n*(i+(m1-m)), Ct + n*i );
  }
  mx_transp(n, C);
  return;
}

int pvec_read(double * Va, double * Vb,
    double * Ca, double * Cb, const char s[], basis * bo){

  FILE * f;
  if( !(f = fopen(s, "r"))){
    return 0;
  }

  uint32_t n;
  if( !fread(&n, sizeof(n), 1, f) || n!=bo->M ){
    fclose(f);
    return 0;
  }

  size_t vsize = sizeof(double)*n;
  size_t csize = sizeof(double)*n*n;
  double * Cta = malloc(csize);
  double * Ctb = malloc(csize);

  if( !fread(Va, vsize, 1, f) || !fread(Cta, csize, 1, f) ||
      !fread(Vb, vsize, 1, f) || !fread(Ctb, csize, 1, f) ){
    free(Cta);
    free(Ctb);
    fclose(f);
    return 0;
  }

  vec_from_p(n, Ca, Cta, bo);
  vec_from_p(n, Cb, Ctb, bo);

  free(Cta);
  free(Ctb);
  fclose(f);
  return 1;
}

int pvec_write(double * Va, double * Vb,
    double * Ca, double * Cb, const char s[], basis * bo){

  FILE * f;
  if( !(f = fopen(s, "w"))){
    return 0;
  }

  uint32_t n = bo->M;
  size_t vsize = sizeof(double)*n;
  size_t csize = sizeof(double)*n*n;
  double * Cta = malloc(csize);
  double * Ctb = malloc(csize);

  vec_to_p(n, Cta, Ca, bo);
  vec_to_p(n, Ctb, Cb, bo);

  if( !fwrite(&n, sizeof(n), 1, f) ||
      !fwrite(Va, vsize, 1, f) || !fwrite(Cta, csize, 1, f) ||
      !fwrite(Vb, vsize, 1, f) || !fwrite(Ctb, csize, 1, f) ){
    free(Cta);
    free(Ctb);
    fclose(f);
    return 0;
  }

  free(Cta);
  free(Ctb);
  fclose(f);
  return 1;
}

