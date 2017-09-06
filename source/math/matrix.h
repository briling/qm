
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vecn.h"

#define symsize(M) (((M)*(M)+(M))/2)

#ifndef MPOS_IS
#define MPOS_IS
static inline unsigned int mpos(unsigned int i, unsigned int j){
/* A[i+j*(j+1)/2], i <= j, 0 <= j < N */
  return (i)+(((j)*((j)+1))>>1);
}
#define MPOSIF(i,j)   ((i)<=(j)? mpos((i),(j)):mpos((j),(i)))
#endif

void     mx_id          (unsigned int n, double * a);

void     mx_print       (unsigned int n, double * a, FILE * f);
void     mx_rect_print  (unsigned int n, unsigned int m, double * a, FILE * f);
void     mx_sym_print   (unsigned int n, double * a, FILE * f);
void     mx_nosym_print (unsigned int n, double * a, FILE * f);
double * mx_read      (unsigned int   n, FILE   * f);
double * mx_sym_read  (unsigned int   n, FILE   * f);

void     mx_transp    (unsigned int   n, double * a);
void     mx_transpcp  (unsigned int   n, double * p, double * a);

void     mx_multtrmx   (unsigned int n, double * p, double * a, double * b);

void     jacobi      (double * a, double * b, double * d, unsigned int n, double eps, unsigned int rot, FILE * f);



