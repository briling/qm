
#include "matrix.h"

void mx_id(unsigned int n, double * a){
  unsigned int i,j;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      a[i*n+j] = (i==j ? 1.0 : 0.0);
    }
  }
}

void mx_print(unsigned int n, double * a, FILE   * f){
  unsigned int i,j;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      fprintf(f, "%25.15lf", a[i*n+j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
  fflush(f);
}

void mx_rect_print(unsigned int n, unsigned int m, double * a, FILE   * f){
  unsigned int i,j;
  for(i=0; i<n; i++){
    for(j=0; j<m; j++){
      fprintf(f, "%25.15lf", a[i*m+j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
  fflush(f);
}

void mx_sym_print(unsigned int n, double * a, FILE   * f){
  unsigned int i,N;
  N = (n*n-n)/2+n;
  for (i=0; i<N; i++){
    fprintf(f,"% .8e\t", a[i]);
    if (!((i+1) % 6)){
      fprintf(f,"\n");
    }
  }
  fprintf(f,"\n\n");
  fflush(f);
}

void mx_nosym_print(unsigned int n, double * a, FILE   * f){
  unsigned int i,j;
  for(i=0; i<n; i++){
    for(j=0; j<=i; j++){
      fprintf(f, "%25.15lf", a[mpos(j,i)]);
    }
    for(j=i+1; j<n; j++){
      fprintf(f, "%25.15lf", a[mpos(i,j)]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
  fflush(f);
}

double * mx_read(unsigned int n, FILE * f){
  unsigned int i,j;
  double * a;
  a = (double *)malloc(n*n*sizeof(double));
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      fscanf(f, "%lf", a+(i*n+j));
    }
  }
  return a;
}

double * mx_sym_read(unsigned int n, FILE * f){
  unsigned int i,N;
  double * b;
  N = (n*n-n)/2 + n;
  b = (double *)malloc(N*sizeof(double));
  for (i=0; i<N; i++){
    fscanf(f, "%lf", b+i);
  }
  return b;
}

void mx_transp(unsigned int n, double * a){
  unsigned int  i,j;
  double t;
  for(i=0; i<n; i++){
    for(j=i+1; j<n; j++){
      t = a[i*n+j];
      a[i*n+j] = a[j*n+i];
      a[j*n+i] = t;
    }
  }
  return;
}

void mx_transpcp (unsigned int   n, double * p, double * a){
  unsigned int  i,j;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      p[j*n+i] = a[i*n+j];
    }
  }
  return;
}

void mx_multtrmx(unsigned int n, double * p, double * a, double * b){
  /* A*Bt */
  unsigned int i,j,k;
  double t;
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      t=0.0;
      for(k=0; k<n; k++){
        t += a[i*n+k] * b[j*n+k];
      }
      p[i*n+j] = t;
    }
  }
  return;
}

