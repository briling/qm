#include "common.h"

void g_print(int n, double * G, const char s[], FILE * f){
  for(int i=0; i<n; i++){
    fprintf(f, "%s %s% .8e % .8e % .8e\n", s, (i?"  ":"G="), G[3*i], G[3*i+1], G[3*i+2]);
  }
  fprintf(f, "\n");
  fflush(f);
}

