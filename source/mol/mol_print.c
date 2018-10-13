#include "mol.h"

void mol_print2(mol * m, FILE * f){
  fprintf(f, "Atomic Coordinates:\n");
  for(int i=0; i<m->n; i++){
    fprintf(f, "%d\t% .8lf\t% .8lf\t% .8lf\n",
        m->q[i], (m->r[3*i])*BA, (m->r[3*i+1])*BA, (m->r[3*i+2])*BA);
  }
  fprintf(f,"#\n\n");
  fflush (f);
}

void mol_print_m(mol * m, int bohr, const char s[], FILE * f){
  fprintf(f, "%s$molecule\n", s);
  if(m->z){
    fprintf(f, "%s charge=%d\n", s, m->z);
  }
  if(m->mult != 1){
    fprintf(f, "%s mult=%d\n", s, m->mult);
  }
  if(bohr){
    fprintf(f, "%s unit=b\n", s);
  }
  fprintf(f, "%s cartesian\n", s);

  for(int i=0; i<m->n; i++){
    if(bohr){
      fprintf(f, "%s%4d%15.8lf%15.8lf%15.8lf",
          s, m->q[i], m->r[3*i], m->r[3*i+1], m->r[3*i+2]);
    }
    else{
      fprintf(f, "%s%4d%15.8lf%15.8lf%15.8lf",
          s, m->q[i], m->r[3*i]*BA, m->r[3*i+1]*BA, m->r[3*i+2]*BA);
    }
    if(m->m[i] > 0){
      fprintf(f, "   mass=%lf", m->m[i]);
    }

    if (m->s[i][0] != 0){
      fprintf(f, "   type=%s", m->s[i]);
    }

    int a = m->l[i];
    int b = m->l[i+1];
    if(b-a > 0){
      fprintf(f, "   k=");
      for(int j=a; j<b; j++){
        fprintf(f, "%d", m->k[j]+1);
        if (m->b[j]!=1){
          fprintf(f, "(%d)", m->b[j]);
        }
        if(j!=b-1){
          fprintf(f, ",");
        }
      }
    }
    fprintf(f, "\n");
  }
  fprintf(f, "%s$end\n\n", s);
  fflush (f);
  return;
}

