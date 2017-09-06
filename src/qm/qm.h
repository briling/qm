#include "qm_t.h"
#include "mol.h"

qmdata * qmdata_read  (FILE * f);
void     qmdata_print (FILE * f, qmdata * qmd);

basis * basis_fill(int M, mol * ml, int * Lq);
int   * basis_al(mol * m, int * Lq);
int     norb(mol * m, int * L);

int pvec_read(double * Va, double * Vb, double * Ca, double * Cb, const char s[], basis * bo);

void mo_table(int N, double * V, double * C, basis * b, FILE * fo);
void spin2(int M, int Na, int Nb, double * Ca, double * Cb, FILE * fo);
void population(double * Da, double * Db, int * alo, mol * m, qmdata * qmd, FILE * fo);

