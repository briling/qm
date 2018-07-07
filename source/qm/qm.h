#include "qm_t.h"
#include "mol.h"

qmdata * qmdata_read  (FILE * f);
void     qmdata_print (FILE * f, qmdata * qmd);

basis * basis_fill(int M, mol * ml, int * Lq);
int   * basis_al(mol * m, int * Lq);
int     norb(mol * m, int * L);

int pvec_read (double * Va, double * Vb, double * Ca, double * Cb, const char s[], basis * bo);
int pvec_write(double * Va, double * Vb, double * Ca, double * Cb, const char s[], basis * bo);

void mo_table(int N, double * V, double * C, basis * b, FILE * fo);
void population(double * Da, double * Db, int * alo, mol * m, qmdata * qmd, FILE * fo);
double * Sab_fill(int M, double * Ca, double * Cb);

void s2uhf(int M, int Na, int Nb, double * Sab, FILE * fo);

double scf(int Na, int Nb, double E0,
    double * Ca, double * Cb, double * Va, double * Vb,
    double * Da, double * Db, double * Dmp,
    int maxit, double dDmax, int * alo, int * alv,
    double * H, double * Hmp, double * mmmm, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd, FILE * fo);
double scf_diis(int Na, int Nb, double E0,
    double * Ca, double * Cb, double * Va, double * Vb,
    double * Da, double * Db, double * Dmp,
    int maxit, int memit, double dDmax, int * alo, int * alv,
    double * H, double * Hmp, double * mmmm, double * pmmm,
    basis * bo, basis * bv, mol * m, qmdata * qmd, FILE * fo);

