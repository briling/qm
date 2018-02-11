#include "qm.h"
#include "matrix.h"
#include "eq.h"
#include "2el.h"
#include "tools.h"
#include "mytime.h"

#define VERSION    "v180211"
#define print_def  1
#define dDmax_def  1e-13
#define maxit_def  64

int main(int argc, char * argv[]){

  FILE * fm;
  FILE * fp;
  if (!((argc >= 1) && (fp = fopen(argv[1], "r"))))  GOTOHELL;
  if (!((argc >= 2) && (fm = fopen(argv[2], "r"))))  GOTOHELL;

  double dDmax = dDmax_def;
  int    maxit = maxit_def;
  int    print = print_def;
  char   vi[256] = {0};
  char   vo[256] = {0};
  FILE * fo = stdout;
  for(int i=3; i<argc; i++){
    if( sscanf (argv[i], "conv:%lf", &dDmax ) ) continue;
    if( sscanf (argv[i], "it:%d",    &maxit ) ) continue;
    if( sscanf (argv[i], "print:%d", &print ) ) continue;
    if( sscanf (argv[i], "read:%s",  &vi    ) ) continue;
    if( sscanf (argv[i], "write:%s", &vo    ) ) continue;
    if(! (fo = fopen(argv[i], "w"))){
      fo = stdout;
    }
  }
  fprintf(fo, "\n"VERSION"\n");
  fprintf(fo, "conv:%e\n", dDmax);
  fprintf(fo, "it:%d\n",   maxit);

  qmdata * qmd = qmdata_read(fp);
  fclose(fp);
#if 0
  qmdata_print(fo, qmd);
#endif
  mol * m = mol_read(fm);
  fclose(fm);
  if(!m){
    GOTOHELL;
  }

  int N = 0;
  for(int i=0; i<m->n; i++){
    N += nel(m->q[i], qmd);
  }
  N -= m->z;
  int Nu = m->mult-1;
  if( (N-Nu)%2 ) {
    fprintf(stderr, "\tN = %d, mult = %d !\n", N, m->mult);
    return 1;
  }
  int Nb = (N-Nu)/2;
  int Na = N-Nb;
  int Mo = norb(m, qmd->Lo);
  int Mv = norb(m, qmd->Lv);
  basis * bo  = basis_fill(Mo, m, qmd->Lo);
  basis * bv  = basis_fill(Mv, m, qmd->Lv);
  int   * alo = basis_al(m, qmd->Lo);
  int   * alv = basis_al(m, qmd->Lv);

  fprintf(fo, "\n");
  fprintf(fo, "  N   = %d,", N);
  fprintf(fo, "  Na  = %d,", Na);
  fprintf(fo, "  Nb  = %d\n", Nb);
  fprintf(fo, "  Mo  = %d,", Mo);
  fprintf(fo, "  Mv  = %d\n", Mv);
  fprintf(fo, "\n");

  double time_sec = myutime();

  // integrals -----------------------------------------------------------------

  double * f    = malloc(sizeof(double)*symsize(Mo));
  double * fmp  = malloc(sizeof(double)*Mo*Mv);
  double * H    = malloc(sizeof(double)*symsize(Mo));
  double * Hmp  = malloc(sizeof(double)*Mo*Mv);
  double * rij  = malloc(sizeof(double)*symsize(m->n));
  euler  * z    = malloc(sizeof(euler )*(m->n)*(m->n));
  distang(rij, z, m);
  double * mmmm = mmmm0_fill(alo, rij, z, bo, m, qmd);
  double * pmmm = pmmm_fill (alo, alv, z, bo, bv, m, qmd);
  f_eq25_mm(f,   z, alo,      bo,     m, qmd);
  f_eq25_mp(fmp, z, alo, alv, bo, bv, m, qmd);
  H_eq22_mm(f,   H,   alo,      mmmm, bo,     m, qmd);
  H_eq22_mp(fmp, Hmp, alo, alv, pmmm, bo, bv, m, qmd);
  mmmm6_add(alo, mmmm, rij, bo, m, qmd);
#if 0
  mmmm_check(mmmm, bo, m, qmd);
  pmmm_check(pmmm, bo, bv, m, qmd);
#endif
  free(z);
  free(rij);

  // init ----------------------------------------------------------------------

  double * Ca   = malloc(sizeof(double)*Mo*Mo);
  double * Cb   = malloc(sizeof(double)*Mo*Mo);
  double * Va   = malloc(sizeof(double)*Mo);
  double * Vb   = malloc(sizeof(double)*Mo);
  double * Da   = malloc(sizeof(double)*symsize(Mo));
  double * Db   = malloc(sizeof(double)*symsize(Mo));
  double * Dmp  = malloc(sizeof(double)*Mo*Mv);

  if(vi[0] && pvec_read(Va, Vb, Ca, Cb, vi, bo)){
    fprintf(fo, " read coefficients from '%s'\n\n", vi);
  }
  else{
    double * Fw = malloc(sizeof(double)*symsize(Mo));
    mx_id(Mo, Ca);
    veccp(symsize(Mo), Fw, f);
    jacobi(Fw, Ca, Va, Mo, 1e-15, 20, NULL);
    eigensort(Mo, Va, Ca);
    veccp(Mo*Mo, Cb, Ca);
    free(Fw);
  }

  // ---------------------------------------------------------------------------

  double E0 = E0_eq2(m, qmd);
  scf(Na, Nb, E0, Ca, Cb, Va, Vb, Da, Db, Dmp, maxit, dDmax, alo, alv, H, Hmp, mmmm, pmmm, bo, bv, m, qmd, fo);

  double dip[3] = {0.0};
  dipole(Da, Db, Dmp, dip, alo, alv, bo, bv, m, qmd);
  fprintf(fo, " dipole: %+10lf %+10lf %+10lf\n", dip[0], dip[1], dip[2]);
  spin2(Mo, Na, Nb, Ca, Cb, fo);

  if(print > 1){
    fprintf(fo, "  alpha:\n");
    mo_table(Na, Va, Ca, bo, fo);
    fprintf(fo, "  beta:\n");
    mo_table(Nb, Vb, Cb, bo, fo);
  }

  if(print > 2){
    population(Da, Db, alo, m, qmd, fo);
  }

  time_sec = myutime()-time_sec;
  fprintf(fo, "\nTIME: %.2lf sec  (%.2lf min)\n\n", time_sec, time_sec/60.0);

  if(vo[0] && pvec_write(Va, Vb, Ca, Cb, vo, bo)){
    fprintf(fo, " wrote coefficients to '%s'\n\n", vo);
  }

#if 0
  Deff_test(Na, Nb, Ca, Cb, Hmp, Dmp, alo, alv, pmmm, bo, bv, m, qmd);
#endif
#if 0
  Heff_test(Na, Nb, Ca, Cb, H, Hmp, alo, alv, mmmm, pmmm, bo, bv, m, qmd);
#endif

  // ---------------------------------------------------------------------------

  free(alo);
  free(alv);
  free(bo);
  free(bv);
  free(m);
  free(mmmm);
  free(pmmm);
  free(qmd);

  free(Ca);
  free(Cb);
  free(Va);
  free(Vb);
  free(Da);
  free(Db);
  free(Dmp);
  free(f);
  free(fmp);
  free(H);
  free(Hmp);

  fclose(fo);

  return 0;
}

