#include "mol.h"
#include "vec3.h"
#include "matrix.h"
#include "qm.h"
#include "eq.h"
#include "2el.h"
#include "tools.h"
#include "mytime.h"

#define dDmax_def  1e-13
#define K_def      64

int main(int argc, char * argv[]){

  FILE * fm;
  FILE * fp;
  if (!((argc >= 1) && (fp = fopen(argv[1], "r"))))  GOTOHELL;
  if (!((argc >= 2) && (fm = fopen(argv[2], "r"))))  GOTOHELL;

  double dDmax = dDmax_def;
  int    K     = K_def;
  char   vi[256] = {0};
  char   vo[256] = {0};
  FILE * fo = stdout;
  for(int i=3; i<argc; i++){
    if( sscanf (argv[i], "tol:%lf",  &dDmax ) ) continue;
    if( sscanf (argv[i], "it:%d",    &K     ) ) continue;
    if( sscanf (argv[i], "read:%s",  &vi    ) ) continue;
    if( sscanf (argv[i], "write:%s", &vo    ) ) continue;
    if(! (fo = fopen(argv[i], "w"))){
      fo = stdout;
    }
  }

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
  double E0 = E0_eq2(m, qmd);

  // ---------------------------------------------------------------------------
  double * f    = malloc(sizeof(double)*symsize(Mo));
  double * H    = malloc(sizeof(double)*symsize(Mo));
  double * Fa   = malloc(sizeof(double)*symsize(Mo));
  double * Fb   = malloc(sizeof(double)*symsize(Mo));
  double * Da   = malloc(sizeof(double)*symsize(Mo));
  double * Db   = malloc(sizeof(double)*symsize(Mo));
  double * oldD = malloc(sizeof(double)*symsize(Mo));
  double * Fw   = malloc(sizeof(double)*symsize(Mo));
  double * Ca   = malloc(sizeof(double)*Mo*Mo);
  double * Cb   = malloc(sizeof(double)*Mo*Mo);
  double * Va   = malloc(sizeof(double)*Mo);
  double * Vb   = malloc(sizeof(double)*Mo);
  double * fmp  = malloc(sizeof(double)*Mo*Mv);
  double * Hmp  = malloc(sizeof(double)*Mo*Mv);
  double * Fmpa = malloc(sizeof(double)*Mo*Mv);
  double * Fmpb = malloc(sizeof(double)*Mo*Mv);
  double * Xa   = malloc(sizeof(double)*Mo*Mv);
  double * Xb   = malloc(sizeof(double)*Mo*Mv);
  double * FaXa = malloc(sizeof(double)*Mo*Mo);
  double * FbXb = malloc(sizeof(double)*Mo*Mo);
  double * sa   = malloc(sizeof(double)*Mo);
  double * sb   = malloc(sizeof(double)*Mo);
  double * F2a  = malloc(sizeof(double)*symsize(Mo));
  double * F2b  = malloc(sizeof(double)*symsize(Mo));
  double * FA   = malloc(sizeof(double)*symsize(Mo));
  double * FB   = malloc(sizeof(double)*symsize(Mo));
  double * rij  = malloc(sizeof(double)*symsize(m->n));
  euler  * z    = malloc(sizeof(euler )*(m->n)*(m->n));
  // ---------------------------------------------------------------------------

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

  if(vi[0] && pvec_read(Va, Vb, Ca, Cb, vi, bo)){
    fprintf(fo, " read coefficients from '%s'\n\n", vi);
  }
  else{
    mx_id(Mo, Ca);
    veccp(symsize(Mo), Fw, f);
    jacobi(Fw, Ca, Va, Mo, 1e-15, 20, NULL);
    eigensort(Mo, Va, Ca);
    veccp(Mo*Mo, Cb, Ca);
  }

  double E1 = 0.0;
  double E2 = 0.0;
  double oldE;
  vecset(symsize(Mo), oldD, 0.0);
  int k=0;
  while(k++<K){
    oldE = E1+E2;
    D_eq9 (Na, Mo, Ca, Da);
    D_eq9 (Nb, Mo, Cb, Db);
    F_eq4 (Da, Db, H, Fa, Fb, alo, mmmm, bo, m, qmd);
    F2_8_7_14_15_6(Da, Db, Hmp, Fmpa, Fmpb, Xa, Xb, FaXa, FbXb, sa, sb, F2a, F2b, alo, alv, pmmm, bo, bv, m, qmd);
    Heff(Da, Db, Xa, Xb, FaXa, FbXb, sa, sb, Fa, Fb, F2a, F2b, Fmpa, Fmpb, FA, FB, alo, alv, pmmm, bo, bv, m, qmd);
    mx_id(Mo, Ca);
    veccp(symsize(Mo), Fw, FA);
    jacobi(Fw, Ca, Va, Mo, 1e-15, 20, NULL);
    eigensort(Mo, Va, Ca);
    mx_id(Mo, Cb);
    veccp(symsize(Mo), Fw, FB);
    jacobi(Fw, Cb, Vb, Mo, 1e-15, 20, NULL);
    eigensort(Mo, Vb, Cb);
    E1 = E1_eq3(Mo, H, Da, Db, Fa, Fb);
    E2 = E2_eq5(Mo, Da, Db, F2a, F2b);

    double dD = 0.0;
    for(int i=0; i<symsize(Mo); i++){
      double ab = Da[i]+Db[i];
      double d  = oldD[i] - ab;
      dD += d*d;
      oldD[i] = ab;
    }

    double dE = E1+E2-oldE;
    fprintf(fo, " it %3d     E = % 17.10lf    dE = % 17.10lf    dD = % 5.2e\n", k, E0+E1+E2, k==1?0.0:dE, dD);
    if(dD < dDmax){
      fprintf(fo, "converged\n");
      break;
    }
  }

  fprintf(fo, "\n");
  fprintf(fo, " (E0   = %20.10lf)\n", E0);
  fprintf(fo, " (E0+1 = %20.10lf)\n", E0+E1);
  fprintf(fo, " (E2   = %20.10lf)\n", E2);
  fprintf(fo, "  E    = %20.10lf\n",  E0+E1+E2);
  spin2(Mo, Na, Nb, Ca, Cb, fo);

  fprintf(fo, "  alpha:\n");
  mo_table(Na, Va, Ca, bo, fo);
  fprintf(fo, "  beta:\n");
  mo_table(Nb, Vb, Cb, bo, fo);

  population(Da, Db, alo, m, qmd, fo);

  time_sec = myutime()-time_sec;
  fprintf(fo, "\nTIME: %.2lf sec  (%.2lf min)\n\n", time_sec, time_sec/60.0);

  if(vo[0] && pvec_write(Va, Vb, Ca, Cb, vo, bo)){
    fprintf(fo, " wrote coefficients to '%s'\n\n", vo);
  }

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
  free(Da);
  free(Db);
  free(F2a);
  free(F2b);
  free(FA);
  free(FB);
  free(Fa);
  free(FaXa);
  free(Fb);
  free(FbXb);
  free(Fmpa);
  free(Fmpb);
  free(Fw);
  free(H);
  free(Hmp);
  free(Va);
  free(Vb);
  free(Xa);
  free(Xb);
  free(f);
  free(fmp);
  free(oldD);
  free(rij);
  free(sa);
  free(sb);
  free(z);

  fclose(fo);

  return 0;
}

